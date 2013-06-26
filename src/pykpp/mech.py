import os

from warnings import warn
from datetime import datetime, timedelta
from copy import deepcopy

from numpy import *
from scipy.constants import *
import scipy.integrate as itg

from plot import plot as _plot
from stdfuncs import *
import stdfuncs
from parse import _parsefile, _reactionstoic, _stoicreaction, _allspc

today = datetime.today()

def addtojac(rxni, reaction, jac, allspcs):
    """
    From parsed reaction, add a jacobian
    
    rxni - reaction list index
    reaction - reaction tokens
        dict(reactants = [(stoic, spc), ...],
             products = [(stoic, spc), ...])
    
    Future work should reduce the jacobian size by banding or LU decomp.
    """
    
    # Create a list of reactant species names
    rcts = [spc for stc, spc in reaction['reactants']]
    
    # Create a list of expressions that are multiplied by the 
    # rate constant (e.g., from k*y[0]*y[1] -> y[0]*y[1])
    spcorder = [('y[%d]**(%s)' % (allspcs.index(spc_), stc_)).replace('**(1.)', '').replace('**(+ 1.)', '') for stc_, spc_ in reaction['reactants']]
    
    # Iterate over all stoichiometries and spcs names for
    # reactions and products in reaction
    for role, stc, spc in [('r', stc_, spc_) for stc_, spc_ in reaction['reactants']] + [('p', stc_, spc_) for stc_, spc_ in reaction['products']]:
        # Species index in ordered species
        ispc = allspcs.index(spc)
        
        # For each reactant, calculate derivative
        # based on removing that reactant from the 
        # rate
        for rstc, rspc in reaction['reactants']:
            # Get index of reactant in ordered species
            irspc = allspcs.index(rspc)
            
            # Create a copy of the ordered species
            tmpspcorder = [i_ for i_ in spcorder]
            
            # Remove this reactant
            tmpspcorder.pop(rcts.index(rspc))
            
            # Create an expression to calculate derivative
            # of reaction
            expr = ' * '.join(['%s%s * rate_const[%d]' % ('-' if role == 'r' else '+', stc, rxni)] +  tmpspcorder)
            
            if expr[1:] in jac[ispc][irspc] and expr not in jac[ispc][irspc]:
                # The oposite signed expression is already in the jacobian
                # the sum of the two is 0, so this removes both expressiosn
                sign = expr[:1]
                if sign == '+':
                    sign = '-'
                elif sign == '-':
                    sign = '+'
                else:
                    raise ValueError('Huh?')
                jac[ispc][irspc] = jac[ispc][irspc].replace(sign + expr[1:], '')
            else:
                # Add expression to the jacobian
                jac[ispc][irspc] += expr
            

class Mech(object):
    """
    Mech object must be able to interpret kpp inputs
    and translate the results to create a reaction
    function (Mech.dy), a jacobian (Mech.ddy), and 
    run the model (Mech.run)
    """
    def __init__(self, path, verbose = False, keywords = ['hv', 'PROD', 'EMISSION'], timeunit = 'local', add_default_funcs = True):
        """
        path     - path to kpp inputs
        verbose  - add printing
        keywords - ignore certain keywords from reactants in
                   calculate reaction rates
        timeunit - 'local' or 'utc'
        add_default_funcs - if True, then add Update_THETA, Update_SUN, and Update_M if
                            THETA, SUN, or M appear in the rate constant expressions
        
        Special functions (e.g., Update_World) for physical environment
        can be added through PY_INIT in the definiton file e.g.,
        
        #MODEL      small_strato
        #INTEGRATOR lsoda
        #INLINE PY_INIT        
        TSTART = (12*3600)
        TEND = TSTART + (3*24*3600)
        DT = 0.25*3600  
        TEMP = 270
        updated_phys = -inf
        temp_times = array([TSTART, TEND])
        temps = array([270., 290.])
        press_times = array([TSTART, TEND])
        pressures = array([1013.25, 1050.])

        def NewUpdater(mech, world):
            global updated_phys
            global temp_times
            global temps
            global press_times
            global pressures
            t = world['t']

            if abs(t - updated_phys) > 1200.:
                world['TEMP'] = interp(t, temp_times, temps)
                world['P'] = interp(t, press_times, pressures)
                updated_phys = t
            Update_SUN(world)
            # Note that Update_RATE should be called or
            # rates will never be evaluated
            Update_RATE(mech, world) 
        
        add_world_updater(NewUpdater)
        #ENDINLINE
        
        see Mech.resetworld for special values that can be added
        to F90_INIT or INITVALUES to control the environment
        """
        self.add_default_funcs = add_default_funcs
        self.outputirr = False
        # Parse kpp inputs and store the results
        self.mechname, dummy = os.path.splitext(os.path.basename(path))

        self._parsed = _parsefile(path)
        
        # For rxn dict-like objects in EQUATIONS
        # remove reactants that are keywords
        for rxn in self._parsed['EQUATIONS']:
            remove = []
            for spci, (stc, spc) in enumerate(rxn['reactants']):
                if spc.strip() in keywords:
                    warn('Ignoring %s' % spc)
                    remove.append(spci)
            [rxn['reactants'].pop(spci) for spci in remove[::-1]]
            remove = []
            for spci, (stc, spc) in enumerate(rxn['products']):
                if spc.strip() in keywords:
                    warn('Ignoring %s' % spc)
                    remove.append(spci)
            [rxn['products'].pop(spci) for spci in remove[::-1]]
        
        # Find all spcs either explicitly defined or
        # referenced in reactions
        self.allspcs = _allspc(self._parsed)

        # Store a copy of the number of species
        nspcs = len(self.allspcs)
        
        # Print number of species and reactions 
        # and details of each if verbose = True
        print self.summary(verbose = verbose)
        
        # Create a world namespace to store
        # variables for calculating chemistry
        self.world = world = {}
        
        # Execute INIT code in the context of world
        if 'INIT' in self._parsed:
            exec(self._parsed['INIT'][-1], None, world)

        nocfactorkeys = tuple(world.keys()) + ('CFACTOR', 'TEMP', 'P', 'StartDate', 'StartJday', 'Latitude_Degrees', 'Latitude_Radians', 'Longitude_Degrees', 'Longitude_Radians')
        # Execute the INITVALUES code in the context
        # of the world
        exec(self._parsed['INITVALUES'][0], None, world)
        
        # Multiply all input values by CFACTOR
        cfactor = self.cfactor = world.get('CFACTOR', 1.)
        
        for k, v in world.iteritems():
            if k not in nocfactorkeys:
                world[k] = v * cfactor
        
        # Exec RCONST code in the context of stdfuncs
        if 'RCONST' in self._parsed:
            exec(self._parsed['RCONST'][-1], stdfuncs.__builtins__['globals'](), stdfuncs.__builtins__['locals']())
            
        self.Update_World = stdfuncs.Update_World

        
        if 'MONITOR' in self._parsed.keys():
            monitor = self._parsed['MONITOR'][0].replace(' ', '').split(';')
            monitor = [k_ for k_ in monitor if k_ != '']
            self.monitor = tuple([(self.allspcs.index(k_) if k_ in self.allspcs else None, k_) for k_ in monitor])
        else:
            self.monitor = tuple([ik for ik in enumerate(self.allspcs)])

        if 'LOOKAT' in self._parsed.keys():
            if self._parsed['LOOKAT'][0] == 'ALL':
                lookat = tuple(self.allspcs)
            else:
                lookat = self._parsed['LOOKAT'][0].replace(' ', '').split(';')
                lookat = [k_ for k_ in lookat if k_ != '']
                if 'ALL' in lookat:
                    lookat.pop(lookat.index('ALL'))
                    lookat = lookat + self.allspcs
            self.lookat = tuple(lookat)
        else:
            self.lookat = tuple(self.allspcs)
        
        # Store a temporary copy of the stoichiometry in each reaction
        dy_stoic = {}
        all_spec = world.get('ALL_SPEC', 0.)
        for spc in self.allspcs:
            if spc not in world:
                world[spc] = all_spec
            dy_stoic[spc] = _reactionstoic(spc, self._parsed['EQUATIONS'])
        
        # Add each reaction rate expression to the dy
        # for the species in the reacton
        dy_exp = ['' for i_ in range(nspcs)]
        for si, spc in enumerate(self.allspcs):
            for ri, stoic in dy_stoic[spc].iteritems():
                dy_exp[si] += (' + %f * rates[%d]' % (stoic, ri)).replace('+ -1.000000 * ', '-').replace('+ 1.000000 * ', '+')
        dy_exp = ['0' if dy_ == '' else dy_ for dy_ in dy_exp]
        self.dy_exp_str = 'array([' + ', '.join(dy_exp) + '])'
        self.dy_exp = compile(self.dy_exp_str, 'dy_exp', 'eval')
        # Create an empty copy of rate constant expressions (rate_const_exp)
        # and rate expressions (rate_exp; e.g., rate_const_exp * y[0] * y[1])
        # and rate derivatives with respect to a species (self.drate_exp)
        rate_const_exp = []
        rate_exp = []
        drate_exp = [['' for i_ in range(nspcs)] for j_ in range(nspcs)]
        for rxni, reaction in enumerate(self._parsed['EQUATIONS']):
            rate_const_exp.append(reaction['rate'])
            spcorder = [('y[%d]**(%s)' % (self.allspcs.index(spc_), stc_)).replace('**(1.)', '').replace('**(1.)', '') for stc_, spc_ in reaction['reactants']]
            rate_exp.append(' * '.join(['rate_const[%d]' % rxni] + spcorder))
            addtojac(rxni, reaction, drate_exp, self.allspcs)
        drate_exp = [['0' if col == '' else col for col in row] for row in drate_exp]
        self.drate_exp_str = 'array([' + ','.join(['[' + ','.join(row) + ']' for row in drate_exp]) + '])'
        self.drate_exp = compile(self.drate_exp_str, 'drate_exp', 'eval')
        self.rate_const_exp_str = 'array([' + ', '.join(rate_const_exp) + '])'
        self.rate_const_exp = compile(self.rate_const_exp_str, 'rate_const_exp', 'eval')
        self.rate_exp_str = 'array([' + ', '.join(rate_exp) + '])'
        self.rate_exp = compile(self.rate_exp_str, 'rate_exp', 'eval')
        self.parsed_world = world
        self.timeunit = timeunit
        
        usetuv = 'TUV_J' in self.rate_const_exp_str
        usetheta = 'THETA' in self.rate_const_exp_str
        usesun = 'SUN' in self.rate_const_exp_str
        self.usetheta = usetheta
        self.usesun = usesun
        if self.add_default_funcs:
            add_func = [(usetheta, Update_THETA),
                      (usesun, Update_SUN),
                      ('M' in self.rate_const_exp_str or 'N2' in self.rate_const_exp_str or 'O2' in self.rate_const_exp_str, Update_M)]
        
        for check, func in add_func:
            if check and not func in self.Update_World.updaters:
                self.Update_World.updaters.append(func)

        self.resetworld()

    def resetworld(self):
        """
        Resets the world to the parsed state and then adds default
        values that control the SUN and THETA (solar angle)
        
        
        Special values can be added to INITVALUES or F90_INIT
        to control solar properties
            StartDate = datetime.today()
            StartJday = int(StartDate.strftime('%j'))
            Latitude_Degrees = 0. [if Latitude_Radians is provided: degrees(Latitude_Radians)]
            Latitude_Radians = radians(Latitude_Degrees)
            SunRise = 4.5 or noon - degrees(arccos(-tan(LatRad) * tan(SolarDeclination))) / 15.
            SunSet  = 19.5 or noon - degrees(arccos(-tan(LatRad) * tan(SolarDeclination))) / 15.
        """
        world = self.world = deepcopy(self.parsed_world)
        if 'Latitude_Radians' not in world:
            LatDeg = world.pop('Latitude_Degrees', 45.)
            LatRad = world.setdefault('Latitude_Radians', radians(LatDeg))
        elif 'Latitude_Degrees' not in world:
            LatRad = world.setdefault('Latitude_Radians', radians(45.))
            
        if 'Longitude_Radians' not in world:
            LonDeg = world.pop('Longitude_Degrees', 0.)
            LonRad = world.setdefault('Longitude_Radians', radians(LonDeg))
        elif 'Longitude_Degrees' not in world:
            LonRad = world.setdefault('Longitude_Radians', radians(0.))
        

        world.setdefault('P', 101325.0)
        if (self.usetheta or self.usesun) and 'StartDate' not in world and 'StartJday' not in world:
            warn('Using SunRise and SunSet of 4.5 and 19.5 (approximately JulianDay 145 and Latitude 45 degrees N)')
            world.setdefault('SunRise', 4.5)
            world.setdefault('SunSet', 19.5)
            world.setdefault('StartJday', 145)
        else:
            if 'StartJday' not in world:
                StartDate = world.get('StartDate', 'datetime.today()')
                StartDate = world['StartDate'] = eval(StartDate)
                StartJday = world.setdefault('StartJday', int(StartDate.strftime('%j')))
            else:
                StartJday = world['StartJday']
                year = StartJday // 1000
                jday = StartJday % 1000
                if year < 1:
                    year = today.year
                    warn('Assuming current year %d' % year)
                StartDate = datetime(year, 1, 1) + timedelta(jday - 1)
            SolarDeclination = solar_declination(StartJday + ((world['TSTART'] + world['TEND']) / 2. / 3600.) // 24)
            half_day = degrees(arccos(-tan(LatRad) * tan(SolarDeclination))) / 15.
            if self.timeunit == 'local':
                solar_noon = solar_noon_local(LonDeg)
            elif self.timeunit == 'utc':
                solar_noon = solar_noon_utc(LonDeg)
            else:
                raise ValueError('timeunit must be either "local" or "utc"')
            world.setdefault('SunRise', solar_noon - half_day)
            world.setdefault('SunSet', solar_noon + half_day)
            if self.usetheta:
                world['SolarDeclination_Radians'] = SolarDeclination

    def summary(self, verbose = False):
        """
        return a string with species number, reaction number, and 
        string representations of each
        """
        if verbose:
            spcstr = ', '.join(self.allspcs)
            rxnstr = '\n' + '\n'.join(self.get_rxn_strs())
            initstr = '\nInitial Variables:\n' + self._parsed['INIT'][0]
            initvaluesstr = '\nInitial Concentrations:\n' + self._parsed['INITVALUES'][0]
        else:
            spcstr = ''
            rxnstr = ''
            initstr = ''
            initvaluesstr = ''
        return 'Species: %d %s\nReactions: %d%s%s%s' % (len(self.allspcs), spcstr, len(self.get_rxn_strs()), rxnstr, initstr, initvaluesstr)
    
    def get_rxn_strs(self, fmt = None):
        """
        Generate reaction string representations from parsed
        reaction  objects
        """
        if not fmt is None:
            fmt = fmt[1:]
        return [' + '.join(['*'.join([spc] if fmt is None else [format(eval(stc), fmt), spc]) for stc, spc in rxn['reactants']]) + ' = ' + ' + '.join(['*'.join([spc] if fmt is None else [format(eval(stc), fmt), spc]) for stc, spc in rxn['products']]) + ': ' + rxn['rate'] + ';' for rxn in self._parsed['EQUATIONS']]
    
    def print_rxns(self, fmt = None):
        """
        Print reaction stings joined by line returns
        """
        print '\n'.join(self.get_rxn_strs(fmt = fmt))
        
    def print_spcs(self, y, t):
        print '%.1f%%; {T=%.3E' % ((t - self.world['TSTART']) / (self.world['TEND'] - self.world['TSTART']) * 100, t),
        for spci, spc in self.monitor:
            if spci is None:
                try:
                    val = eval(spc, None, self.world)
                except:
                    val = nan
            else:
                val = y[spci] / self.cfactor
            print '%s=%.3E' % (spc, val),
        print '}'

    def output(self, outpath = None):
        lookat = self.lookat
        if not 't' in lookat:
            lookat = ('t',) + lookat
        if outpath is None:
            outpath = self.mechname + '.pykpp.dat'

        outfile = file(outpath, 'w')
        outfile.write(','.join([l_.ljust(8) for l_ in lookat]) + '\n')
        cfactor = self.cfactor
        for ti, time in enumerate(self.world['t']):
            outvals = []
            for k in lookat:
                v = self.world.get(k, nan)
                if hasattr(v, '__iter__'):
                    v = v[ti]
                if k in self.allspcs:
                    v = v / cfactor
                outvals.append(v)
            outfile.write(','.join(['%.8e' % v_ for v_ in outvals]) + '\n')
        outfile.seek(0, 0)
        return outfile

    def dy(self, y, t):
        """
        Function that returns the change in species
        with respect to time
        """
        time_since_print = t - getattr(self, 'monitor_time', 0) 
        self.world['t'] = t
        self.world['y'] = y
        self.Update_World(self, self.world)
        if time_since_print >= self.world['DT']:
            self.print_spcs(y, t)
            self.monitor_time = t

        rate_const = self.rate_const
        rates = eval(self.rate_exp)
            
        out = eval(self.dy_exp)
                    
        return out[:]

    def ddy(self, y, t):
        """
        Returns the derivative of the species change rate with
        respect to species
        """
        world = self.world
        world['t'] = t
        rate_const = self.rate_const
        tmp = {}
        tmp.update(locals())
        nspcs = len(self.allspcs)
        out = eval(self.drate_exp)
        #out = zeros((nspcs, nspcs))
        #for ri, row in enumerate(self.drate_exp):
        #    for ci, colexpr in enumerate(row):
        #        if colexpr != '':
        #            out[ri, ci] = eval(colexpr, None, tmp)
                
        return out

    def run(self, solver = None, tstart = None, tend = None, dt = None, jac = True, **solver_keywords):
        """
        Load solvers with Mech object function (Mech.dy), jacobian (Mech.ddy),
        and mechanism specific options
        """
        solver_trans = dict(kpp_lsode = 'lsoda')
        solvers = ('lsoda', 'vode', 'zvode', 'dopri5', 'dop853')
        if solver is None:
            parsed_solver = self._parsed['INTEGRATOR'][0]
            solver = solver_trans.get(parsed_solver, parsed_solver)
        if not solver in solvers:
            print 
            print 
            warn('Solver %s may not exist; if you get an error try one of (%s)' % (solver, ', '.join(solvers)))
            print 
            print 
        from time import time
        run_time0 = time()
        maxstepname = dict(lsoda = 'mxstep')
        default_solver_params = dict(atol = 1e-3, rtol = 1e-4, maxstep = 1000, hmax = self.world['DT'])
        for k, v in default_solver_params.iteritems():
            if k == 'maxstep':
                k = maxstepname.get(solver, 'max_step')
            solver_keywords.setdefault(k, v)
        if not hasattr(self, 'monitor_time'):
            self.monitor_time = -inf
        old_monitor_time = self.monitor_time
            
        if tstart is None: tstart = self.world.get('TSTART', 12*3600)
        if tend is None: tend = self.world.get('TEND', (12 + 24)*3600)
        if dt is None: dt = self.world.get('DT', 3600)
        
        y0 = array([eval(spc, None, self.world) for spc in self.allspcs])
        self.world['t'] = tstart
        self.Update_World(self, self.world)
        if solver == 'lsoda':
            ts = arange(tstart, tend + dt, dt)
            if solver_keywords.get('col_deriv', False):
                self.drate_exp = compile(self.drate_exp_str + '.T', 'drate_exp', 'eval')
                
            # With full details
            try:
                Y, infodict = itg.odeint(self.dy, y0, ts.copy(), Dfun = self.ddy if jac else None, mxords = 2, mxordn = 2, full_output = True, **solver_keywords)
                self.infodict = infodict
            except ValueError as e:
                raise ValueError(str(e) + '\n\n ------------------------------- \n If running again, you must reset the world (mech.resetworld())')

            # Without full details
            #Y = itg.odeint(self.dy, y0, ts, Dfun = self.ddy if jac else None, mxords = 2, mxordn = 2, **solver_keywords)
            
            # New method that forces steps to be independent
            # of each other; this makes the solver dumber.
            #Y = y0.copy()
            #infodicts = []
            #for start_end in ts.repeat(2, 0)[1:-1].reshape(-1, 2):
            #    Y_, infodict = itg.odeint(self.dy, y0, start_end, Dfun = self.ddy if jac else None, mxords = 2, mxordn = 2, full_output = True, **solver_keywords)
            #    Y = vstack([Y, Y_[-1]])
            #    infodicts.append(infodict)
            #self.infodict = infodicts
            
            if solver_keywords.get('col_deriv', False):
                self.drate_exp = compile(self.drate_exp_str, 'drate_exp', 'eval')
            
        else:
            solver_keywords.pop('hmax')
            # New method
            Y = y0
            ts = array([tstart])
            r = itg.ode(lambda t_, y_: self.dy(y_, t_), jac = lambda t_, y_: self.ddy(y_, t_) if jac else None)
            r.set_integrator(solver, method = 'bdf', **solver_keywords)
            r.set_initial_value(y = y0, t = tstart)
            while r.t < tend:
                nextt = r.t + dt
                while r.t < nextt:
                    try:
                        r.integrate(nextt)
                    except ValueError as e:
                        raise ValueError(str(e) + '\n\n ------------------------------- \n If running again, you must reset the world (mech.resetworld())')
                Y = vstack([Y, r.y])
                ts = append(ts, r.t)
            #    print r.t, r.y

        if getattr(self, 'monitor_time', -inf) != ts[-1]:
            self.print_spcs(Y[-1], ts[-1])
        run_time1 = time()
        self.world.update(dict(zip(self.allspcs, Y.T)))
        self.world['t'] = ts
        self.world['Y'] = Y
        self.monitor_time = old_monitor_time
        return run_time1 - run_time0
    
    def get_rates(self, **kwds):
        self.world.update(kwds)
        self.Update_World(self, self.world)
        self.world['y'] = y = array([eval(spc, None, self.world) for spc in self.allspcs])
        rate_const = self.world['rate_const'] = eval(self.rate_const_exp, None, self.world)
        rates = eval(self.rate_exp, None, self.world)
        
        return zip(self.get_rxn_strs(), rate_const, rates)
    
    def check_rates(self, **kwds):
        for lbl, a, val in self.get_rates(**kwds):
            print '%10.2e,%10.2e,%s' % (a, val, lbl)
        
    def print_world(self, out_keys = None, format = '%.8e', verbose = False):
        t = self.world['t']
        cfactor = self.world.get('CFACTOR', 1)
        if out_keys is None:
            out_keys = ['t'] + [k for i, k in self.monitor]
        if verbose:
            for ti in arange(t.size):
                print '%.1f%%.' % (float(ti) / (t.size - 1) * 100),
                for k in out_keys:
                    print ('%s=' + format) % (k, self.world[k][ti] / (1 if k == 't' else cfactor)),
                print
        else:
            print ','.join(out_keys)
            for ti in arange(t.size):
                print ','.join([format % self.world[k][ti] for k in out_keys])
    def plot(self, **kwds):
        return _plot(self, self.world, **kwds)