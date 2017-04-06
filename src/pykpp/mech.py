from __future__ import print_function, unicode_literals
import os
import sys

from warnings import warn

from datetime import datetime, timedelta
from copy import deepcopy

from numpy import *
import numpy as np
from scipy.constants import *
import scipy.integrate as itg

from .plot import plot as _plot
from pykpp import stdfuncs
from .stdfuncs import *
from . import stdfuncs
from .parse import _parsefile, _reactionstoic, _stoicreaction, _allspc, _prune_meta_species

today = datetime.today()

def addtojac(rxni, reaction, jac, allspcs):
    """
    From parsed reaction, add a jacobian element
    
    rxni - reaction list index
    reaction - reaction tokens
        dict(reactants = [(stoic, spc), ...],
             products = [(stoic, spc), ...])
    
    Future work should reduce the jacobian size by banding or LU decomp.
    
    jac[i,j] = d f[i] / d y[j]
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
    def add_world_updater(self, updater):
        self.updaters.append(updater)
    
    def add_func_updater(self, *args, **keywords):
        self.add_world_updater(func_updater(*args, **keywords))
    
    def __init__(self, path = None, mechname = None, keywords = ['hv', 'PROD', 'EMISSION'], incr = 300, monitor_incr = 3600, timeunit = 'local', add_default_funcs = True, doirr = False, verbose = False):
        """
        path              - path to kpp inputs
        mechname          - name for output
        keywords          - ignore certain keywords from reactants in
                            calculate reaction rates
        timeunit          - 'local' or 'utc'
        add_default_funcs - if True, then add Update_THETA, Update_SUN, and Update_M if
                            THETA, SUN, or M appear in the rate constant expressions.
        incr              - default functions will be called when simulated time has
                            progressed more than incr seconds
        monitor_incr      - If monitor_incr is not None, add a monitor function
                            updaters, so it has a special optional keyword
        doirr             - save reaction rates for further analysis
        verbose           - add printing
        
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
        self.verbose = verbose
        self.add_default_funcs = add_default_funcs
        self.outputirr = False
        self.updaters = []
        self.constraints = []
        self.path = path
        self.incr = incr
        self.monitor_incr = monitor_incr
        # Create a world namespace to store
        # variables for calculating chemistry
        self.world = world = {}
        
        # Add functions to world for easier parsing
        self.world['add_world_updater'] = self.add_world_updater
        self.world['add_func_updater'] = self.add_func_updater
        self.world['add_constraint'] = self.add_constraint
        
        # Set mechanism name
        self.set_mechname(mechname)
        
        if path is None:
            # A default mechanism has no reactions
            self._parsed = dict(EQUATIONS = [])
        else:
            self._parsed = _parsefile(path)
        
        # For rxn dict-like objects in EQUATIONS
        # remove reactants that are keywords
        for rxn in self._parsed['EQUATIONS']:
            _prune_meta_species(rxn, *keywords)
        
        # Find all spcs either explicitly defined or
        # referenced in reactions
        self.allspcs = _allspc(self._parsed)
        
        # Enable IRR - place holder
        self.doirr = doirr
        self.rate_history = []
        
        # Store a copy of the number of species
        nspcs = len(self.allspcs)
        
        # Print number of species and reactions 
        # and details of each if verbose = True
        print(self.summary(verbose = self.verbose))
        
        self.set_spcidx()

        # Execute INIT code in the context of world
        if 'INIT' in self._parsed:
            for initcode in self._parsed['INIT']:
                exec(initcode, None, world)

        # Do not apply CFACTOR to these key words even 
        # if they are defined in INITVALUES
        nocfactorkeys = tuple(world.keys()) + ('CFACTOR', 'TEMP', 'P', 'StartDate', 'StartJday', 'Latitude_Degrees', 'Latitude_Radians', 'Longitude_Degrees', 'Longitude_Radians')
        
        # Execute the INITVALUES code in the context
        # of the world
        if 'INITVALUES' in self._parsed:
            exec(self._parsed['INITVALUES'][0], None, world)
        
        # Multiply all input species by CFACTOR
        cfactor = self.start_cfactor = world.get('CFACTOR', 1.)
        for k, v in world.items():
            if k in self.allspcs and k not in nocfactorkeys:
                world[k] = v * cfactor
        
        # Add utility functions to the world
        if 'UTIL' in self._parsed:
            for utilcode in self._parsed['UTIL']:
                exec(utilcode, None, world)

        # Add RATES functions to the world
        if 'RATES' in self._parsed:
            for ratecode in self._parsed['RATES']:
                exec(ratecode, None, world)
                
        # Prepare to monitor expressions
        if 'MONITOR' in self._parsed.keys():
            monitor = self._parsed['MONITOR'][0].replace(' ', '').split(';')
            monitor = [k_ for k_ in monitor if k_ != '']
            self.monitor_expr = tuple([(self.allspcs.index(k_) if k_ in self.allspcs else None, k_) for k_ in monitor])
        else:
            self.monitor_expr = tuple([ik for ik in enumerate(self.allspcs)])
        
        # Prepare to archive expressions
        if 'LOOKAT' in self._parsed.keys():
            if self._parsed['LOOKAT'] == 'ALL':
                lookat = tuple(self.allspcs)
            else:
                lookat = self._parsed['LOOKAT'].replace(' ', '').split(';')
                lookat = [k_ for k_ in lookat if k_ != '']
                if 'ALL' in lookat:
                    lookat.pop(lookat.index('ALL'))
                    lookat = lookat + self.allspcs
            self.lookat = tuple(lookat)
        else:
            self.lookat = tuple(self.allspcs)
        
        # initialize non-defined species as ALL_SPC keyword
        self.set_all_spc()
        
        # Create dy expression and expression string (for review)
        self.set_dy_exp()
        
        # Create rate and drate expression and expression strings
        self.set_rate_exp()
        
        self.parsed_world = world
        del world['add_world_updater']
        del world['add_func_updater']
        del world['add_constraint']
        self.timeunit = timeunit
        
        # set over all time increment
        #self.set_incr(incr)
        
        # add default funcs is appropriate
        self.add_funcs()
        
        # Working on reorder and efficiency
        # so far lsoda chokes on packed format
        # so far lsoda chokes on baned format
        self.packed = False
        self.banded = False
        self.doreorder = False
        if self.doreorder or self.banded or self.packed:
            # Initialize clean world
            self.resetworld()
            neworder, self.ml, self.mu = self.reorder()
            self.resetworld()
            
            
        # Initialize clean world
        self.resetworld()

    def set_incr(self, incr = None):
        """
        Sets incr property
        
        incr - time increment (in rate unit); if None, select minimum time from updaters and constraints
        """
        if incr is None:
            self.incr = min([func.incr for func in self.updaters + self.constraints])
        else:
            self.incr = incr            
    
    def add_funcs(self):
        """
        Add default functions (e.g., Update_THETA, Update_SUN, Update_M, Monitor)
        if they are needed.
        """
        usetuv = 'TUV_J' in self.rate_const_exp_str
        usetheta = 'THETA' in self.rate_const_exp_str
        usesun = 'SUN' in self.rate_const_exp_str
        self.usetheta = usetheta
        self.usesun = usesun
        if not self.monitor_incr is None:
            self.add_world_updater(func_updater(Monitor, incr = self.monitor_incr, allowforce = False, verbose = self.verbose))
        if not self.add_default_funcs: return
        add_func = [(usetheta, Update_THETA),
                  (usesun, Update_SUN),
                  ('M' in self.rate_const_exp_str or 'N2' in self.rate_const_exp_str or 'O2' in self.rate_const_exp_str, Update_M)]
        for check, func in add_func:
            if check and not func in self.updaters:
                self.add_world_updater(func_updater(func, incr = self.incr, verbose = self.verbose))

    def set_all_spc(self):
        """
        Any undefined species are set to ALL_SPEC keyword
        or 0 if ALL_SPEC is undefined
        """
        world = self.world
        all_spec = world.get('ALL_SPEC', 0.)
        for spc in self.allspcs:
            if spc not in world:
                world[spc] = all_spec
    
    def set_rate_exp(self):
        """
        This function serves three major goals:
        
        1) Define the rate_const_exp and rate_const_exp_str
           rate_const_exp - is evaluated to calculate all rate coefficients
        2) Define the rate_exp and rate_exp_str
           rate_exp - is evaluated to calculate all reaction rates
        3) Define drate_exp, fill_drate_exp (and _str's)
           drate_exp - can be evaluated to calculate the derivative of the 
                       rate_exp with respect to y. This is currently 
                       deprecated in favor of fill_drate_exp
           fill_drate_exp - can be executed to fill the drate_per_dspc array,
                       which must exist. This fills the same roll as drate_exp,
                       but is much faster because it only allocates the memory once.
        """
        # Store a copy of the number of species
        nspcs = len(self.allspcs)
        nrxns = len(self._parsed['EQUATIONS'])
        world = self.world
        # Create an empty copy of rate constant expressions (rate_const_exp)
        # and rate expressions (rate_exp; e.g., rate_const_exp * y[0] * y[1])
        # and rate derivatives with respect to a species (self.drate_exp)
        rate_const_exp = []
        rate_exp = []
        drate_exp = [['' for i_ in range(nspcs)] for j_ in range(nspcs)]
        self.fill_rate_const_exp_str = '\n'.join(self._parsed.get('RCONST', '')) + '\n'
        for rxni, reaction in enumerate(self._parsed['EQUATIONS']):
            rate_const_exp.append(reaction['rate'])
            spcorder = [('y[%d]**(%s)' % (self.allspcs.index(spc_), stc_)).replace('**(1.)', '').replace('**(1.)', '') for stc_, spc_ in reaction['reactants']]
            rate_exp.append(' * '.join(['rate_const[%d]' % rxni] + spcorder))
            addtojac(rxni, reaction, drate_exp, self.allspcs)
        
        drate_exp = [['0' if col == '' else col for col in row] for row in drate_exp]
        self.drate_exp_str = 'array([' + ','.join(['[' + ','.join(row) + ']' for row in drate_exp]) + '])'
        self.drate_per_dspc = np.zeros((nspcs, nspcs), dtype = 'd')
        self.fill_drate_exp_str = '\n'.join(['\n'.join(['drate_per_dspc[%i,%i] = %s' % (j, i, v) for i, v in enumerate(row) if v != '0']) for j, row in enumerate(drate_exp)])
        self.fill_drate_exp = compile(self.fill_drate_exp_str, 'drate_exp', 'exec')
        self.drate_exp = compile(self.drate_exp_str, 'drate_exp', 'eval')
        self.rate_const_exp_str = 'array([' + ',\n '.join(rate_const_exp) + '])'
        self.rate_const_exp = compile(self.rate_const_exp_str, 'rate_const_exp', 'eval')
        self.fill_rate_const_exp = compile(self.fill_rate_const_exp_str + 'rate_const[:] = ' + self.rate_const_exp_str, 'fill_rate_const_exp', 'exec')
        self.rate_const = np.zeros((nrxns,), dtype = 'd')
        self.rate_exp_str = 'array([' + ',\n '.join(rate_exp) + '])'
        self.fill_rate_exp_str = '\n'.join(['rates[%d] = %s' % ri_rassign for ri_rassign in enumerate(rate_exp)])
        self.fill_rate_exp = compile(self.fill_rate_exp_str, 'fill_rate_exp', 'exec')
        self.rate_exp = compile(self.rate_exp_str, 'rate_exp', 'eval')

    def set_dy_exp(self):
        # Store a copy of the number of species
        nspcs = len(self.allspcs)
        world = self.world
        # Store a temporary copy of the stoichiometry in each reaction
        dy_stoic = {}
        for spc in self.allspcs:
            dy_stoic[spc] = _reactionstoic(spc, self._parsed['EQUATIONS'])
        
        # Add each reaction rate expression to the dy
        # for the species in the reacton
        dy_exp = ['' for i_ in range(nspcs)]
        for si, spc in enumerate(self.allspcs):
            for ri, stoic in dy_stoic[spc].items():
                dy_exp[si] += (' + %f * rates[%d]' % (stoic, ri)).replace('+ -1.000000 * ', '-').replace('+ 1.000000 * ', '+')
        dy_exp = ['0' if dy_ == '' else dy_ for dy_ in dy_exp]
        self.dy_exp_str = 'array([' + ',\n '.join(dy_exp) + '])'
        self.fill_dy_exp_str = '\n'.join(['dy[%d] = %s' % yi_assign for yi_assign in enumerate(dy_exp)])
        self.dy_exp = compile(self.dy_exp_str, 'dy_exp', 'eval')
        self.fill_dy_exp = compile(self.fill_dy_exp_str, 'fill_dy_exp', 'exec')
        
    def set_spcidx(self):
        # Create references to spc in conc array
        for si, spcn in enumerate(self.allspcs):
            self.world['ind_'+spcn] = si
        
        self.world['NSPCS'] = len(self.allspcs)


    def set_mechname(self, mechname):
        if mechname is None:
            # Heuristically determine the mechanism name
            if isinstance(self.path, (str, unicode)):
                # Parse kpp inputs and store the results
                self.mechname, dummy = os.path.splitext(os.path.basename(self.path))
            else:
                self.mechname = 'unknown'
        else:
            self.mechname = mechname

    def Update_World(self, world = None, forceupdate = False):
        """
        world - dictionary defining current state
        forceupdate - boolean that forces all functions to update
        Calls 
            + user added functions (added using add_world_updater)
            update_func_world(world)
            Update_RATE(mech, world)
    
        *Note: Update_SUN, Update_THETA, and Update_M can be added by
               the mechanism heuristically.
    
        see these functions for more details
        """
        if world is None:
            world = self.world
        time2update = False
        for func in self.updaters + self.constraints:
            updatedf = func(self, world, force = forceupdate)
            time2update = time2update or updatedf
            
        if time2update or (world['t'] - self.last_updated) > self.incr:
            self.last_updated = world['t']
            if self.verbose: print(self.last_updated, 'Updating world')
            update_func_world(self, world)
            self.update_rate_const()
    
    def update_rate_const(self):
        self.world['rate_const'] = self.rate_const
        exec(self.fill_rate_const_exp, None, self.world)
            
    def resetworld(self):
        """
        Resets the world to the parsed state and then adds default
        values that control the SUN and THETA (solar angle)
        
        
        Special values can be added to INITVALUES or F90_INIT
        to control solar properties
            StartDate = datetime.today()
            StartJday = int(StartDate.strftime('%j'))
            Latitude_Degrees = 45. [if Latitude_Radians is provided: degrees(Latitude_Radians)]
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
            warn('WARN: Using SunRise and SunSet of 4.5 and 19.5 (approximately JulianDay 145 and Latitude 45 degrees N)')
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
            if 'TSTART' not in world:
                warn('Assuming solar noon TSTART')
                world['TSTART'] = 12. * 3600.
            if 'TEND' not in world:
                world['TEND'] = 12. * 3600.
                warn('Assuming solar noon TEND')
                
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
        self.world.setdefault('DT', self.incr)
        self.world.setdefault('MONITOR_DT', self.world['DT'])
        for obj in self.updaters + self.constraints:
            obj.reset()
        self.world['history'] = {}
        self.world['allspcs'] = self.allspcs
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
        print('\n'.join(self.get_rxn_strs(fmt = fmt)))
        
    def print_monitor(self, y, t):
        """
        
        """
        out = '{t:%.0f' % (t,)
        for spci, spc in self.monitor_expr:
            if spci is None:
                try:
                    val = eval(spc, None, self.world)
                except:
                    val = nan
            else:
                val = y[spci] / self.start_cfactor
            out += ',%s:%.2G' % (spc, val)
        print (out + '}')

    def output(self, outpath = None):
        if outpath is None:
            outpath = self.mechname + '.pykpp.tsv'
        outfile = open(outpath, 'w')
        self._archive.seek(0, 0)
        outfile.write(self._archive.read())
        outfile.close()

    def archive(self):
        lookat = self.lookat
        if not 't' in lookat:
            lookat = ('t',) + lookat
        if not 'P' in lookat and 'P' in self.world:
            lookat = ('P',) + lookat
        if not 'TEMP' in lookat and 'TEMP' in self.world:
            lookat = ('TEMP',) + lookat
        if not 'CFACTOR' in lookat and 'CFACTOR' in self.world:
            lookat = ('CFACTOR',) + lookat
        if not hasattr(self, '_archive'):
            from io import StringIO
            self._archive = StringIO()
            self._archive.write('\t'.join([l_ for l_ in lookat]) + '\n')
        
        outfile = self._archive
        cfactor = self.start_cfactor
        outvals = []
        for k in lookat:
            try:
                v = eval(k, None, self.world)
            except:
                v = self.world.get(k, nan)
            if k in self.allspcs:
                v = v / cfactor
            outvals.append(v)
        outfile.write(u'\t'.join(['%.8e' % v_ for v_ in outvals]) + '\n')
        return outfile

    #def monitor(self, y, t):
    #    """
    #    Print state at t
    #    """
    #    time_since_print = t - getattr(self, 'monitor_time', 0) 
    #    if time_since_print >= self.world['MONITOR_DT']:
    #        self.print_monitor(y, t)
    #        self.monitor_time = t
    #        if self.doirr:
    #            self.rate_history.append(rates)
        
    def dy(self, y, t):
        """
        Function that returns the change in species
        with respect to time
        """
        self.world['t'] = t
        self.update_y_from_y(y)
        self.Update_World(self.world)
        rate_const = self.rate_const
        rates = eval(self.rate_exp)
        #rates = self.arates
        #exec(self.fill_rate_exp)
        #dy = out = self.ady
        #exec(self.fill_dy_exp)
        out = eval(self.dy_exp)
        return out[:]
    
    def reorder(self):
        self.Update_World()
        nspcs = len(self.allspcs)
        y = self.update_y_from_world()
        y[:]*=0
        y[:]+=1

        t = 0
        drate_per_dspc = self.drate_per_dspc
        rate_const = self.rate_const
        rate_const[:] *= 0
        rate_const[:] += 1
        exec(self.fill_drate_exp)
        inmatrix = self.drate_per_dspc
        nonzero = inmatrix[:] != 0
        density = nonzero.sum(0) + nonzero.sum(1)
        densityasc = density.argsort()
        neworder = np.concatenate([densityasc[::-1][::2][::-1], densityasc[::-1][1::2]])
        from_left = nonzero.argmax(1)
        from_right = nonzero[:, ::-1].argmax(1)
        from_bottom = nonzero.argmax(0)
        from_top = nonzero[::-1].argmax(0)
        neworder = (from_bottom - from_top).argsort()
        neworder = (from_left - from_right).argsort()
        testmatrix = inmatrix[neworder][:, neworder].copy()
        """Returns ml and mu, the lower and upper band sizes of a."""
        a = testmatrix
        nrows, ncols = a.shape
        ml = 0
        for k in range(-nrows+1, 0):
            if np.diag(a, k).any():
                ml = -k
                break
        mu = 0
        for k in range(nrows-1, 0, -1):
            if np.diag(a, k).any():
                mu = k
                break
        subd = ml; superd = mu
        self.allspcs = [self.allspcs[i] for i in neworder]
        del self.world['y']
        self.set_spcidx()
        self.set_dy_exp()
        self.set_rate_exp()
        self.ml = subd
        self.mu = superd
        return neworder, subd, superd
        
        drate_per_dspc = self.drate_per_dspc
        rate_const = self.rate_const
        drate_per_dspc = self.drate_per_dspc
        exec(self.fill_drate_exp)
        for d in range(subd, nspcs):
            assert((np.diagonal(drate_per_dspc, -d) == 0).all())
            
        for d in range(superd, nspcs):
            assert((np.diagonal(drate_per_dspc, d) == 0).all())
            
        return neworder, subd, superd    
    
    def ddy(self, y, t):
        """
        Returns the derivative of the species change rate with
        respect to species
        """
        self.world['t'] = t
        self.update_y_from_y(y)
        rate_const = self.rate_const
        #out = eval(self.drate_exp)
        out = drate_per_dspc = self.drate_per_dspc
        exec(self.fill_drate_exp)
        if self.banded:
            diagonals = []
            for di, d in enumerate(range(-self.ml, self.mu+1)):
                if d < 0:
                    before = [0] * -d
                    after = []
                elif d > 0:
                    after = [0] * d
                    before = []
                else:
                    after = []
                    before = []
                diagonals.append(np.append(np.append(before, np.diagonal(out, d)), after)[None])
            out = np.concatenate(diagonals, axis = 0)
            print(1)
        elif self.packed:
            nr = 2*self.ml+self.mu+1
            jac_packed = np.zeros((nr, len(y)), dtype = out.dtype)
            for (i, j), v in np.ndenumerate(out):
                myr = i-j+self.mu
                if j >= (i-self.ml) and j <= (i+self.mu):
                    jac_packed[myr, j] = v
                else:
                    assert(v == 0)
            out = jac_packed
        return out
    
    def get_y(self, world = None):
        """
        Return a vector of species values in order of self.allspcs
        Optional:
            world - special dictionary to create vector (defaults to self.world)
        """
        if world is None:
            world = self.world
        return array([eval(spc, None, world) for spc in self.allspcs])
    
    def update_y_from_y(self, y):
        """
        Update y if it exists or create it, if it does not
        """
        if 'y' in self.world:
            self.world['y'][:] = y
        else:
            self.world['y'] = y.copy()
        
        return self.world['y']
    
    def update_y_from_world(self):
        """
        Create y vector in world with values for each species in allspcs
        """
        self.update_y_from_y(self.get_y())
        return self.world['y']
    
    def update_world_from_y(self, y = None):
        """
        Update world state from y
        """
        world = self.world
        if y is None:
            y = world['y']
        else:
            self.update_y_from_y(y)
        
        for si, spc in enumerate(self.allspcs):
            world[spc] = y[si]
        
            
    def set_spcs(self, **spcs):
        """
        Set each spcs to values in world and world['y']
        """
        self.world.update(kwds)
        return self.update_y_from_world()
        
    def integrate(self, t0, t1, y0, solver = None, jac = True, verbose = False, **solver_keywords):
        """
        Arguments:
            t0 - starting time in (unit of kinetic rates)
            t1 - ending time in (unit of kinetic rates)
            y0 - vector of concentrations (unit of kinetic rates)
            solver - (optional; default None) string indicating solver approach
                * None if not provided, defaults to integrator set in mechanism defintion
                * 'odeint' use odeint with lsoda from ODEPACK (scipy.integrate.odeint)
                * 'lsoda', 'vode', 'zvode', 'dopri5', 'dopri853' use named solver with ode object (see scipy.integrate.ode)
            jac - (optional; default True) Use jacobian to speed up solution (default: True)
            verbose - (optional; default False)add additional printing
            solver_keywords - (optional) solver keywords are specific to each solver; see scipy.integrate.ode for more details on keywords
                * with odeint, defaults to dict(atol = 1e-3, rtol = 1e-4, maxstep = 1000, hmax = self.world['DT'], mxords = 2, mxordn = 2)
                * with ode, no defaults are set

        Returns:
            ts - array of times (start, end)
            Y - array of states with dimensions (time, spc)
        """
        if self.banded or self.packed:
            solver_keywords['uband'] = self.mu
            solver_keywords['lband'] = self.ml
        
        solver_trans = dict(kpp_lsode = 'odeint')
        solvers = ('lsoda', 'vode', 'zvode', 'dopri5', 'dop853', 'odeint')
        if solver is None:
            parsed_solver = self._parsed['INTEGRATOR'][0]
            solver = solver_trans.get(parsed_solver, parsed_solver)
        if not solver in solvers:
            print() 
            print() 
            warn('Solver %s may not exist; if you get an error try one of (%s)' % (solver, ', '.join(solvers)))
            print()
            print()
        elif self.verbose or verbose:
            print('Solver:', solver)
        if solver == 'odeint':
            maxstepname = dict(odeint = 'mxstep')
            default_solver_params = dict(atol = 1e-3, rtol = 1e-4, maxstep = 1000, hmax = self.world['DT'], mxords = 2, mxordn = 2)
            for k, v in default_solver_params.items():
                if k == 'maxstep':
                    k = maxstepname.get(solver, 'max_step')
                solver_keywords.setdefault(k, v)
            if solver_keywords.get('col_deriv', False):
                self.drate_exp = compile(self.drate_exp_str + '.T', 'drate_exp', 'eval')
                
            # With full details
            if len(self.constraints) > 0:
                warn("Constraints affect rates, but are not seen in outputs with odeint")
            if verbose:
                print(solver, solver_keywords)
            try:
                ts = array([t0, t1])
                Y = y0
                dy = self.dy
                ddy = self.ddy
                Ystep, infodict = itg.odeint(dy, y0, ts, Dfun = ddy if jac else None, full_output = True, **solver_keywords)
                y0 = Ystep[-1]
                Y = vstack([Y, y0])
                self.infodict = infodict
                self.update_world_from_y(Ystep[-1])
                if self.infodict['message'] != "Integration successful.":
                    print(self.infodict['message'])
            except ValueError as e:
                raise ValueError(str(e) + '\n\n ------------------------------- \n If running again, you must reset the world (mech.resetworld())')

            if solver_keywords.get('col_deriv', False):
                self.drate_exp = compile(self.drate_exp_str, 'drate_exp', 'eval')
            
        else:
            # New method
            Y = y0
            ts = array([t0])
            def ody(t, y):
                return self.dy(y, t)
            
            def oddy(t, y):
                return self.ddy(y, t) 
            
            if self.verbose:
                print(solver, solver_keywords)
            
            if jac:
                r = itg.ode(ody, jac = oddy)
            else:
                r = itg.ode(ody)
            
            r.set_initial_value(y = y0, t = t0)
            #def solout(t, y):
            #    self.world['t'] = t
            #    self.update_world_from_y(y)
            #    self.update_y_from_y(y)
            #    self.Update_World()
            #r.set_solout(solout)
            nextt = t1
            r.stiff = 1
            while r.t < nextt:
                r.set_integrator(solver, **solver_keywords)
                try:
                    r.integrate(nextt)
                except ValueError as e:
                    raise ValueError(str(e) + '\n\n ------------------------------- \n If running again, you must reset the world (mech.resetworld())')
                if r.t < nextt:
                    self.Update_World(self.world, forceupdate = True)
                    self.update_world_from_y(r.y.copy())
                    if len(self.constraints) > 0:
                        self.update_y_from_world()
                        self.Update_World(self.world, forceupdate = True)
                    r.set_initial_value(self.world['y'], r.t)
                if (self.verbose or verbose) and not r.successful():
                    warn('Partial step from %.0f to %.0f instead of %.0f' % (t0, r.t, nextt))

            Y = vstack([Y, r.y])
            ts = append(ts, r.t)
            del r
        
        return ts, Y
        
    def run(self, solver = None, tstart = None, tend = None, dt = None, jac = True, **solver_keywords):
        """
        Load solvers with Mech object function (Mech.dy), jacobian (Mech.ddy),
        and mechanism specific options
        """
        cfactor = self.world['CFACTOR']
        from time import time
        run_time0 = time()
        if 'atol' in self.world:
            solver_keywords.setdefault('atol', self.world['atol'])
        if 'rtol' in self.world:
            solver_keywords.setdefault('rtol', self.world['rtol'])
        
        if tstart is None: tstart = self.world.get('TSTART', 12*3600)
        if tend is None: tend = self.world.get('TEND', (12 + 24)*3600)
        if dt is None: dt = self.world.get('DT', 3600)
        
        t = self.world['t'] = tstart
        self.last_updated = t
        self.Update_World(self.world, forceupdate = True)
        self.archive()
        while t < tend:
            # Get y vector from world state
            y0 = self.get_y()
            ts, Y = self.integrate(t, t+dt, y0, solver = solver, jac = jac, **solver_keywords)
            # Set new time to last time of integration
            t = ts[-1];
            # sync world variables and y 
            self.update_world_from_y(Y[-1])
            # Run updater functions
            self.Update_World(forceupdate = True)
            self.archive()

        run_time1 = time()
        self.world['history'].update(dict(zip(self.allspcs, Y.T)))
        self.world.update(dict(zip(self.allspcs, Y[-1])))
        self.world['history']['t'] = ts
        self.world['history']['CFACTOR'] = cfactor
        self.world['Y'] = Y
        self.output()
        return run_time1 - run_time0
    
    def add_constraint(self, func):
        """
        Constraints are used to prevent concentration changes. 
        Usually, the intent is to alter reaction rates.
        
        func - a function that takes a mechanism object 
               and a world dictionary. The function should 
               edit the species in world and in the world y
               vector to ensure
               
        def func_example_nox_constraint(mech, world):
            TOTALNOx = world['TOTALNOx']
            y = world['y']
            ind_NO, ind_NO2, ind_HO2NO2, ind_NO3, ind_N2O5 = eval('ind_NO, ind_NO2, ind_HO2NO2, ind_NO3, ind_N2O5', None, world)
            fNOx = TOTALNOx/y[[ind_NO, ind_NO2, ind_HO2NO2, ind_NO3, ind_N2O5]].sum()
            world['NO'] = y[ind_NO] = y[ind_NO] * fNOx
            world['NO2'] = y[ind_NO2] = y[ind_NO2] * fNOx
            world['NO3'] = y[ind_NO3] = y[ind_NO3] * fNOx
            world['N2O5'] = y[ind_N2O5] = y[ind_N2O5] * fNOx
            world['HO2NO2'] = y[ind_HO2NO2] = y[ind_HO2NO2] * fNOx
        
        add_constraint(func_example_nox_constraint)
        """
        self.constraints.append(func)
    
    def get_rates(self, **kwds):
        self.world.update(kwds)
        self.Update_World(self.world)
        self.update_y_from_world()
        self.update_rate_const()
        rate_const = self.world['rate_const'] = self.rate_const
        rates = eval(self.rate_exp, None, self.world)
        
        return zip(self.get_rxn_strs(), rate_const, rates)
    
    def check_rates(self, **kwds):
        for lbl, a, val in self.get_rates(**kwds):
            print('%10.2e,%10.2e,%s' % (a, val, lbl))
        
    def print_world(self, out_keys = None, format = '%.8e', verbose = False):
        t = self.world['t']
        cfactor = self.world.get('CFACTOR', 1)
        if out_keys is None:
            out_keys = ['t'] + [k for i, k in self.monitor_expr]
        if verbose:
            for ti in arange(t.size):
                print('%.1f%%.' % (float(ti) / (t.size - 1) * 100),)
                for k in out_keys:
                    print(('%s=' + format) % (k, self.world[k][ti] / (1 if k == 't' else cfactor)),)
                print()
        else:
            print(','.join(out_keys))
            for ti in arange(t.size):
                print(','.join([format % self.world[k][ti] for k in out_keys]))

    def plot(self, **kwds):
        return _plot(self, self.world, **kwds)