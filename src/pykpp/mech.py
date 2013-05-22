from numpy import *
from scipy.constants import *
import scipy.integrate as itg
from stdfuncs import *
from parse import _parsefile, _reactionstoic, _stoicreaction, _allspc

def Update_RATE(mech, world):
    """
    Update_RATE is the default rate updater that is called when
    rates are updated for the integrator solution
    
    1) Call update_func_world(world)
    2) Set world['rate_const'] equal to evaluated world['rate_const_exp']
    """
    update_func_world(world)
    mech.rate_const = map(lambda x: eval(x, None, world), mech.rate_const_exp)

def Update_SUN(world):
    """
    Update_SUN is the default updater for the solar zenith
    angle based on 'SunRise' (default: 4.5) and 'SunSet' (default: 19.5)
    """
    t = world['t']
    SunRise = world.get('SunRise', 4.5)
    SunSet  = world.get('SunSet', 19.5)
    Thour = t/3600.0
    Tlocal = Thour - (int(Thour)/24)*24.

    if (Tlocal >= SunRise) and (Tlocal <= SunSet):
       Ttmp = (2.0 * Tlocal - SunRise - SunSet) / (SunSet - SunRise)
       if (Ttmp > 0):
          Ttmp =  Ttmp * Ttmp
       else:
          Ttmp = -Ttmp * Ttmp

       SUN = ( 1.0 + cos(pi * Ttmp) )/2.0
    else:
       SUN = 0.0
    world['SUN'] = SUN

def Update_World(mech, world):
    """
    Calls Update_SUN(world) and Update_RATE(mech, world)
    
    see these functions for more details
    """
    Update_SUN(world)
    Update_RATE(mech, world)

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
    spcorder = [('y[%d]**(%s)' % (allspcs.index(spc_), stc_)).replace('**(1.)', '') for stc_, spc_ in reaction['reactants']]
    
    # Iterate over all stoichiometries and spcs names for
    # reactions and products in reaction
    for stc, spc in reaction['reactants'] + reaction['products']:
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
            expr = ' * '.join(['%s%s*rate_const[%d]' % ('-' if spc in rcts else '+', stc, rxni)] +  tmpspcorder)
            
            # Add expression to the jacobian
            jac[ispc][irspc] += expr


class Mech(object):
    """
    Mech object must be able to interpret kpp inputs
    and translate the results to create a reaction
    function (Mech.dy), a jacobian (Mech.ddy), and 
    run the model (Mech.run)
    """
    def __init__(self, path, verbose = False, keywords = ['hv']):
        """
        path     - path to kpp inputs
        verbose  - add printing
        keywords - ignore certain keywords from reactants in
                   calculate reaction rates
        """
        
        # Parse kpp inputs and store the results
        self._parsed = _parsefile(path)
        
        # For rxn dict-like objects in EQUATIONS
        # remove reactants that are keywords
        for rxn in self._parsed['EQUATIONS']:
            remove = []
            for spci, (stc, spc) in enumerate(rxn['reactants']):
                if spc.strip() in keywords:
                    remove.append(spci)
            [rxn['reactants'].pop(spci) for spci in remove[::-1]]
        
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
        self.world = {}
        
        # If no WORLDUPDATER is provided,
        # use the default
        if 'WORLDUPDATER' in self._parsed:
            tmp = {}
            exec(self._parsed['WORLDUPDATER'][0], None, tmp)
            self.Update_World = tmp[tmp.keys()[0]]
        else:
            self.Update_World = Update_World

        # Execute the INITVALUES code in the context
        # of the world
        exec(self._parsed['INITVALUES'][0], None, self.world)
        
        # Multiply all input values by CFACTOR
        cfactor = self.world.get('CFACTOR', 1.)
        
        for k, v in self.world.iteritems():
            if k not in ('CFACTOR', 'TEMP', 'P'): self.world[k] = v * cfactor
        
        # Execute INIT code in the context of world
        exec(self._parsed['INIT'][0], None, self.world)
        
        
        # Store a temporary copy of the stoichiometry in each reaction
        dy_stoic = {}
        for spc in self.allspcs:
            if spc not in self.world:
                self.world[spc] = self.world.get('ALL_SPEC', 0.) * cfactor
            dy_stoic[spc] = _reactionstoic(spc, self._parsed['EQUATIONS'])
        
        # Add each reaction rate expression to the dy
        # for the species in the reacton
        self.dy_exp = ['0' for i_ in range(nspcs)]
        for si, spc in enumerate(self.allspcs):
            for ri, stoic in dy_stoic[spc].iteritems():
                self.dy_exp[si] += (' + %f * rates[%d]' % (stoic, ri)).replace('+ -1.000000 * ', '-').replace('+ 1.000000 * ', '+')
        
        # Create an empty copy of rate constant expressions (rate_const_exp)
        # and rate expressions (rate_exp; e.g., rate_const_exp * y[0] * y[1])
        # and rate derivatives with respect to a species (self.drate_exp)
        rate_const_exp = self.rate_const_exp = []
        self.rate_exp = []
        self.drate_exp = [['' for i_ in range(nspcs)] for j_ in range(nspcs)]
        for rxni, reaction in enumerate(self._parsed['EQUATIONS']):
            rate_const_exp.append(reaction['rate'])
            spcorder = [('y[%d]**(%s)' % (self.allspcs.index(spc_), stc_)).replace('**(1.)', '') for stc_, spc_ in reaction['reactants']]
            self.rate_exp.append(' * '.join(['rate_const[%d]' % rxni] + spcorder))
            addtojac(rxni, reaction, self.drate_exp, self.allspcs)
    
    def summary(self, verbose = False):
        """
        return a string with species number, reaction number, and 
        string representations of each
        """
        if verbose:
            spcstr = ', '.join(self.allspcs)
            rxnstr = '\n' + '\n'.join(self.get_rxn_strs())
        else:
            spcstr = ''
            rxnstr = ''
        return 'Species: %d%s\nReactions: %d%s' % (len(self.allspcs), spcstr, len(self.get_rxn_strs()), rxnstr)
    
    def get_rxn_strs(self):
        """
        Generate reaction string representations from parsed
        reaction  objects
        """
        return [' + '.join(['*'.join(stcspc) for stcspc in rxn['reactants']]) + '->' + ' + '.join(['*'.join(stcspc) for stcspc in rxn['products']]) + ': ' + rxn['rate'] + ';' for rxn in self._parsed['EQUATIONS']]
    
    def print_rxns(self):
        """
        Print reaction stings joined by line returns
        """
        print '\n'.join(self.get_rxn_strs())
        
    def dy(self, y, t):
        """
        Function that returns the change in species
        with respect to time
        """
        world = self.world
        allspcs = self.allspcs
        nspcs = len(allspcs)
        world['t'] = t
        out = zeros(nspcs, dtype = 'd')
        self.Update_World(self, world)
        rate_const = self.rate_const
        tmp = {}
        tmp.update(locals())
        tmp.update(world)
        rates = world['rates'] = map(lambda expr: eval(expr, None, tmp), self.rate_exp)
        for si, dy_exp in enumerate(self.dy_exp):
            out[si] = eval(dy_exp)
                    
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
        out = zeros((nspcs, nspcs))
        for ri, row in enumerate(self.drate_exp):
            for ci, colexpr in enumerate(row):
                if colexpr != '':
                    out[ri, ci] = eval(colexpr, None, tmp)
                
        return out[:]

    def run(self, solver = None, tstart = None, tend = None, dt = None, jac = True, **solver_keywords):
        """
        Load solvers with Mech object function (Mech.dy), jacobian (Mech.ddy),
        and mechanism specific options
        """
        
        if solver is None:
            solver = self._parsed['INTEGRATOR']
        from time import time
        run_time0 = time()
        maxstepname = dict(lsoda = 'mxstep')
        default_solver_params = dict(atol = 1e-3, rtol = 1e-4, maxstep = 1000)
        for k, v in default_solver_params.iteritems():
            if k == 'maxstep':
                k = maxstepname.get(solver, 'max_step')
            solver_keywords.setdefault(k, v)
        
        if tstart is None: tstart = self.world.get('TSTART', 12*3600)
        if tend is None: tend = self.world.get('TEND', (12 + 24)*3600)
        if dt is None: dt = self.world.get('DT', 3600)
        
        y0 = array([eval(spc, None, self.world) for spc in self.allspcs])
        self.world['t'] = tstart
        self.Update_World(self, self.world)
        
        if solver == 'lsoda':
            # Old Method
            ts = arange(tstart, tend + dt, dt)
            Y, infodict = itg.odeint(self.dy, y0, ts, Dfun = self.ddy if jac else None, mxords = 2, full_output = True, **solver_keywords)
            self.infodict = infodict
        else:
            # New method
            Y = y0
            ts = array([tstart])
            r = itg.ode(lambda t_, y_: self.dy(y_, t_), jac = lambda t_, y_: self.ddy(y_, t_) if jac else None)
            r.set_integrator(solver, method = 'bdf', **solver_keywords)
            r.set_initial_value(y = y0, t = tstart)
            while r.t < tend:
                r.integrate(r.t+dt)
                Y = vstack([Y, r.y])
                ts = append(ts, r.t)
            #    print r.t, r.y
        run_time1 = time()
        self.world.update(dict(zip(self.allspcs, Y.T)))
        self.world['t'] = ts
        self.world['Y'] = Y
        return run_time1 - run_time0
        
    def print_world(self, out_keys = None, format = '%.8e', verbose = False):
        t = self.world['t']
        if out_keys is None and 'MONITOR' in self._parsed.keys():
            out_keys = self._parsed['MONITOR'][0].replace(' ', '').split(';')
            out_keys = [k_ for k_ in out_keys if k_ != '']
        cfactor = self.world.get('CFACTOR', 1)
        if out_keys is None:
            out_keys = ['t'] + [k for k, v in self.world.iteritems() if hasattr(v, 'shape') and v.shape == t.shape and k != 't']
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
        