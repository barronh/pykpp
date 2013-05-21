from numpy import *
import scipy.integrate as itg
from stdfuncs import *
from parse import _parsefile, _reactionstoic, _stoicreaction, _allspc

def Update_RATE(world):
    rate_const_exp = world['rate_const_exp']
    world['rate_const'] = map(lambda x: eval(x, None, world), rate_const_exp)
    #set_printoptions(precision = 4)
    #print world['t'], (array(world['rate_const']) - oldrate) / oldrate * 100
    #set_printoptions(precision = None)

def Update_SUN(world):
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

def Update_World(world):
    Update_SUN(world)
    Update_RATE(world)

def addtojac(rxni, reaction, jac, allspcs):
    rcts = [spc for stc, spc in reaction['reactants']]
    spcorder = [('y[%d]**(%s)' % (allspcs.index(spc_), stc_)).replace('**(1.)', '') for stc_, spc_ in reaction['reactants']]
    for stc, spc in reaction['reactants'] + reaction['products']:
        ispc = allspcs.index(spc)
        for rstc, rspc in reaction['reactants']:
            irspc = allspcs.index(rspc)
            tmpspcorder = [i_ for i_ in spcorder]
            tmpspcorder.pop(rcts.index(rspc))
            expr = ' * '.join(['%s%s*rate_const[%d]' % ('-' if spc in rcts else '+', stc, rxni)] +  tmpspcorder)
            #print ispc, irspc, 0, expr
            #print ispc, irspc, 1, jac[ispc][irspc]
            jac[ispc][irspc] += expr
            #print ispc, irspc, 2, jac[ispc][irspc]
    
class Mech(object):
    def __init__(self, path, verbose = False, keywords = ['hv']):
        self._parsed = _parsefile(path)
        for rxn in self._parsed['EQUATIONS']:
            remove = []
            for spci, (stc, spc) in enumerate(rxn['reactants']):
                if spc.strip() in keywords:
                    remove.append(spci)
            [rxn['reactants'].pop(spci) for spci in remove[::-1]]
        self.allspcs = _allspc(self._parsed)


        nspcs = len(self.allspcs)
        self.summary()
        self.world = {}
        if 'worldupdater' in self._parsed:
            tmp = {}
            exec(self._parsed['WORLDUPDATER'][0], None, tmp)
            self.Update_World = tmp[tmp.keys()[0]]
        else:
            self.Update_World = Update_World
        for k, v in self._parsed['INITVALUES'].asList():
            self.world[k] = eval(v, None, self.world)
                    
        self.dy_stoic = {}

        for spc in self.allspcs:
            if spc not in self.world:
                self.world[spc] = self.world.get('DEFAULTCONC', 0.)
            self.dy_stoic[spc] = _reactionstoic(spc, self._parsed['EQUATIONS'])
        self.dy_exp = ['0' for i_ in range(nspcs)]
        for si, spc in enumerate(self.allspcs):
            for ri, stoic in self.dy_stoic[spc].iteritems():
                self.dy_exp[si] += (' + %f * rates[%d]' % (stoic, ri)).replace('+ -1.000000 * ', '-').replace('+ 1.000000 * ', '+')

        rate_const_exp = self.world['rate_const_exp'] = []
        self.rate_exp = []
        self.drate_exp = [['' for i_ in range(nspcs)] for j_ in range(nspcs)]
        for rxni, reaction in enumerate(self._parsed['EQUATIONS']):
            rate_const_exp.append(reaction['rate'])
            spcorder = [('y[%d]**(%s)' % (self.allspcs.index(spc_), stc_)).replace('**(1.)', '') for stc_, spc_ in reaction['reactants']]
            self.rate_exp.append(' * '.join(['rate_const[%d]' % rxni] + spcorder))
            addtojac(rxni, reaction, self.drate_exp, self.allspcs)
        #print 'Jac'
        #for ri, row in enumerate(self.drate_exp):
        #    for ci, colexpr in enumerate(row):
        #        print ri, ci, colexpr
    def summary(self, verbose = False):
        if verbose:
            spcstr = ', '.join(self.allspcs)
            rxnstr = '\n' + '\n'.join(self.get_rxn_strs)
        else:
            spcstr = ''
            rxnstr = ''
        print 'Species: %d%s\nReactions: %d%s' % (len(self.allspcs), spcstr, len(self.get_rxn_strs()), rxnstr)
    def get_rxn_strs(self):
        return [' + '.join(['*'.join(stcspc) for stcspc in rxn['reactants']]) + '->' + ' + '.join(['*'.join(stcspc) for stcspc in rxn['products']]) + ': ' + rxn['rate'] + ';' for rxn in self._parsed['EQUATIONS']]
    
    def print_rxns(self):
        print '\n'.join(self.get_rxn_strs())
        
    def dy(self, y, t):
        world = self.world
        allspcs = self.allspcs
        nspcs = len(allspcs)
        world['t'] = t
        out = zeros(nspcs, dtype = 'd')
        self.Update_World(world)
        rate_const = world['rate_const']
        tmp = {}
        tmp.update(locals())
        tmp.update(world)
        rates = world['rates'] = map(lambda expr: eval(expr, None, tmp), self.rate_exp)
        for si, dy_exp in enumerate(self.dy_exp):
            out[si] = eval(dy_exp)
                    
        return out[:]

    def ddy(self, y, t):
        world = self.world
        world['t'] = t
        rate_const = world['rate_const']
        tmp = {}
        tmp.update(locals())
        nspcs = len(self.allspcs)
        out = zeros((nspcs, nspcs))
        for ri, row in enumerate(self.drate_exp):
            for ci, colexpr in enumerate(row):
                if colexpr != '':
                    out[ri, ci] = eval(colexpr, None, tmp)
                
        return out[:]

    def run(self, solver = 'lsoda', startt = None, endt = None, dt = None, jac = True, **solver_keywords):
        from time import time
        run_time0 = time()
        maxstepname = dict(lsoda = 'mxstep')
        default_solver_params = dict(atol = 1e-3, rtol = 1e-4, maxstep = 1000)
        for k, v in default_solver_params.iteritems():
            if k == 'maxstep':
                k = maxstepname.get(solver, 'max_step')
            solver_keywords.setdefault(k, v)
        
        if startt is None: startt = self.world.get('TSTART', 12*3600)
        if endt is None: endt = self.world.get('TEND', (12 + 24)*3600)
        if dt is None: dt = self.world.get('DT', 3600)
        
        y0 = array([eval(spc, None, self.world) for spc in self.allspcs])
        self.world['t'] = startt
        self.Update_World(self.world)
        
        if solver == 'lsoda':
            # Old Method
            ts = arange(startt, endt + dt, dt)
            Y, infodict = itg.odeint(self.dy, y0, ts, Dfun = self.ddy if jac else None, mxords = 2, full_output = True, **solver_keywords)
            self.infodict = infodict
        else:
            # New method
            Y = y0
            ts = array([startt])
            r = itg.ode(lambda t_, y_: self.dy(y_, t_), jac = lambda t_, y_: self.ddy(y_, t_) if jac else None)
            r.set_integrator(solver, method = 'bdf', **solver_keywords)
            r.set_initial_value(y = y0, t = startt)
            while r.t < endt:
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
        if out_keys is None:
            out_keys = ['t'] + [k for k, v in self.world.iteritems() if hasattr(v, 'shape') and v.shape == t.shape and k != 't']
        if verbose:
            for ti in arange(t.size):
                print '%.1f%%.' % (float(ti) / (t.size - 1) * 100),
                for k in out_keys:
                    print ('%s=' + format) % (k, self.world[k][ti] / self.world.get('CFACTOR', 1)),
                print
        else:
            print ','.join(out_keys)
            for ti in arange(t.size):
                print ','.join([format % self.world[k][ti] for k in out_keys])
        