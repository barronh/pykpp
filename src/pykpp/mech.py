from numpy import *
import scipy.integrate as itg
from parse import _parsefile, _reactionstoic, _stoicreaction, _allspc

class Mech(object):
    def __init__(self, path):
        self._parsed = _parsefile(path)
        self.allspcs = _allspc(self._parsed)
        self.world = dict(self._parsed['initvalues'].asList())
        self.dy_stoic = {}
        self._solver = itg.odeint
        self._rate_exp_updated = -inf
        self._rate_exp_dt = 3600. * .25
        for spc in self.allspcs:
            if spc not in self.world:
                self.world[spc] = 0.
            else:
                self.world[spc] = eval(self.world[spc])
            self.dy_stoic[spc] = _reactionstoic(spc, self._parsed['reactions'])

        self.rate_exp = []
        for reaction in self._parsed['reactions']:
            self.rate_exp.append(' * '.join([reaction['rate']] + ['y[%d]**(%s)' % (self.allspcs.index(spc), stc) for stc, spc in reaction['reactants']]).replace('**(1.)', ''))

    def Update_SUN(self, t):
    
        SunRise = self.world.get('SunRise', 4.5)
        SunSet  = self.world.get('SunSet', 19.5)
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
        self.world['SUN'] = SUN

    def dy(self, y, t):
        world = self.world
        allspcs = self.allspcs
        world['y'] = y
        world['t'] = t
        out = zeros(len(allspcs))
        self.Update_SUN(t)
        rates = map(lambda expr: eval(expr, None, world), self.rate_exp)
        for si, spc in enumerate(allspcs):
            for ri, stoic in self.dy_stoic[spc].iteritems():
                out[si] += rates[ri] * stoic
                    
        return out[:]

    def run(self, startt = 0, endt = 3600 * 24, dt = 3600):
        y0 = array([eval(spc, None, self.world) for spc in self.allspcs])
        t = arange(startt, endt + dt, dt)
        Y = self._solver(self.dy, y0, t)
        self.world.update(dict(zip(self.allspcs, Y.T)))
        self.world['t'] = t

    def print_world(self, format = '%.8e', verbose = False):
        t = self.world['t']
        tout_keys = ['t'] + [k for k, v in self.world.iteritems() if hasattr(v, 'shape') and v.shape == t.shape and k != 't']
        if verbose:
            for ti in arange(t.size):
                print '%.1f%%.' % (float(ti) / t.size * 100),
                for k in tout_keys:
                    print ('%s=' + format) % (k, self.world[k][ti]),
                print
        else:
            print ','.join(tout_keys)
            for ti in arange(t.size):
                print ','.join([format % self.world[k][ti] for k in tout_keys])
        