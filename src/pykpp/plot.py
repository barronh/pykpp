#!/usr/bin/env python
from matplotlib import use
from pylab import *
from matplotlib.mlab import csv2rec
from matplotlib.colors import cnames

def plot(mech, world, fig = None, axes = [0.1, 0.3, 0.8, 0.6], ax_props = dict(xlabel = 'hour', ylabel = 'ppb'), **kwds):
    if fig is None:
        fig = figure()
    ax = fig.add_axes(axes)
    setp(ax, **ax_props)
    t = world['t'] / 3600.
    for varkey, varprop in kwds.iteritems():
        CFACTOR = varprop.get('CFACTOR', world['CFACTOR'])
        if 'CFACTOR' in varprop:
            del varprop['CFACTOR']
        varexpr = varprop.get('expr', varkey)
        if 'expr' in varprop:
            del varprop['expr']
        varlabel = varprop.setdefault('label', varkey)
        var = eval(varexpr, None, world) / CFACTOR
        ax.plot(t, var, **varprop)
    
    
    
    handles = [l for l in ax.lines if l.get_label()[:1] != '_']
    labels = [l.get_label() for l in handles]
    fig.legend(handles, labels, loc = 'upper center', bbox_to_anchor = (0.5, -.1), ncol = min(3, np.ceil(len(sys.argv[2:]) / 2.)))
    if 'path' in kwds:
        fig.savefig(kwds['path'])
