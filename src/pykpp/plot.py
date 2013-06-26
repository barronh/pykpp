#!/usr/bin/env python
from matplotlib import use
from pylab import *
from matplotlib.mlab import csv2rec
from matplotlib.colors import cnames

def plot_from_file(path, **kwds):
    """
    See plot for more details.
    
    Note that kwds that are explicit in plot
    will be removed from kwds before plot
    evaluates these keywords
    """
    data = csv2rec(path)
    data = dict([(k, data[k]) for k in data.dtype.names])
    return plot(None, world = data, **kwds)

def plot(mech, world, fig = None, axes = [0.1, 0.3, 0.8, 0.6], ax_props = dict(xlabel = 'hour', ylabel = 'ppb'), path = None, **kwds):
    """
    
    mech - Not used unless now kwds are provided
    world - dictionary to find data in
    fig  - optional figure to append
    axes - axes bounding box or axes object to be added to fig
    ax_props - properties of the axes
    kwds - each named keyword is a set of options to plot;
           if kwds == {}; then kwds = dict([(k, {}) for i, k in mech.monitor])
    """
    if fig is None:
        fig = figure()
    ax = fig.add_axes(axes)
    setp(ax, **ax_props)
    t = world['t'] / 3600.
    if kwds == {}:
        kwds = dict([(k, {}) for i, k in mech.monitor])
    
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
    fig.legend(handles, labels, loc = 'lower center', bbox_to_anchor = (0.5, 0), ncol = min(3, len(kwds)))
    if not path is None:
        fig.savefig(path)
