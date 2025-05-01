#!/usr/bin/env python
import matplotlib.pyplot as plt
import pandas as pd
from warnings import warn


def plot_from_file(path, **kwds):
    """
    See plot for more details.

    Note that kwds that are explicit in plot
    will be removed from kwds before plot
    evaluates these keywords
    """
    data = pd.read_csv(path, sep=',')
    data = dict([(k, data[k]) for k in data.dtype.names if k != 'CFACTOR'])
    return plot(None, world=dict(history=data), **kwds)


def plot(mech, world, fig=None, ax=None, ax_props=None, path=None, **kwds):
    """
    mech - Not used unless now kwds are provided
    world - dictionary to find data in
    fig  - optional figure to append
    ax - axes bounding box or axes object to be added to fig
    ax_props - properties of the axes
    kwds - each named keyword is a set of options to plot;
           if kwds == {}; then kwds = dict([(k, {}) for i, k in mech.monitor])
    """
    if ax is None:
        ax = (0.1, 0.3, 0.8, 0.6)
    if ax_props is None:
        ax_props = dict(xlabel='hour', ylabel='unknown')
    if isinstance(ax, plt.matplotlib.axes.Axes):
        fig = ax.figure
    elif fig is None:
        fig = plt.figure()

    ax = fig.add_axes(ax)
    plt.setp(ax, **ax_props)
    t = world['history']['t'] / 3600.
    if kwds == {}:
        kwds = dict([(k, {}) for i, k in mech.monitor])

    for varkey, varprop in kwds.items():
        CFACTOR = varprop.get(
            'CFACTOR',
            world['history'].get('CFACTOR', world.get('CFACTOR', 1))
        )
        if 'CFACTOR' in varprop:
            del varprop['CFACTOR']
        varexpr = varprop.get('expr', varkey)
        if 'expr' in varprop:
            del varprop['expr']
        varprop.setdefault('label', varkey)
        try:
            var = eval(varexpr, None, world['history']) / CFACTOR
            ax.plot(t, var, **varprop)
        except Exception:
            warn('Skipping %s' % varexpr)

    handles = [_l for _l in ax.lines if _l.get_label()[:1] != '_']
    labels = [_l.get_label() for _l in handles]
    fig.legend(
        handles, labels, loc='lower center', bbox_to_anchor=(0.5, 0),
        ncol=min(3, len(kwds))
    )
    if path is not None:
        fig.savefig(path)
    return fig
