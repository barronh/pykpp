from matplotlib.mlab import csv2rec, rec2csv
import numpy as np
import os
import shutil
from StringIO import StringIO
kpphome = os.environ.get('KPP_HOME', '.')
_allmodels = ['small', 'cbm4']

def runmodels(models = _allmodels, pykpp = True, kpp = True, verbose = False):
    for model in models:
        os.system('rm -rf %(model)s' % locals())
        os.mkdir(model)
        path = os.path.join(kpphome, 'examples', model + '_f90.kpp')
        if not os.path.exists(path):
            path = os.path.join(kpphome, 'models', model + '.def')
        if not os.path.exists(path):
            raise IOError('No kpp or def file could be found for %(model)s' % locals())
        modeldef = model + '.kpp'
        shutil.copy(path, os.path.join(model, modeldef))
        if kpp:
            exit_code = os.system('cd %(model)s && kpp %(modeldef)s && make -f Makefile_%(model)s && ./%(model)s.exe' % locals())
            if exit_code != 0:
                raise Exception('KPP failed; see above')
        if pykpp: os.system('cd %(model)s && python -m pykpp %(modeldef)s' % locals())


def makediffs(models = _allmodels, verbose = False):
    for model in models:
        kpp = csv2rec(os.path.join(model, model + '.dat'), delimiter = ' ')
        pykpp = csv2rec(os.path.join(model, model + '.pykpp.dat'), delimiter = ',')
        diff = pykpp.copy()
        pct = pykpp.copy()
        keys = set(kpp.dtype.names).intersection(pykpp.dtype.names)
        notkeys = set(pykpp.dtype.names).difference(kpp.dtype.names)
        notkeys.remove('t')
        for k in notkeys:
            diff[k] = np.nan
            pct[k] = np.nan
    
        for k in keys:
            diff[k] = pykpp[k] - kpp[k][:]
            pct[k] = diff[k] / kpp[k][:] * 100
        diff['t'] = pykpp['t'] - (kpp['time'] * 3600. + pykpp['t'][0])
        pct['t'] = diff['t'] / (kpp['time'] * 3600. + pykpp['t'][0]) * 100
        
        rec2csv(diff, os.path.join(model, model + '.diff.csv'), delimiter = ',')
        rec2csv(pct, os.path.join(model, model + '.pct.csv'), delimiter = ',')

def checkmodels(models = _allmodels, verbose = False):
    for model in models:
        kpp = csv2rec(os.path.join(model, model + '.dat'), delimiter = ' ')
        pykpp = csv2rec(os.path.join(model, model + '.pykpp.dat'), delimiter = ',')
        diff = csv2rec(os.path.join(model, model + '.diff.csv'), delimiter = ',')
        pct = csv2rec(os.path.join(model, model + '.pct.csv'), delimiter = ',')
        keys = set(kpp.dtype.names).intersection(pykpp.dtype.names)
        keys.add('time')
        output = StringIO('')
        print  >> output, model, 'Check' 
        print  >> output, '%20s,%s,%s,%s,%s,%10s,%10s,%10s,%10s' % ('Key', 'Type  ', 'Time pct', 'Time    ', 'Idx', 'KPP', 'PYKPP', 'Absolute', 'Percent')
        trans = dict(time = 't')
        check = True
        for k in keys:
            thisdiff = diff[trans.get(k, k)]
            thispct = pct[trans.get(k, k)]
            thispykpp = pykpp[trans.get(k, k)]
            thiskpp = kpp[k]
            if k == 'time':
                thiskpp *= 3600.
                thiskpp += thispykpp[0]
            thismaxdiff = thispykpp.max() - thiskpp.max()
            thismaxpct = thismaxdiff / thiskpp.max() * 100
            absidx = np.abs(thisdiff).argmax()
            pctidx = np.abs(thispct).argmax()
            abspass = ('PASS' if np.abs(thispct[absidx]) < 2 else 'FAIL') if np.abs(thisdiff[absidx]) > .001 else 'PASS'
            pctpass = ('PASS' if np.abs(thispct[pctidx]) < 2 else 'FAIL') if np.abs(thisdiff[pctidx]) > .001 else 'PASS'
            endpass = ('PASS' if np.abs(thispct[-1]) < 2 else 'FAIL') if np.abs(thisdiff[-1]) > .001 else 'PASS'
            maxpass = ('PASS' if np.abs(thismaxpct) < 2 else 'FAIL') if np.abs(thismaxdiff) > .001 else 'PASS'
            check = check and all([abspass == 'PASS', pctpass == 'PASS', endpass == 'PASS', maxpass == 'PASS'])
            print >> output, '%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%' % (k.upper(), 'MaxAbs', 100 * absidx / len(pykpp['t']), pykpp['t'][absidx], absidx, thiskpp[absidx], thispykpp[absidx], thisdiff[absidx], thispct[absidx]), abspass
            print >> output, '%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%' % (k.upper(), 'MaxPct', 100 * pctidx / len(pykpp['t']), pykpp['t'][pctidx], absidx, thiskpp[pctidx], thispykpp[pctidx], thisdiff[pctidx], thispct[pctidx]), pctpass
            print >> output,  '%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%' % (k.upper(), 'Max   ', -1, -1, -1, thiskpp.max(), thispykpp.max(), thismaxdiff, thismaxpct), maxpass
            print >> output,  '%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%' % (k.upper(), 'Ending', 100, pykpp['t'][-1], -1, thiskpp[-1], thispykpp[-1], thisdiff[-1], thispct[-1]), endpass
        if not check:
            output.seek(0, 0)
            print output.read()
        print model, 'PASS' if check else 'FAIL'

#for model in _allmodels:
#    runmodels([model])
#    makediffs([model])
#    checkmodels([model])
if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage("""Usage: python -m pykpp.test [-v]""")

    parser.add_option("-v", "--verbose", dest = "verbose", action = "store_true", default = False, help = "Show extended output")
    options, args = parser.parse_args()
    runmodels(verbose = options.verbose)
    makediffs(verbose = options.verbose)
    checkmodels(verbose = options.verbose)