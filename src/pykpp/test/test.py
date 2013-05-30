from matplotlib.mlab import csv2rec, rec2csv
import numpy as np

for model in ['cbm4']:
    exit_code = os.system('kpp %(model)s && make -f Makefile_%(model)s && ./%(model)s' % globals())
    if exit_code != 0:
        raise Exception('KPP failed; see above')
    os.system('python -m pykpp %(model)s')
    
    kpp = csv2rec('%(model)s.dat' % globals(), delimiter = ' ')
    pykpp = csv2rec('$(model)s.pykpp.dat' % globals(), delimiter = ',')
    diff = pykpp.copy()
    pct = pykpp.copy()
    keys = set(kpp.dtype.names).intersection(pykpp.dtype.names)
    notkeys = set(pykpp.dtype.names).difference(kpp.dtype.names)
    notkeys.remove('t')
    print notkeys
    for k in notkeys:
        diff[k] = np.nan
        pct[k] = np.nan
    
    for k in keys:
        diff[k] = pykpp[k] - kpp[k][:pykpp[k].shape[0]]
        pct[k] = diff[k] / kpp[k][:pykpp[k].shape[0]] * 100
    
    rec2csv(diff, '%(model)s.diff.csv' % globals(), delimiter = ',')
    rec2csv(pct, '%(model)s.pct.csv' % globals(), delimiter = ',')
