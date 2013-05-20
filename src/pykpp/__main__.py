from optparse import OptionParser
from mech import Mech

parser = OptionParser()
parser.set_usage("""Usage: python -m pykpp mech.def

Where mech.def is a single file that contains initial
values and reaction definitions according to the KPP
format.
""")

(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_usage()
    exit()
    
for arg in args:
    out = Mech(arg, verbose = True)
    runtime = out.run(solver = 'lsoda', jac = False)
    print 'Solved in %f seconds' % runtime
    out.print_world(format = '%.4e', verbose = True)
    #for rxni, (rxn, rate) in enumerate(zip(out.get_rxn_strs(), out.world['rate_const'])):
    #    print '%d: %.2e,' % (rxni + 1, rate)
