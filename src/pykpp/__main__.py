from optparse import OptionParser
from mech import Mech

parser = OptionParser()
parser.set_usage("""Usage: python -m pykpp mech.def

Where mech.pykpp is a single file that contains initial
values and reaction definitions according to the KPP
format.
""")

parser.add_option("-a", "--atol", dest = "atol", type = "float", default = 1e-3, help = "absolute tolerance")

parser.add_option("-r", "--rtol", dest = "rtol", type = "float", default = 1e-4, help = "relative tolerance")

parser.add_option("", "--tstart", dest = "tstart", type = "string", default = '0.', help = "Start time")

parser.add_option("", "--tend", dest = "tend", type = "string", default = '10.', help = "End time")

parser.add_option("", "--dt", dest = "dt", type = "string", default = '1', help = "Time step")

parser.add_option("-j", "--jacobian", dest = "jacobian", action="store_true", default = False, help = "Enable use of jacobian")

parser.add_option("-k", "--keywords", dest = "keywords", type="string", default = 'hv', help = "List of keywords to be ignored in reactants or products (comma delimited; default='hv')")

parser.add_option("-s", "--solver", dest = "solver", default = 'lsoda', help = "solver (default: lsoda; vode; zvode; dopri5; dop853)")

(options, args) = parser.parse_args()

if len(args) == 0:
    parser.print_help()
    exit()
    
for arg in args:
    out = Mech(arg, verbose = True, keywords = [k_.strip() for k_ in options.keywords.split(',')])
    runtime = out.run(solver = options.solver, jac = options.jacobian, atol = options.atol, rtol = options.rtol)
    print 'Solved in %f seconds' % runtime
    out.print_world(format = '%.4e', verbose = True)
    #for rxni, (rxn, rate) in enumerate(zip(out.get_rxn_strs(), out.world['rate_const'])):
    #    print '%d: %.2e,' % (rxni + 1, rate)
