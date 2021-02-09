from __future__ import print_function
import os
from glob import glob
from .mech import Mech

from argparse import ArgumentParser, Action, RawDescriptionHelpFormatter

modelpaths = [os.path.basename(p) for p in glob(os.path.join(os.path.dirname(__file__), 'models', '*.kpp'))]
parser = ArgumentParser(description = 'pykpp Argument Parsing', usage = """Usage: python -m pykpp mech.def

Where mech.pykpp is a single file that contains initial
values and reaction definitions according to the KPP
format.
""", formatter_class = RawDescriptionHelpFormatter)

parser.add_argument('paths', nargs='+', help='path to a pykpp definition file')
    
parser.add_argument("-v", "--verbose", dest = "verbose", action = "store_true", default = False, help = "Show extended output")

parser.add_argument("-a", "--atol", dest = "atol", type = float, default = 1e-3, help = "absolute tolerance")

parser.add_argument("-r", "--rtol", dest = "rtol", type = float, default = 1e-4, help = "relative tolerance")

parser.add_argument("--tstart", dest = "tstart", type = float, default = None, help = "Start time")

parser.add_argument("--tend", dest = "tend", type = float, default = None, help = "End time")

parser.add_argument("--norun", dest = "norun", action = 'store_true', default = False, help = "Don't run")

parser.add_argument("--dt", dest = "dt", type = float, default = None, help = "Time step")

parser.add_argument("--no-jacobian", dest = "jacobian", action="store_false", default = True, help = "disable use of jacobian")

parser.add_argument("--no-autofuncs", dest = "add_default_funcs", action="store_false", default = True, help = "disable use of automatically identified functions.")

parser.add_argument("-k", "--keywords", dest = "keywords", type=str, default = 'hv,PROD,EMISSION', help = "List of keywords to be ignored in reactants or products (comma delimited; default='hv')")

parser.add_argument("-o", "--outpath", dest = "outpath", type=str, default = None, help = "Output path.")

parser.add_argument("-m", "--monitor", dest = "monitor", type=str, default = None, help = "Extra monitor values (comma separated string).")

parser.add_argument("-s", "--solver", dest = "solver", default = None, help = "solver (default: lsoda; vode; zvode; dopri5; dop853) with optional keywords comma delimited (e.g., lsoda,atol=1e-2,rtol=1e-5,max_order_s=2,max_order_ns=2,nsteps=1000)")

parser.add_argument("-c", "--code", dest = "code", default = "", help = "code to solve (exec) after model is compiled and run (unless --norun); out is a keyword for the mechanism that has been generated")

parser.epilog = """

functions for world updating can be explored using pydoc pykpp.updaters

functions for reaction rates and world updating can be explored using pydoc pykpp.stdfuncs

run existing models (pykpp modelfile.kpp):\n\t""" + '\n\t'.join(modelpaths)
def main(globals={}, options=None):
    if options is None:
        options = parser.parse_args()

    if len(options.paths) == 0:
        parser.print_help()
        exit()
    
    outs = []
    for arg in options.paths:
        out = Mech(arg, verbose = options.verbose, keywords = [k_.strip() for k_ in options.keywords.split(',')], add_default_funcs = options.add_default_funcs)
        if options.monitor is not None:
            out.monitor = tuple([(None if k not in out.allspcs else out.allspcs.index(k), k) for k in options.monitor.split(',')]) + out.monitor
        if not options.norun:
            if not options.solver is None:
                solver = options.solver.split(',')[0]
                solver_keywords = eval('dict(' + ','.join(options.solver.split(',')[1:])+ ')')
            else:
                solver = options.solver
                solver_keywords = {}
            runtime = out.run(tstart = options.tstart, tend = options.tend, dt = options.dt, solver = solver, jac = options.jacobian, **solver_keywords)
            print('Solved in %f seconds' % runtime)
            
        outs.append(out)

    if os.path.exists(options.code):
        options.code = file(options.code, 'r').read()
    exec(options.code)
    globals.update(locals())

if __name__ == '__main__':
    main(globals=globals())
