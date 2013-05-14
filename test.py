import scipy.integrate as itg
from mech import Mech
out = Mech('small_strato.eqn')
out.run(startt = 12*3600., endt = (12+24*3) * 3600., dt = 3600.)
out.print_world(verbose = True, format = '%.4e')