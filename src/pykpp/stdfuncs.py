#__all__ = ['solar_declination', 'solar_noon_local', 'solar_noon_utc', \
#           'TUV_J', 'h2o_from_rh_and_temp']

__doc__ = """
stdfuncs 

stdfuncs provides an natural environment and functions for evaluating rate 
constants from rate expressions in the context of that environment.

Environment:
    The uninitialized environment includes a default temperature (298. K), pressure (101325 Pa), M, O2, N2, and H2.

Rate Evaluation:
    Because many models use different forms (e.g., signs), several models native
    forms are provided:

    GEOS-Chem associates different forms with a letter.  The basic form is provided 
    by GEOS_STD, and then other forms are provided by GEOS_\S+ where \S+ is one or
    more non-white space (e.g., GEOS_A or GEOS_B).

    MOZART4 uses normal mathematical expressions and several special functions.
    Special functions include a troe form (MZ4_TROE) and numbered user expressions
    (MZ4_USR\d+) where \d+ is one or more digits.

    CMAQ uses special characters to key to numbered forms documented in the Science
    documentation appendix. CMAQ_1to4 is the most basic where all 4 forms can be
    represented by supplying a 0 for one or more parameters. CMAQ_\d+ where \d+ is 5
    to 10 are special functions (10 is the troe fall off equation).

    CHIMERE uses basic expressions and two special functions. Special functions are
    two forms of the TROE equation (CHIMERE_TROE and CHIMERE_MTROE).


"""

from datetime import datetime
import numpy as np
from warnings import warn
from numpy import *
from scipy.constants import *
from matplotlib.mlab import csv2rec
from tuv.tuv5pt0 import TUV_J
import funcs

from funcs.geoschem import *
from funcs.mcm import *
from funcs.am3 import *
from funcs.mozart4 import *
from funcs.chimere import *
from funcs.cmaq import *
from funcs.camx import *
from funcs.kpp import *
from updaters import *
del __version__
__all__ = ['TUV_J', 'update_func_world', 'solar_noon_local', 'solar_noon_utc', 'h2o_from_rh_and_temp', 'initstdenv', 'solar_noon']
import scipy.constants
__all__ += scipy.constants.__all__
import funcs.geoschem
__all__ += funcs.geoschem.__all__
import funcs.mcm
__all__ += funcs.mcm.__all__
import funcs.am3
__all__ += funcs.am3.__all__
import funcs.mozart4
__all__ += funcs.mozart4.__all__
import funcs.chimere
__all__ += funcs.chimere.__all__
import funcs.cmaq
__all__ += funcs.cmaq.__all__
import funcs.camx
__all__ += funcs.camx.__all__
import funcs.kpp
__all__ += funcs.kpp.__all__
import updaters
__all__ += updaters.__all__
__all__ += ['datetime']
try:
    boltz  = Boltzmann / centi**2 * kilo # in erg/K
except:
    boltz  = constants.k / centi**2 * kilo # in erg/K

def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    globals().update(world)
    for modkey in dir(funcs):
        mod = getattr(funcs, modkey)
        if not hasattr(mod, 'update_func_world'): continue
        mod.update_func_world(mech, world)

def solar_noon_local(LonDegE):
    """
    Assumes that time is Local Time and, therefore, returns 12.
    """
    return 12.

def solar_noon_utc(LonDegE):
    """
    Returns solar noon in UTC based on 15degree timezones 0+-7.5
    LonDegE - degrees longitude (-180, 180)
    """
    _timezone = array([-180, -172.5, -157.5, -142.5, -127.5, -112.5,  -97.5,  -82.5,  -67.5,  -52.5, -37.5,  -22.5,   -7.5,    7.5,   22.5,   37.5,   52.5,   67.5, 82.5,   97.5,  112.5,  127.5,  142.5,  157.5,  172.5,  180]).repeat(2, 0)[1:-1].reshape(-1, 2)
    for i, (low, high) in enumerate(_timezone):
        if LonDegE >= low:
            if LonDegE <= high:
                return 12 -(-12 + i)

solar_noon = solar_noon_local
    
def h2o_from_rh_and_temp(RH, TEMP):
    """
    Return H2O in molecules/cm**3 from RH (0-100) and
    TEMP in K
    """
    TC = TEMP - 273.15
    frh = RH / 100.
    svp_millibar = 6.11 * 10**((7.5 * TC)/(TC+237.3))
    svp_pa = svp_millibar * 100
    vp_pa = svp_pa * frh
    molecule_per_cubic_m = vp_pa * Avogadro / R / TEMP
    molecule_per_cubic_cm = molecule_per_cubic_m * centi**3
    #print RH, TEMP, molecule_per_cubic_cm
    return molecule_per_cubic_cm

def initstdenv(TEMP = 298., P = 101325., ):
    """
    Initialize a std environemnt
    """
    defenv = dict(TEMP = TEMP, P = P)
    Update_M(None, defenv)
    update_func_world(None, defenv)
    del defenv

