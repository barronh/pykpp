__all__ = ['Update_World', 'Update_RATE', 'Update_SUN', 'Update_THETA', 'Update_M', 'add_world_updater', 'solar_declination', 'solar_noon_local', 'solar_noon_utc', \
           'ARR', 'ARR2', 'EP2', 'EP3', 'FALL', 'DP3', 'k_3rd', 'k_arr', \
           'MZ4_TROE', 'MZ4_USR1', 'MZ4_USR2', 'MZ4_USR3', 'MZ4_USR4', 'MZ4_USR5', 'MZ4_USR6', 'MZ4_USR7', 'MZ4_USR8', 'MZ4_USR9', 'MZ4_USR10', 'MZ4_USR11', 'MZ4_USR12', 'MZ4_USR14', 'MZ4_USR21', 'MZ4_USR22', 'MZ4_USR23', 'MZ4_USR24', \
           'RACM_TROE', 'RACM_TROE_EQUIL', 'RACM_THERMAL', 'RACM_THERMAL_T2', \
           'GEOS_STD', 'GEOS_P', 'GEOS_Z', 'GEOS_Y', 'GEOS_X', 'GEOS_C', 'GEOS_K', 'GEOS_V', 'GEOS_E', 'FYRNO3', 'JHNO4_NEAR_IR', 'GEOS_KHO2', 'GEOS_A', 'GEOS_B', 'GEOS_JO3', 'GEOS_G', \
           'CMAQ_1to4', 'CMAQ_5', 'CMAQ_6', 'CMAQ_7', 'CMAQ_8', 'CMAQ_9', 'CMAQ_10', 'CMAQ_10D', 'OH_CO', \
           'TUV_J', 'update_func_world']

from numpy import *
from scipy.constants import *
from tuv.tuv5pt0 import TUV_J
try:
    boltz  = Boltzmann / centi**2 * kilo # in erg/K
except:
    boltz  = constants.k / centi**2 * kilo # in erg/K

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
    _timezone = np.array([-180, -172.5, -157.5, -142.5, -127.5, -112.5,  -97.5,  -82.5,  -67.5,  -52.5, -37.5,  -22.5,   -7.5,    7.5,   22.5,   37.5,   52.5,   67.5, 82.5,   97.5,  112.5,  127.5,  142.5,  157.5,  172.5,  180]).repeat(2, 0)[1:-1].reshape(-1, 2)
    for i, (low, high) in enumerate(_timezone):
        if LonDegE >= low:
            if LonDegE <= high:
                return 12 -(-12 + i)

solar_noon = solar_noon_local
    
def solar_declination(N):
    """
    N - julian day 1-365 (1 = Jan 1; 365 = Dec 31)
    Returns solar declination in radians
    
    wikipedia.org/wiki/Declination_of_the_Sun
    dec_degrees = -23.44 * cos_degrees(360./365 * (N + 10))
    dec_radians = pi / 180. * -23.44) * cos_radians(pi / 180. * 360./365 * (N + 10))
    """
    return -0.40910517666747087 * cos(0.017214206321039961 * (N + 10.))

def Update_SUN(mech, world):
    """
    Updates world dectionary to contain
      SUN - scaling variable between 0 and 1 (following Sandu et al.)

    if t < SunRise or t > SunSet: SUN = 0.
    
    hour = time since noon
    squared = abs(hour) * hour
    
    """
    t = world['t']
    SunRise = world['SunRise']
    SunSet  = world['SunSet']
    Thour = t/3600.0
    Tlocal = Thour % 24.

    if (Tlocal >= SunRise) and (Tlocal <= SunSet):
        Ttmp = (2.0 * Tlocal - SunRise - SunSet) / (SunSet - SunRise)
        if (Ttmp > 0.):
            Ttmp =  Ttmp * Ttmp
        else:
            Ttmp = -Ttmp * Ttmp
        
        SUN = ( 1.0 + cos(pi * Ttmp) )/2.0
    else:
        SUN = 0.0
    world['SUN'] = SUN

def Update_THETA(mech, world):
    """
    Adds solar zenith angle (THETA; angle from solar noon) in degrees
    to the world dictionary based on time
    
    THETA = arccos(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(houra))
    """
    phi = world['Latitude_Radians']
    t = world['t']
    if 'SolarDeclination_Radians' in world:
        dec = world['SolarDeclination_Radians']
    else:
        StartJday = world['StartJday']
        N = StartJday + t // 24
        dec = solar_declination(N)
    Tlocal = (t / 3600.) % 24.
    houra = radians((Tlocal - 12.) * 15.)
    THETA = arccos(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(houra))
    world['THETA'] = degrees(THETA)

def Update_RATE(mech, world):
    """
    Update_RATE is the default rate updater that is called when
    rates are updated for the integrator solution
    
    1) Call update_func_world(world)
    2) Set world['rate_const'] equal to evaluated world['rate_const_exp']
    """
    mech.rate_const = eval(mech.rate_const_exp, None, world)

def Update_M(mech, world):
    """
    Adds M, O2, N2, and H2 to world based on Pressure (P=Pascals; R=m**3 * Pascals/K/mol
    """
    try:
        M = float(world['M'])
    except:
        M = eval('P / (R / centi**3) / TEMP * Avogadro', None, world)
    world['M'] = M
    world['O2'] = 0.20946 * M
    world['N2'] = 0.78084 * M
    world['H2'] = 0.00000055 * M
    

def Update_World(mech, world):
    """
    Calls 
        + user added functions (added using add_world_updater)
        update_func_world(world)
        Update_RATE(mech, world)
    
    *Note: Update_SUN, Update_THETA, and Update_M can be added by
           the mechanism heuristically.
    
    see these functions for more details
    """
    for func in Update_World.updaters:
        func(mech, world)

    update_func_world(world)
    Update_RATE(mech, world)

Update_World.updaters = []

def add_world_updater(func):
    """
    Add func to be called with mech and world
    to update the world environment
    """
    Update_World.updaters.append(func)

def ARR( A0, B0, C0 ):
    """
    A0, B0 and C0 - numeric values used to
    calculate a reaction rate (1/s) based on 
    the Arrhenius equation in the following
    form:
    
    A0 * exp(-B0/TEMP) * (TEMP / 300.)**(C0)
    
    Returns a rate in per time
    """
    out = A0 * exp(-B0 / TEMP) * (TEMP / 300.0)**(C0)
    return out

def ARR2( A0, B0):
    """
    A0 and B0 - numeric values used to
    calculate a reaction rate (1/s) based on 
    the Arrhenius equation in the following
    form:
    
    A0 * exp(-B0/TEMP)
    
    Returns a rate in per time
    """
    out = A0 * exp( B0 / TEMP )
    return out

def EP2(A0, C0, A2, C2, A3, C3):
    """
    A0, C0, A2, C2, A3, and C3 - numeric 
    values used to calculate 3 rates (K0, K2,
    and K3), each of the form A * exp(-C / TEMP),
    to return a rate of the following form:
    
    K0 + K3 * M * 1e6 / (1. + K3 * M * 1e6 / K2)

    Returns a rate in per time
    """
    K0 = A0 * exp(-C0 / TEMP)
    K2 = A2 * exp(-C2 / TEMP)
    K3 = A3 * exp(-C3 / TEMP)
    K3 = K3 * M * 1.0E6
    return K0 + K3 / (1.0 + K3 / K2 )

def EP3(A1, C1, A2, C2):
    """
    A1, C1, A2, and C2 - numeric 
    values used to calculate 3 rates (K0, K2,
    and K3), each of the form A * exp(-C / TEMP),
    to return a rate of the following form:
    
    K1 + K2 * M * 1e6

    Returns a rate in per time
    """
    K1 = A1 * exp(-C1 / TEMP)
    K2 = A2 * exp(-C2 / TEMP)
    return K1 + K2 * (1.0E6 * M)

def FALL ( A0, B0, C0, A1, B1, C1, CF):
    """
    Troe fall off equation
    
    A0, B0, C0, A1, B1, C1 - numeric values
    to calculate 2 reaction rates (K0, K1) using ARR
    function; returns a rate in the following form
    
    K0M = K0 * M * 1e6
    KR = K0 / K1
    
    Returns (K0M / (1.0 + KR))* CF**(1.0 / (1.0 + (log10(KR))**2))
    """
    
    K0 = ARR(A0, B0, C0)
    K1 = ARR(A1, B1, C1)
    K0 = K0 * M * 1.0E6
    K1 = K0 / K1
    return (K0 / (1.0 + K1))* CF**(1.0 / (1.0 + (log10(K1))**2))

def k_3rd(temp, cair, k0_300K, n, kinf_300K, m, fc):
    """
    """
    zt_help = 300. / temp
    k0_T    = k0_300K   * zt_help**(n) * cair # k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        # k_inf at current T
    k_ratio = k0_T / kinf_T
    return k0_T / (1. + k_ratio) * fc**(1. / (1. + log10(k_ratio)**2))

  
DP3 = EP3

def k_arr(k_298, tdep, temp):
    """
    """
    return k_298 * exp(tdep * (1. / temp - 3.3540E-3)) # 1/298.15=3.3540e-3

def MZ4_TROE( A0, B0, A1, B1, factor):
  """
  Troe fall off equation as calculated in
  MOZART4
  """
  #REAL(kind=dp) A0, B0, factor, A1, B1
  #REAL(kind=dp) ko, kinf, xpo
  ko = A0 * (300.e0 / TEMP)**B0
  kinf = A1 * (300.e0 / TEMP)**B1
  xpo  = ko * M / kinf
  MZ4_TROE_TMP = ko / (1. + xpo)
  xpo  = log10( xpo )
  xpo  = 1. / (1. + xpo*xpo)
  return MZ4_TROE_TMP * factor**xpo


def MZ4_USR1():
    """
    USR1 reaction rate as defined in MOZART4
    
    Returns 6.e-34 * (300.e0/TEMP)**2.4
    """
    return 6.e-34 * (300.e0/TEMP)**2.4


def MZ4_USR2():
    """
    USR2 reaction rate as defined in MOZART4
    
    Returns MZ4_TROE(8.5e-29, 6.5e0, 1.1e-11, 1.e0, .6e0)
    """
    return MZ4_TROE(8.5e-29, 6.5e0, 1.1e-11, 1.e0, .6e0)

def MZ4_USR3():
    """
    USR3 reaction rate as defined in MOZART4

    Returns MZ4_USR2() * 3.333e26 * exp( -10990.e0/TEMP )
    """
    return MZ4_USR2() * 3.333e26 * exp( -10990.e0/TEMP )

def MZ4_USR4():
    """
    USR4 reaction rate as defined in MOZART4
    
    Returns MZ4_TROE(2.0e-30, 3.0e0, 2.5e-11, 0.e0, .6e0)
    """
    return MZ4_TROE(2.0e-30, 3.0e0, 2.5e-11, 0.e0, .6e0)

def MZ4_USR5():
    """
    USR5 reaction rate as defined in MOZART4
    
    TINV = 1/TEMP
    ko = M * 6.5e-34 * exp( 1335.*tinv )
    ko = ko / (1. + ko/(2.7e-17*exp( 2199.*tinv )))

    Returns ko + 2.4e-14*exp( 460.*tinv )
    
    """
    #REAL(kind=dp) KO, TINV

    TINV = 1/TEMP
    ko = M * 6.5e-34 * exp( 1335.*tinv )
    ko = ko / (1. + ko/(2.7e-17*exp( 2199.*tinv )))
    return ko + 2.4e-14*exp( 460.*tinv )

def MZ4_USR6():
    """
    USR6 reaction rate as defined in MOZART4
    
    Returns MZ4_TROE(1.8e-31, 3.2e0, 4.7e-12, 1.4e0, .6e0)
    """
    return MZ4_TROE(1.8e-31, 3.2e0, 4.7e-12, 1.4e0, .6e0)

def MZ4_USR7():
    """
    USR7 reaction rate as defined in MOZART4
    
    Returns MZ4_USR6() * exp( -10900./TEMP )/ 2.1e-27
    """
    return MZ4_USR6() * exp( -10900./TEMP )/ 2.1e-27

def MZ4_USR8():
    """
    USR8 reaction rate as defined in MOZART4
    
    Returns 1.5e-13 * (1. + 6.e-7 * boltz * M * TEMP)
    """
    #real, parameter ::  boltz    = 1.38044e-16      ! erg/k
    return 1.5e-13 * (1. + 6.e-7 * boltz * M * TEMP)

def MZ4_USR9():
    """
    USR9 reaction rate as defined in MOZART4
    
    REAL(kind = dp) ko, kinf, fc, tinv
    tinv = 1.0/TEMP
    ko   = 2.3e-13 * exp( 600.0*tinv )
    kinf = 1.7e-33 * M * exp( 1000.0*tinv )
    fc   = 1.0 + 1.4e-21 * H2O * exp( 2200.0*tinv )

    Returns (ko + kinf) * fc
    """
    #REAL(kind = dp) ko, kinf, fc, tinv
    tinv = 1.0/TEMP
    ko   = 2.3e-13 * exp( 600.0*tinv )
    kinf = 1.7e-33 * M * exp( 1000.0*tinv )
    fc   = 1.0 + 1.4e-21 * H2O * exp( 2200.0*tinv )

    return (ko + kinf) * fc

def MZ4_USR10():
    """
    USR10 reaction rate as defined in MOZART4
    
    Returns MZ4_TROE(8.e-27, 3.5e0, 3.e-11, 0e0, .5e0)
    """
    return MZ4_TROE(8.e-27, 3.5e0, 3.e-11, 0e0, .5e0)

def MZ4_USR11():
    """
    USR11 reaction rate as defined in MOZART4
    
    Returns MZ4_TROE(8.5e-29, 6.5e0, 1.1e-11, 1.0, .6e0)
    """
    return MZ4_TROE(8.5e-29, 6.5e0, 1.1e-11, 1.0, .6e0)

def MZ4_USR12():
    """
    USR12 reaction rate as defined in MOZART4
    
    Return MZ4_USR11() * 1.111e28 * exp( -14000.0 / TEMP )
    """
    return MZ4_USR11() * 1.111e28 * exp( -14000.0 / TEMP )

def MZ4_USR14():
    """
    USR14 reaction rate as defined in MOZART4
    
    Return 1.1e-11 * 300.e0/ TEMP / M
    """
    return 1.1e-11 * 300.e0/ TEMP / M

def MZ4_USR15():
    """
    USR15 reaction rate as defined in MOZART4
    
    Returns MZ4_USR14() * 1.111e28 *  exp( -14000.e0 / TEMP )
    """
    return MZ4_USR14() * 1.111e28 *  exp( -14000.e0 / TEMP )

def MZ4_USR21():
    """
    USR21 reaction rate as defined in MOZART4
    
    Returns TEMP**2 * 7.69e-17 * exp( 253.e0/TEMP )
    """
    return TEMP**2 * 7.69e-17 * exp( 253.e0/TEMP )

def MZ4_USR22():
    """
    USR22 reaction rate as defined in MOZART4

    Returns 3.82e-11 * exp( -2000.0/TEMP ) + 1.33e-13
    """
    return 3.82e-11 * exp( -2000.0/TEMP ) + 1.33e-13

def MZ4_USR23():
    """
    USR23 reaction rate as defined in MOZART4
        
    fc = 3.0e-31 *(300.0/TEMP)**3.3e0
    ko = fc * M / (1.0 + fc * M / 1.5e-12) 
    Returns ko * .6e0**(1. + (log10(fc * M / 1.5e-12))**2.0)**(-1.0)
    """
    #REAL(kind=dp) ko, fc
    
    fc = 3.0e-31 *(300.0/TEMP)**3.3e0
    ko = fc * M / (1.0 + fc * M / 1.5e-12) 
    return ko * .6e0**(1. + (log10(fc * M / 1.5e-12))**2.0)**(-1.0)

def MZ4_USR24():
    """
    USR24 reaction rate as defined in MOZART4
    
    #REAL(kind=dp) ko
    ko = 1.0 + 5.5e-31 * exp( 7460.0/TEMP ) * M * 0.21e0
    return 1.7e-42 * exp( 7810.0/TEMP ) * M * 0.21e0 / ko
    """

    #REAL(kind=dp) ko
    ko = 1.0 + 5.5e-31 * exp( 7460.0/TEMP ) * M * 0.21e0
    return 1.7e-42 * exp( 7810.0/TEMP ) * M * 0.21e0 / ko

def RACM_TROE(A0, B0, A1, B1):
    """
    RACM_TROE equation as defined in the RACM SBOX model

    K0 = (A0 * (TEMP/300.0)**(-B0))
    K1 = (A1 * (TEMP/300.0)**(-B1))
    K0 = K0 * M
    K1 = K0 / K1

    Returns (K0 / (1.0 + K1))*   \
       CF**(1.0 / (1.0 / N + (log10(K1))**2))
    """
    #REAL A0, B0, A1, B1, C1
    #REAL(kind=dp) K0, K1
    #REAL(kind=dp), PARAMETER :: CF = 0.6_dp, N = 1._dp
    K0 = (A0 * (TEMP/300.0)**(-B0))
    K1 = (A1 * (TEMP/300.0)**(-B1))
    K0 = K0 * M
    K1 = K0 / K1
    return (K0 / (1.0 + K1))*   \
       CF**(1.0 / (1.0 / N + (log10(K1))**2))
       
def RACM_TROE_EQUIL(A0, B0, A1, B1, A2, C2):
    """
    RACM Troe equilibrium equation as defined in the RACM SBOX model

    #REAL A0, B0, A1, B1, A2, C2
    return RACM_TROE( A0, B0, A1, B1) *  (1./A2 * exp(-C2 / TEMP))
    """
    #REAL A0, B0, A1, B1, A2, C2
    return RACM_TROE( A0, B0, A1, B1) *  (1./A2 * exp(-C2 / TEMP))

def RACM_THERMAL(A0, B0):
    """
    RACM Thermal equation as defined in the RACM SBOX model

    #REAL A0, B0
    #   RACM2 reaction rates have the form K = A * EXP(-B / T)
    #   
    #   Translation adds a 0 C
    return (A0 * exp(-B0 / TEMP))
    """
    #REAL A0, B0
    #   RACM2 reaction rates have the form K = A * EXP(-B / T)
    #   
    #   Translation adds a 0 C
    return (A0 * exp(-B0 / TEMP))

def RACM_THERMAL_T2(A0, B0):
    """
    RACM Thermal T2 equation as defined in the RACM SBOX model

    #REAL A0, B0
    #REAL, PARAMETER :: C0 = 0.
    #   
    #   Translation adds a 0 C
    return (A0)*TEMP**2*exp(-(B0)/TEMP)
    """
    #REAL A0, B0
    #REAL, PARAMETER :: C0 = 0.
    #   
    #   Translation adds a 0 C
    return (A0)*TEMP**2*exp(-(B0)/TEMP)

def GEOS_STD(A0, B0, C0):
    """
    GEOS-Chem standard reaction rate with 
    the form K = A * (300 / T)**B * EXP(C / T)

    Returns A0 * (300. / TEMP)**B0 * exp(C0 / TEMP)
    """
    #REAL A0, B0, C0

    return A0 * (300. / TEMP)**B0 * exp(C0 / TEMP)
  
def GEOS_P(A0, B0, C0, A1, B1, C1, \
                                FCV, FCT1, FCT2):
    """
    GEOS-Chem pressure dependent TROE falloff equation
    
    if (FCT2 != 0.000000e+00):
      CF = exp(-TEMP / FCT1) + exp(-FCT2 / TEMP) 
    elif (FCT1 != 0.000000e+00):
      CF = exp(-TEMP / FCT1)
    else:
      CF = FCV

    K0M = GEOS_STD(A0, B0, C0) * M

    K1 = GEOS_STD(A1, B1, C1)
    K1 = K0M / K1

    return (K0M / (1.0 + K1))*   \
           (CF)**(1.0 / (1.0 + (log10(K1))**2))

    """
    
    #REAL A0, B0, C0, A1, B1, C1 ,CF
    #REAL FCV, FCT1, FCT2
    #REAL(kind=dp) K0M, K1

    if (FCT2 != 0.000000e+00):
      CF = exp(-TEMP / FCT1) + exp(-FCT2 / TEMP) 
    elif (FCT1 != 0.000000e+00):
      CF = exp(-TEMP / FCT1)
    else:
      CF = FCV
    

    K0M = GEOS_STD(A0, B0, C0) * M

    K1 = GEOS_STD(A1, B1, C1)
    K1 = K0M / K1

    return (K0M / (1.0 + K1))*   \
           (CF)**(1.0 / (1.0 + (log10(K1))**2))
  
def GEOS_Z(A0, B0, C0, A1, B1, C1, A2, B2, C2):
    """
    GEOS-Chem Z reaction rate form

    K0 = GEOS_STD(A0, B0, C0)
    K1 = GEOS_STD(A1, B1, C1)*M
    K2 = GEOS_STD(A2, B2, C2)

    Returns (K0 + K1) * (1 + H2O * K2)
    """
    #REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    #REAL(kind=dp) K0, K1, K2

    K0 = GEOS_STD(A0, B0, C0)
    K1 = GEOS_STD(A1, B1, C1)*M
    K2 = GEOS_STD(A2, B2, C2)

    return (K0 + K1) * (1 + H2O * K2)

def GEOS_Y(A0, B0, C0):
    """
    GEOS-Chem reaction form Y
    
    A0, B0, and C0 are numeric inputs that are ignored
    IGNORES INPUTS per v08-02-04 update
    """
    #REAL A0, B0, C0
    #REAL(kind=dp) K0
    #REAL(kind=dp) KHI1,KLO1,XYRAT1,BLOG1,FEXP1,KHI2,KLO2,XYRAT2,BLOG2
    #REAL(kind=dp) FEXP2,KCO1,KCO2,KCO

    #IGNORES INPUTS per v08-02-04 update
    #K0 = GEOS_STD(A0, B0, C0)
    #GEOS_Y = K0 * (1 + .6 * (PRESS * 100.) / 101325.)
    KLO1 = 5.9e-33 * (300 / TEMP)**(1.4e0) 
    KHI1 = 1.1e-12 * (300 / TEMP)**(-1.3e0)
    XYRAT1 = KLO1 * M /KHI1
    BLOG1 = log10(XYRAT1)
    FEXP1 = 1.e0 / (1.e0 + BLOG1 * BLOG1)
    KCO1 = KLO1 * M * 0.6**FEXP1 /(1.e0 + XYRAT1)
    KLO2 = 1.5e-13 * (300 / TEMP)**(-0.6e0)
    KHI2 = 2.1e09 * (300 / TEMP)**(-6.1e0)
    XYRAT2 = KLO2 * M /KHI2
    BLOG2 = log10(XYRAT2)
    FEXP2 = 1.e0 / (1.e0 + BLOG2 * BLOG2)
    KCO2 = KLO2 * 0.6**FEXP2 / (1.e0 + XYRAT2)
    return KCO1 + KCO2

def GEOS_X(A0, B0, C0, A1, B1, C1, A2, B2, C2):
    """
    GEOS-Chem reaction rate form Z

    K0 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    K3 = GEOS_STD(A2, B2, C2)
    K3 = K3 * M

    Returns K0 + K3 / (1.0 + K3 / K2 )
    """
    #REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    #REAL(kind=dp) K0, K2, K3            
    K0 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    K3 = GEOS_STD(A2, B2, C2)
    K3 = K3 * M
    return K0 + K3 / (1.0 + K3 / K2 )

def GEOS_C(A0, B0, C0):
    """
    GEOS-Chem reaction form C

    K1 = GEOS_STD(A0, B0, C0)
    
    Returns K1 * (O2 + 3.5e18) / (2.0 * O2 + 3.5e18)
    """
    #REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    #REAL(kind=dp) K1
    K1 = GEOS_STD(A0, B0, C0)
    return K1 * (O2 + 3.5e18) / (2.0 * O2 + 3.5e18)

def GEOS_K(A0, B0, C0):
    """
    GEOS-Chem reaction form K
    
    ** Not implemented returns 0.
    """
    warn("GEOS_K not implemented, returning 0")
    #REAL A0, B0, C0
    return 0

def GEOS_V(A0, B0, C0, A1, B1, C1):
    """
    GEOS-Chem reaction form V
    
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    return K1 / (1 + K2)
    """
    #REAL A0, B0, C0, A1, B1, C1
    #REAL(kind=dp) K1, K2
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    return K1 / (1 + K2)

def GEOS_E(A0, B0, C0, Kf):
    """
    GEOS-Chem reaction form E
    
    K1 = GEOS_STD(A0, B0, C0)

    Returns Kf / K1
    """

    #REAL A0, B0, C0
    #REAL(kind=dp) K1, Kf
    K1 = GEOS_STD(A0, B0, C0)
    return Kf / K1

def FYRNO3(CN):
    """
    GEOS-Chem equation FYRNO3 implemented based on GEOS-Chem
    version 9
    """
    
    Y300 = .826; ALPHA = 1.94E-22
    BETA = .97; XM0 = 0.; XMINF = 8.1
    XF = .411
    
    #REAL*4 CN
    #REAL*4 XCARBN, ZDNUM, TT, XXYN, YYYN, AAA, ZZYN, RARB
    XCARBN = CN
    ZDNUM = M
    TT = TEMP
    
    XXYN = ALPHA * exp(BETA * XCARBN) * ZDNUM * ((300. / TT)**XM0)
    YYYN = Y300 * ((300. / TT)**XMINF)
    AAA = log10(XXYN / YYYN)
    ZZYN = 1. / (1. + AAA / AAA)
    RARB = (XXYN / (1. + (XXYN / YYYN))) * (XF**ZZYN)
    return RARB / (1. + RARB)

def JHNO4_NEAR_IR(HNO4J):
    """
    Adding 1e-5 (1/s) to HNO4 photolysis to
    account for near IR
    """
    if (HNO4J > 0.e0):
      return HNO4J + 1e-5
    else:
        return HNO4J


def GEOS_KHO2(A0, B0, C0):
    """
    Implemented KHO2 based on GEOS-Chem version 9
    """
    STKCF = HO2( SLF_RAD, TEMP, M, SQRT(DBLE(A0)), C(ind_HO2), 8, 0 )
    GEOS_KHO2 = ARSL1K( SLF_AREA, SLF_RAD, M, STKCF, sqrt(TEMP), sqrt((A0)))
    STKCF = HO2( EPOM_RAD, TEMP, M, sqrt((A0)), C(ind_HO2), 10, 0 )
    return GEOS_KHO2 + ARSL1K( EPOM_AREA, EPOM_RAD, M, STKCF, sqrt(TEMP), sqrt(DBLE(A0)))

def GEOS_A(A0, B0, C0, A1, B1, C1 ):
    """
    GEOS-Chem reaction form A
    
    TMP_A0 = A0 * FYRNO3(A1)

    Returns GEOS_STD(TMP_A0, B0, C0)
    """
    #REAL A0, B0, C0, A1, B1, C1
    #REAL TMP_A0
    TMP_A0 = A0 * FYRNO3(A1)
    return GEOS_STD(TMP_A0, B0, C0)

def GEOS_B(A0, B0, C0, A1, B1, C1 ):
    """
    GEOS-Chem reaction form B
    
    TMP_A0 = A0 * ( 1. - FYRNO3(A1) )

    Returns GEOS_STD(TMP_A0, B0, C0)
    """
    #REAL A0, B0, C0, A1, B1, C1
    #REAL TMP_A0
    TMP_A0 = A0 * ( 1. - FYRNO3(A1) )
    return GEOS_STD(TMP_A0, B0, C0)

def GEOS_JO3(O3J):
    """
    GEOS-Chem reaction form ozone photolysis
    
    T3I = 1.0/TEMP
    Returs O3J * \
               1.45e-10 * exp( 89.0 * T3I) * H2O / \
               ( 1.45e-10 * exp( 89.0 * T3I) * H2O + \
                 2.14e-11 * exp(110.0 * T3I) * N2 + \
                 3.20e-11 * exp( 70.0 * T3I) * O2 \
               )
    """
    
    #REAL(kind=dp) O3J, T3I
    T3I = 1.0/TEMP
    return O3J * \
               1.45e-10 * exp( 89.0 * T3I) * H2O / \
               ( 1.45e-10 * exp( 89.0 * T3I) * H2O + \
                 2.14e-11 * exp(110.0 * T3I) * N2 + \
                 3.20e-11 * exp( 70.0 * T3I) * O2 \
               )

def GEOS_G(A0, B0, C0, A1, B1, C1):
    """
    GEOS-Chem reaction form A
    
    
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    
    Returns K1 / ( 1.0 + K1 * O2 )
    """
    #REAL A0, B0, C0, A1, B1, C1
    #REAL(kind=dp) K1, K2
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    return K1 / ( 1.0 + K1 * O2 )

def CMAQ_1to4(A0, B0, C0):
    """
    CMAQ reaction rates form 1-4 have the form K = A * (T/300.0)**B * EXP(-C/T) 
    """
    #REAL A0, B0, C0

    return (A0 * (TEMP/300.0)**B0 * exp(-C0/TEMP))

def CMAQ_5(A0, B0, C0, Kf):
    """
    CMAQ reaction form 5

    K1 = CMAQ_1to4(A0, B0, C0)

    Returns Kf / K1
    """
    #REAL A0, B0, C0
    #REAL(kind=dp) K1, Kf
    K1 = CMAQ_1to4(A0, B0, C0)
    return Kf / K1
  
def CMAQ_6(A0, B0, C0, Kf):
    """
    CMAQ reaction form 6
    
    K1 = CMAQ_1to4(A0, B0, C0)

    Returns Kf * K1
    """
    #REAL A0, B0, C0
    #REAL(kind=dp) K1, Kf
    K1 = CMAQ_1to4(A0, B0, C0)
    return Kf * K1
  
def CMAQ_7(A0, B0, C0):
    """
    CMAQ reaction form 6
    
    K0 = CMAQ_1to4(A0, B0, C0)
    
    Returns K0 * (1 + .6 * PRESS / 101325.) # Pressure is in Pascals
    """
    #REAL A0, B0, C0
    #REAL(kind=dp) K0
    K0 = CMAQ_1to4(A0, B0, C0)
    return K0 * (1 + .6 * PRESS / 101325.) # Pressure is in Pascals
  
def CMAQ_8(A0, C0, A2, C2, A3, C3):
    """
    CMAQ reaction form 8
    
    K0 = (A0) * exp(-(C0) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)
    K3 = (A3) * exp(-(C3) / TEMP)
    K3 = K3 * M

    Returns K0 + K3 / (1.0 + K3 / K2 )
    """
    #REAL A0, C0, A2, C2, A3, C3
    #REAL(kind=dp) K0, K2, K3            
    K0 = (A0) * exp(-(C0) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)
    K3 = (A3) * exp(-(C3) / TEMP)
    K3 = K3 * M
    return K0 + K3 / (1.0 + K3 / K2 )

def CMAQ_9(A1, C1, A2, C2):
    """
    CMAQ reaction rate form 9
    
    K1 = (A1) * exp(-(C1) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)
    return K1 + K2 * M
    """    
    #REAL(kind=dp) A1, C1, A2, C2
    #REAL(kind=dp) K1, K2      
    K1 = (A1) * exp(-(C1) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)
    return K1 + K2 * M

def CMAQ_10 ( A0, B0, C0, A1, B1, C1, CF, N):
    """
    CMAQ reaction rate form 10

    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0 * M
    K1 = K0 / K1

    Returns (K0 / (1.0 + K1))*   \
         (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))
    """
    #REAL A0, B0, C0, A1, B1, C1, CF, N
    #REAL(kind=dp) K0, K1     
    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0 * M
    K1 = K0 / K1
    return (K0 / (1.0 + K1))*   \
         (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))

def CMAQ_10D ( A0, B0, C0, A1, B1, C1, CF, N):
    """
    Same as reaction rate form 10, but implemented
    to provide compatibility for fortran code
    that need a DOUBLE form
    """
    return CMAQ_10(A0, B0, C0, A1, B1, C1, CF, N)

def OH_CO ( A0, B0, C0, A1, B1, C1, CF, N):
    """
    OH + CO reaction rate
    
    *Note: Mostly like CMAQ_10, but slight difference in K1
    
    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0
    K1 = K0 / (K1 / M)
    return (K0 / (1.0 + K1))*   \
         (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))

    """
    #REAL A0, B0, C0, A1, B1, C1, CF, N
    #REAL(kind=dp) K0, K1     
    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0
    K1 = K0 / (K1 / M)
    return (K0 / (1.0 + K1))*   \
         (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))


def update_func_world(world):
    """
    Function to update globals for user defined functions
    """
    globals().update(world)