__all__ = [
    'GEOS_STD', 'GEOS_P', 'GEOS_Z', 'GEOS_Y', 'GEOS_X', 'GEOS_C', 'GEOS_K',
    'GEOS_HR', 'GEOS_V', 'GEOS_E', 'FYRNO3', 'JHNO4_NEAR_IR', 'GEOS_KHO2',
    'GEOS_A', 'GEOS_B', 'GEOS_JO3', 'GEOS_G', 'GEOS_T', 'GEOS_F', 'GEOS_L',
    'GEOS_O', 'GEOS_N', 'GEOS_Q', 'GEOS_JO3', 'GEOS_JO3_2', 'update_func_world'
]
from warnings import warn
from numpy import exp, log10, float64, sqrt

DBLE = float64
N2 = O2 = M = TEMP = INVTEMP = 0
int_HO2 = C = H2O = H2 = 0


def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    global INVTEMP, T300I
    globals().update(world)
    INVTEMP = 1 / world['TEMP']
    T300I = 300. / world['TEMP']


def GEOS_STD(A0, B0, C0):
    """
    GEOS-Chem standard reaction rate with
    the form K = A * (300 / T)**B * EXP(C / T)

    Returns A0 * (300. / TEMP)**B0 * exp(C0 / TEMP)
    """
    global INVTEMP
    # REAL A0, B0, C0

    return A0 * (T300I)**B0 * exp(C0 * INVTEMP)


def GEOS_P(A0, B0, C0, A1, B1, C1,
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

    # REAL A0, B0, C0, A1, B1, C1 ,CF
    # REAL FCV, FCT1, FCT2
    # REAL(kind=dp) K0M, K1

    if (FCT2 != 0.000000e+00):
        CF = exp(-TEMP / FCT1) + exp(-FCT2 * INVTEMP)
    elif (FCT1 != 0.000000e+00):
        CF = exp(-TEMP / FCT1)
    else:
        CF = FCV

    K0M = GEOS_STD(A0, B0, C0) * M

    K1 = GEOS_STD(A1, B1, C1)
    K1 = K0M / K1

    return (K0M / (1.0 + K1)) * \
           (CF)**(1.0 / (1.0 + (log10(K1))**2))


def GEOS_Z(A0, B0, C0, A1, B1, C1, A2, B2, C2):
    """
    GEOS-Chem Z reaction rate form

    K0 = GEOS_STD(A0, B0, C0)
    K1 = GEOS_STD(A1, B1, C1)*M
    K2 = GEOS_STD(A2, B2, C2)

    Returns (K0 + K1) * (1 + H2O * K2)
    """
    # REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    # REAL(kind=dp) K0, K1, K2

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
    # REAL A0, B0, C0
    # REAL(kind=dp) K0
    # REAL(kind=dp) KHI1,KLO1,XYRAT1,BLOG1,FEXP1,KHI2,KLO2,XYRAT2,BLOG2
    # REAL(kind=dp) FEXP2,KCO1,KCO2,KCO

    # IGNORES INPUTS per v08-02-04 update
    # K0 = GEOS_STD(A0, B0, C0)
    # GEOS_Y = K0 * (1 + .6 * (PRESS * 100.) / 101325.)
    KLO1 = 5.9e-33 * (T300I)**(1.4e0)
    KHI1 = 1.1e-12 * (T300I)**(-1.3e0)
    XYRAT1 = KLO1 * M / KHI1
    BLOG1 = log10(XYRAT1)
    FEXP1 = 1.e0 / (1.e0 + BLOG1 * BLOG1)
    KCO1 = KLO1 * M * 0.6**FEXP1 / (1.e0 + XYRAT1)
    KLO2 = 1.5e-13 * (T300I)**(-0.6e0)
    KHI2 = 2.1e09 * (T300I)**(-6.1e0)
    XYRAT2 = KLO2 * M / KHI2
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
    # REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    # REAL(kind=dp) K0, K2, K3
    K0 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    K3 = GEOS_STD(A2, B2, C2)
    K3 = K3 * M
    return K0 + K3 / (1.0 + K3 / K2)


def GEOS_C(A0, B0, C0):
    """
    GEOS-Chem reaction form C

    K1 = GEOS_STD(A0, B0, C0)

    Returns K1 * (O2 + 3.5e18) / (2.0 * O2 + 3.5e18)
    """
    # REAL A0, B0, C0, A1, B1, C1, A2, B2, C2
    # REAL(kind=dp) K1
    K1 = GEOS_STD(A0, B0, C0)
    return K1 * (O2 + 3.5e18) / (2.0 * O2 + 3.5e18)


def GEOS_K(A0, B0, C0):
    """
    GEOS-Chem reaction form K

    ** Not implemented returns 0.
    """
    warn("GEOS_K not implemented, returning 0")
    # REAL A0, B0, C0
    return 0


def GEOS_HR(A0, B0, C0, A1, B1, C1):
    """
    GEOS-Chem reaction form HR

    ** Not implemented returns 0.
    """
    # REAL A0, B0, C0
    XCARBN = A1
    return GEOS_STD(A0, B0, C0) * (1e+0 - exp(-0.245e+0 * XCARBN))


def GEOS_N(A0, B0, C0):
    """
    GEOS-Chem reaction form N

    ** Not implemented returns 0.
    """
    ONE_SEVENTYTHREE = 1e+0 / 73e+0
    GLYC_FRAC = 1e+0 - 11.0729e+0 * exp(-ONE_SEVENTYTHREE * TEMP)
    if GLYC_FRAC < 0:
        GLYC_FRAC = 0
    # REAL A0, B0, C0
    return GEOS_STD(A0, B0, C0) * GLYC_FRAC


def GEOS_O(A0, B0, C0):
    """
    GEOS-Chem reaction form O

    ** Not implemented returns 0.
    """
    ONE_SEVENTYTHREE = 1e+0 / 73e+0
    GLYC_FRAC = 1e+0 - 11.0729e+0 * exp(-ONE_SEVENTYTHREE * TEMP)
    if GLYC_FRAC < 0:
        GLYC_FRAC = 0
    # REAL A0, B0, C0
    return GEOS_STD(A0, B0, C0) * (1 - GLYC_FRAC)


def GEOS_Q(A0):
    """
    GEOS-Chem reaction form Q

    ** Not implemented returns 0.
    """
    RO1DplH2O = 1.63e-10 * exp(60. / TEMP) * H2O
    RO1DplH2 = 1.2e-10 * H2
    RO1DplN2 = 2.15e-11 * exp(110. / TEMP) * N2
    RO1DplO2 = 3.30e-11 * exp(55. / TEMP) * O2
    RO1D = RO1DplH2O + RO1DplH2 + RO1DplN2 + RO1DplO2
    return A0 * RO1DplH2O / RO1D


def GEOS_F(A0, B0, C0):
    ONE_SIXTY = 1. / 60.
    HAC_FRAC = 1e+0 - 23.7e+0 * exp(-ONE_SIXTY * TEMP)
    if HAC_FRAC < 0e+0:
        HAC_FRAC = 0e+0
    return GEOS_STD(A0, B0, C0) * HAC_FRAC


def GEOS_L(A0, B0, C0):
    ONE_SIXTY = 1. / 60.
    HAC_FRAC = 1e+0 - 23.7e+0 * exp(-ONE_SIXTY * TEMP)
    if HAC_FRAC < 0e+0:
        HAC_FRAC = 0e+0
    return GEOS_STD(A0, B0, C0) * (1 - HAC_FRAC)


def GEOS_T(A0, B0, C0):
    """
    GEOS-Chem reaction form T

    ** Not implemented returns 0.
    """
    warn("GEOS_T assumed active.")
    # REAL A0, B0, C0
    return GEOS_STD(A0, B0, C0)


def GEOS_V(A0, B0, C0, A1, B1, C1):
    """
    GEOS-Chem reaction form V

    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    return K1 / (1 + K2)
    """
    # REAL A0, B0, C0, A1, B1, C1
    # REAL(kind=dp) K1, K2
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    return K1 / (1 + K2)


def GEOS_E(A0, B0, C0, Kf):
    """
    GEOS-Chem reaction form E

    K1 = GEOS_STD(A0, B0, C0)

    Returns Kf / K1
    """

    # REAL A0, B0, C0
    # REAL(kind=dp) K1, Kf
    K1 = GEOS_STD(A0, B0, C0)
    return Kf / K1


def FYRNO3(CN):
    """
    GEOS-Chem equation FYRNO3 implemented based on GEOS-Chem
    version 9
    """

    Y300 = .826
    ALPHA = 1.94E-22
    BETA = .97
    XM0 = 0.
    XMINF = 8.1
    XF = .411

    # REAL*4 CN
    # REAL*4 XCARBN, ZDNUM, TT, XXYN, YYYN, AAA, ZZYN, RARB
    XCARBN = CN
    ZDNUM = M
    TT = TEMP

    # XXYN = ALPHA * exp(BETA * XCARBN) * ZDNUM * ((300. / TT)**XM0)
    # YYYN = Y300 * ((300. / TT)**XMINF)
    XXYN = ALPHA * exp(BETA * XCARBN) * ZDNUM * ((T300I)**XM0)
    YYYN = Y300 * ((T300I)**XMINF)
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
    STKCF = HO2(SLF_RAD, TEMP, M, sqrt(DBLE(A0)), C(ind_HO2), 8, 0)
    GEOS_KHO2 = ARSL1K(SLF_AREA, SLF_RAD, M, STKCF, sqrt(TEMP), sqrt((A0)))
    STKCF = HO2(EPOM_RAD, TEMP, M, sqrt((A0)), C(ind_HO2), 10, 0)
    out = (
        GEOS_KHO2 + ARSL1K(
            EPOM_AREA, EPOM_RAD, M, STKCF, sqrt(TEMP), sqrt(DBLE(A0))
        )
    )
    return out


def GEOS_A(A0, B0, C0, A1, B1, C1):
    """
    GEOS-Chem reaction form A

    TMP_A0 = A0 * FYRNO3(A1)

    Returns GEOS_STD(TMP_A0, B0, C0)
    """
    # REAL A0, B0, C0, A1, B1, C1
    # REAL TMP_A0
    TMP_A0 = A0 * FYRNO3(A1)
    return GEOS_STD(TMP_A0, B0, C0)


def GEOS_B(A0, B0, C0, A1, B1, C1):
    """
    GEOS-Chem reaction form B

    TMP_A0 = A0 * ( 1. - FYRNO3(A1) )

    Returns GEOS_STD(TMP_A0, B0, C0)
    """
    # REAL A0, B0, C0, A1, B1, C1
    # REAL TMP_A0
    TMP_A0 = A0 * (1. - FYRNO3(A1))
    return GEOS_STD(TMP_A0, B0, C0)


def GEOS_G(A0, B0, C0, A1, B1, C1):
    """
    GEOS-Chem reaction form A

    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)

    Returns K1 / ( 1.0 + K1 * O2 )
    """
    # REAL A0, B0, C0, A1, B1, C1
    # REAL(kind=dp) K1, K2
    K1 = GEOS_STD(A0, B0, C0)
    K2 = GEOS_STD(A1, B1, C1)
    return K1 / (1.0 + K1 * O2)


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

    # REAL(kind=dp) O3J, T3I
    T3I = 1.0*INVTEMP
    fh2o = (
        1.45e-10 * exp(89.0 * T3I) * H2O
        / (
            1.45e-10 * exp(89.0 * T3I) * H2O
            + 2.14e-11 * exp(110.0 * T3I) * N2
            + 3.20e-11 * exp(70.0 * T3I) * O2
        )
    )
    return O3J * fh2o


def GEOS_JO3_2(O3J):
    T3I = 1.0*INVTEMP
    RO1DplH2O = 1.63e-10 * exp(60.0 * T3I) * H2O
    RO1DplH2 = 1.2e-10 * H2
    RO1DplN2 = 2.15e-11 * exp(110.0 * T3I) * N2
    RO1DplO2 = 3.30e-11 * exp(55.0 * T3I) * O2
    RO1D = RO1DplH2O + RO1DplH2 + RO1DplN2 + RO1DplO2
    return O3J * RO1DplH2 / RO1D
