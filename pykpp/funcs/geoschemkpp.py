__all__ = [
    'HET', 'PHOTOL', 'GCARR', 'GCIUPAC3', 'GCJPL3', 'GCJPLEQ', 'GCJPLPR',
    'GC_DMSOH', 'GC_GLYCOHA', 'GC_GLYCOHB', 'GC_GLYXNO3', 'GC_HACOHA',
    'GC_HACOHB', 'GC_HO2NO3', 'GC_OHCO', 'GC_OHHNO3', 'GC_RO2HO2',
    'GC_RO2NO', 'GC_TBRANCH', 'HO2_H2O', 'OH_O1D', 'update_func_world'
]
from warnings import warn
from numpy import exp, log10, float64

DBLE = float64
N2 = O2 = M = TEMP = INVTEMP = 0
NUMDEN = int_HO2 = C = H2O = H2 = 0


def HET(spc_idx, rct_idx):
    return 0.


def PHOTOL(pidx):
    from ..tuv.tuv5pt0 import TUV_J5pt0 as TUV_J
    gc2tuv = {
        2: "TUV_J(2, THETA)",  # O3 + hv = O + O2
        3: "TUV_J(3, THETA)",  # O3 + hv = O1D + O2
        1: "TUV_J(1, THETA)",  # O2 + hv = 2.000O
        11: "TUV_J(6, THETA)",  # NO2 + hv = NO + O
        9: "TUV_J(5, THETA)",  # H2O2 + hv = OH + OH
        10: "TUV_J(31, THETA)",  # MP + hv = CH2O + HO2 + OH
        7: "TUV_J(17, THETA)",  # CH2O + hv = HO2 + H + CO
        8: "TUV_J(18, THETA)",  # CH2O + hv = H2 + CO
        16: "TUV_J(13, THETA)",  # HNO3 + hv = OH + NO2
        15: "TUV_J(12, THETA)",  # HNO2 + hv = OH + NO
        17: "TUV_J(14, THETA)*.5",  # HNO4 + hv = OH + NO3
        18: "TUV_J(14, THETA)*.5",  # HNO4 + hv = HO2 + NO2
        12: "TUV_J(8, THETA)",  # NO3 + hv = NO2 + O
        13: "TUV_J(7, THETA)",  # NO3 + hv = NO + O2
        14: "TUV_J(10, THETA)",  # N2O5 + hv = NO3 + NO2
        61: "TUV_J(19, THETA)+TUV_J(21, THETA)",  # ALD2 + hv = 0.880MO2 + HO2
                                                  #      + 0.880CO + 0.120MCO3
        62: "TUV_J(20, THETA)",  # ALD2 + hv = CH4 + CO
        59: "TUV_J(40, THETA)+TUV_J(41, THETA)",  # PAN + hv = 0.700MCO3
                                                  #            + 0.700NO2
                                                  #            + 0.300MO2
                                                  #            + 0.300NO3
        70: "TUV_J(22, THETA)",  # RCHO + hv = ETO2 + HO2 + CO
        76: "TUV_J(26, THETA)*.5",  # half?ACET + hv = MCO3 + MO2
        77: "TUV_J(26, THETA)*.5",  # half?ACET + hv = 2.000MO2 + CO
        69: "TUV_J(28, THETA)",  # MEK + hv = 0.850MCO3 + 0.850ETO2 + 0.150MO2
                                 #            + 0.150RCO3
        68: "TUV_J(23, THETA)",  # GLYC + hv = 0.900CH2O + 1.730HO2 + CO
                                 #             + 0.070OH + 0.100MOH
        72: "TUV_J(44, THETA)",  # GLYX + hv = 2.000HO2 + 2.000CO
        73: "TUV_J(45, THETA)*0.",  # none?GLYX + hv = H2 + 2.000CO
        74: "TUV_J(45, THETA)",  # GLYX + hv = CH2O + CO
        71: "TUV_J(46, THETA)",  # MGLY + hv = MCO3 + CO + HO2
        63: "TUV_J(27, THETA)/3",  # third?MVK + hv = PRPE + CO
        64: "TUV_J(27, THETA)/3",  # third?MVK + hv = MCO3 + CH2O + CO + HO2
        65: "TUV_J(27, THETA)/3",  # third?MVK + hv = MO2 + RCO3
        66: "TUV_J(25, THETA)",  # MACR + hv = CO + HO2 + CH2O + MCO3
        75: "TUV_J(29, THETA)+TUV_J(30, THETA)",  # HAC + hv = MCO3+CH2O + HO2
        78: "TUV_J(31, THETA)",  # INPN + hv = OH + HO2 + RCHO + NO2
        79: "TUV_J(31, THETA)",  # PRPN + hv = OH + HO2 + RCHO + NO2
        80: "TUV_J(31, THETA)",  # ETP + hv = OH + HO2 + ALD2
        81: "TUV_J(31, THETA)",  # RA3P + hv = OH + HO2 + RCHO
        82: "TUV_J(31, THETA)",  # RB3P + hv = OH + HO2 + ACET
        83: "TUV_J(31, THETA)",  # R4P + hv = OH + HO2 + RCHO
        84: "TUV_J(31, THETA)",  # PP + hv = OH + HO2 + ALD2 + CH2O
        85: "TUV_J(31, THETA)",  # RP + hv = OH + HO2 + ALD2
        86: "TUV_J(31, THETA)",  # RIP + hv = OH + HO2 + 0.710CH2O + 0.425MVK
                                 #            +  0.285MACR + 0.290HC5
        87: "TUV_J(31, THETA)",  # IAP + hv = OH + HO2 + 0.670CO + 0.190H2
                                 #         +  0.360HAC + 0.260GLYC + 0.580MGLY
        88: "TUV_J(31, THETA)",  # ISNP + hv = OH + HO2 + RCHO + NO2
        89: "TUV_J(31, THETA)",  # VRP + hv = OH + 0.300HO2 + 0.300CH2O
                                 #            + 0.700MCO3 + 0.700GLYC
                                 #            + 0.300MGLY
        90: "TUV_J(31, THETA)",  # MRP + hv = OH+HO2+HAC + 0.500CO + 0.500CH2O
        91: "TUV_J(31, THETA)",  # MAOP + hv = OH + CH2O + MCO3
        98: "TUV_J(36, THETA)",  # chose iso-butane nitrate?
                                 # R4N2 + hv = NO2 + 0.320ACET + 0.190MEK
                                 #            + 0.180MO2 + 0.270HO2 + 0.320ALD2
                                 #            + 0.130RCHO + 0.050A3O2
                                 #            + 0.180B3O2 + 0.320ETO2
        99: "TUV_J(33, THETA)",  # MAP + hv = OH + MO2
        92: "TUV_J(37, THETA)",  # probably not right index?
                                 # MACRN + hv = NO2 + HAC + MGLY + 0.500CH2O
                                 #              + HO2 + 0.500CO
        93: "TUV_J(37, THETA)",  # probably not right index?
                                 # MVKN + hv = GLYC + NO2 + MCO3
        94: "TUV_J(37, THETA)",  # probably not right index
                                 # ?ISOPNB + hv = HC5 + NO2 + HO2
        23: "TUV_J(54, THETA)",  # Br2 + hv = 2.000Br
        28: "TUV_J(55, THETA)",  # BrO + hv = Br + O
        32: "TUV_J(56, THETA)",  # HOBr + hv = Br + OH
        29: "TUV_J(58, THETA)",  # BrNO3 + hv = Br + NO3
        30: "TUV_J(57, THETA)",  # BrNO3 + hv = BrO + NO2
        31: "0.#TUV_J(, THETA)",  # needs update?BrNO2 + hv = Br + NO2
        56: "0.#TUV_J(, THETA)",  # needs update?CHBr3 + hv = 3.000Br
        104: "0.#TUV_J(, THETA)",  # needs update?MPN + hv = CH2O + NO3 + HO2
        105: "0.#TUV_J(, THETA)",  # needs update?MPN + hv = MO2 + NO2
        95: "TUV_J(36, THETA)",  # ISOPND + hv = HC5 + NO2 + HO2
        96: "TUV_J(36, THETA)",  # PROPNN + hv = CH2O + NO2 + CO + MO2
        97: "TUV_J(31, THETA)",  # ATOOH + hv = OH + CH2O + MCO3
        36: "TUV_J(11, THETA)",  # N2O + hv = N2 + O1D
        34: "0.#TUV_J(, THETA)",  # OCS + hv = SO2 + CO
        100: "0.#TUV_J(, THETA)",  # SO4 + hv = SO2 + 2.000OH
        6: "0.#TUV_J(, THETA)",  # NO + hv = O + N
        50: "TUV_J(77, THETA)",  # CH3Br + hv = MO2 + Br
        33: "0.#TUV_J(, THETA)",  # not found?BrCl + hv = Br + Cl
        22: "TUV_J(49, THETA)",  # Cl2 + hv = 2.000Cl
        27: "TUV_J(85, THETA)",  # ClO + hv = Cl + O
        25: "0.#TUV_J(, THETA)",  # not found?OClO + hv = ClO + O
        26: "TUV_J(50, THETA)",  # Cl2O2 + hv = Cl + ClOO
        21: "TUV_J(86, THETA)",  # ClNO2 + hv = Cl + NO2
        19: "TUV_J(52, THETA)",  # ClNO3 + hv = Cl + NO3
        20: "TUV_J(53, THETA)",  # ClNO3 + hv = ClO + NO2
        24: "0.#TUV_J(, THETA)",  # not found? HOCl + hv = Cl + OH
        43: "TUV_J(59, THETA)",  # CH3Cl + hv = MO2 + Cl
        44: "TUV_J(69, THETA)",  # CH3CCl3 + hv = 3.000Cl
        42: "TUV_J(61, THETA)",  # CCl4 + hv = 4.000Cl
        37: "TUV_J(67, THETA)",  # CFC11 + hv = 3.000Cl
        38: "TUV_J(68, THETA)",  # CFC12 + hv = 2.000Cl
        39: "TUV_J(64, THETA)",  # CFC113 + hv = 3.000Cl
        40: "TUV_J(65, THETA)",  # CFC114 + hv = 2.000Cl
        41: "TUV_J(66, THETA)",  # CFC115 + hv = Cl
        47: "TUV_J(70, THETA)",  # HCFC123 + hv = 2.000Cl
        48: "TUV_J(72, THETA)",  # HCFC141b + hv = 2.000Cl
        49: "TUV_J(73, THETA)",  # HCFC142b + hv = 2.000Cl
        46: "TUV_J(76, THETA)",  # HCFC22 + hv = 2.000Cl
        53: "TUV_J(79, THETA)",  # H1301 + hv = Br
        51: "TUV_J(82, THETA)",  # H1211 + hv = Cl + Br
        54: "TUV_J(80, THETA)",  # H2402 + hv = 2.000Br
        55: "0.#TUV_J(, THETA)",  # not found?CH2Br2 + hv = 2.000Br
        102: "TUV_J(51, THETA)",  # ClOO + hv = Cl + O2
        }
    return eval(gc2tuv[pidx])


def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    globals().update(world)
    globals()['NUMDEN'] = world['M']
    globals()['INVTEMP'] = 1 / world['TEMP']
    globals()['T300I'] = 300. / world['TEMP']


def OH_O1D(J, H2O, TEMP, NUMDEN):
    """
    REAL*8 J, H2O, TEMP, NUMDEN
    REAL*8 K1, K2, K3
    REAL*8 N2, O2
    """
    N2 = 0.79e0
    O2 = 0.21e0

    K1 = 1.63e-10 * exp(60e0 / TEMP)
    K2 = 2.15e-11 * exp(110e0 / TEMP)
    K3 = 3.30e-11 * exp(55e0 / TEMP)

    return J * K1 * H2O / (K1 * H2O + K2 * N2 * NUMDEN + K3 * O2 * NUMDEN)


def HO2_H2O(H2O, TEMP):
    # REAL*8 TEMP, H2O
    return 1 + 1.4e-21 * H2O * exp(2200 / TEMP)


def GCARR(A0, B0, C0):
    # REAL A0,B0,C0
    return float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)


def GC_HO2NO3(A0, B0, C0, A1, B1, C1):
    # REAL A0,B0,C0,A1,B1,C1
    # REAL(kind=dp) :: R0,R1
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    R1 = float(A1) * exp(float(C1) / TEMP) * (300. / TEMP)**float(B1)
    return (R0 + R1 * NUMDEN)*(1.e0 + 1.4E-21 * H2O * exp(2200.E+0 / TEMP))


def GC_TBRANCH(A0, B0, C0, A1, B1, C1):
    # Temperature Dependent Branching Ratio
    # REAL A0,B0,C0,A1,B1,C1
    # REAL(kind=dp) :: R0,R1
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    R1 = float(A1) * exp(float(C1) / TEMP) * (300. / TEMP)**float(B1)

    return R0/(1.e0+R1)


def GC_RO2HO2(A0, B0, C0, A1, B1, C1):
    # Carbon Dependence of RO2+HO2
    # REAL A0,B0,C0,A1,B1,C1
    # REAL(kind=dp) :: R0,R1
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    R1 = float(A1) * exp(float(C1) / TEMP) * (300. / TEMP)**float(B1)

    return R0 * (1E0 - exp(-0.245E0 * R1))


def GC_DMSOH(A0, B0, C0, A1, B1, C1):
    # DMS+OH+O2
    # REAL A0,B0,C0,A1,B1,C1
    # REAL(kind=dp) :: R0,R1
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    R1 = float(A1) * exp(float(C1) / TEMP) * (300. / TEMP)**float(B1)
    #    GC_DMSOH = R0/(1e0+R1*0.2095e0)
    return (R0 * NUMDEN * 0.2095e0)/(1e0 + R1 * 0.2095e0)


def GC_GLYXNO3(A0, B0, C0):
    # ---  K = K1*([O2]+3.5D18)/(2*[O2]+3.5D18)
    # --- HO2+2*CO branch of GLYX+OH/NO3
    # REAL A0,B0,C0
    # REAL(kind=dp) R0
    # REAL(kind=dp) O2

    O2 = NUMDEN * 0.2095e0
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    return R0 * (O2 + 3.5E+18)/(2.E+0 * O2 + 3.5E+18)


def GC_OHHNO3(A0, B0, C0, A1, B1, C1, A2, B2, C2):
    # ---  OH + HNO3:   K = K0 + K3[M] / (1 + K3[M]/K2)  ------
    # REAL A0,B0,C0,A1,B1,C1,A2,B2,C2
    # REAL(kind=dp) R0,R1,R2
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    R1 = float(A1) * exp(float(C1) / TEMP) * (300. / TEMP)**float(B1)
    R2 = NUMDEN * (
        float(A2) * exp(float(C2) / TEMP) * (300. / TEMP)**float(B2)
    )
    return R0 + R2 / (1.E0 + R2 / R1)


def GC_GLYCOHA(A0, B0, C0):
    # REAL A0,B0,C0,R0,GLYC_FRAC
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    GLYC_FRAC = 1e+0 - 11.0729e+0 * exp(-(1. / 73.) * TEMP)
    if (GLYC_FRAC < 0e+0):
        GLYC_FRAC = 0e+0
    return R0 * GLYC_FRAC


def GC_GLYCOHB(A0, B0, C0):
    # REAL A0,B0,C0
    # REAL(kind=dp) :: R0,GLYC_FRAC
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    GLYC_FRAC = 1e+0 - 11.0729e+0 * exp(-(1. / 73.) * TEMP)
    if (GLYC_FRAC < 0e+0):
        GLYC_FRAC = 0e+0
    return R0 * (1e0 - GLYC_FRAC)

    # END FUNCTION GC_GLYCOHB


def GC_HACOHA(A0, B0, C0):
    # REAL A0,B0,C0
    # REAL(kind=dp) :: R0,HAC_FRAC
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    HAC_FRAC = 1e+0 - 23.7e+0 * exp(-(1. / 60.) * TEMP)
    if (HAC_FRAC < 0e+0):
        HAC_FRAC = 0e+0
    return R0 * HAC_FRAC


def GC_HACOHB(A0, B0, C0):
    # REAL A0,B0,C0
    # REAL(kind=dp) :: R0,HAC_FRAC
    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    HAC_FRAC = 1e+0 - 23.7e+0 * exp(-(1. / 60.) * TEMP)
    if (HAC_FRAC < 0e+0):
        HAC_FRAC = 0e+0
    return R0 * (1.E0 - HAC_FRAC)


def GC_OHCO(A0, B0, C0):
    # REAL A0,B0,C0,R0
    # REAL KLO1,KLO2,KHI1,KHI2,XYRAT1,XYRAT2,BLOG1,BLOG2,FEXP1,FEXP2
    # REAL KCO1,KCO2,KCO
    # R0 = float(A0) * exp(float(C0)/TEMP) * (300./TEMP)**float(B0)
    # R0 = R0 * (1.E+0 + 0.6e+0*9.871E7*PRESS)
    # new OH+CO rate from JPL2006.
    KLO1 = 5.9E-33 * (300. / TEMP)**(1.4E+0)
    KHI1 = 1.1E-12 * (300. / TEMP)**(-1.3E0)
    XYRAT1 = KLO1 * NUMDEN / KHI1
    BLOG1 = log10(XYRAT1)
    FEXP1 = 1.E+0 / (1.E+0 + BLOG1 * BLOG1)
    KCO1 = KLO1 * NUMDEN * 0.6**FEXP1 / (1.e+0 + XYRAT1)
    KLO2 = 1.5E-13 * (300. / TEMP)**(-0.6E+0)
    KHI2 = 2.1e+09 * (300. / TEMP)**(-6.1E+0)
    XYRAT2 = KLO2*NUMDEN/KHI2
    BLOG2 = log10(XYRAT2)
    FEXP2 = 1.E+0 / (1.E+0 + BLOG2 * BLOG2)
    KCO2 = KLO2 * 0.6**FEXP2 / (1.e+0 + XYRAT2)
    KCO = KCO1 + KCO2
    return KCO


def GC_RO2NO(B, A0, B0, C0, A1, B1, C1):
    # ---  K = K1*(1-FYRNO3(K2,M,T))  ---  abstraction branch of RO2+NO
    # CHARACTER(*) B !Branch Toggle
    # REAL A0,B0,C0,A1,B1,C1
    # REAL(kind=dp) :: R0,R1
    # REAL(kind=dp) :: YYYN, XXYN,  AAA,  RARB, ZZYN
    # REAL(kind=dp) :: XF, ALPHA, Y300, BETA, XMINF, XM0
    # REAL(kind=dp) :: FYRNO3
    Y300 = 0.826
    ALPHA = 1.94e-22
    BETA = 0.97
    XM0 = 0.
    XMINF = 8.1
    XF = 0.411

    R0 = float(A0) * exp(float(C0) / TEMP) * (300. / TEMP)**float(B0)
    R1 = float(A1) * exp(float(C1) / TEMP) * (300. / TEMP)**float(B1)
    # Initialize static variables

    XXYN = ALPHA * exp(BETA * R1) * NUMDEN * ((300. / TEMP)**XM0)
    YYYN = Y300 * ((300. / TEMP)**XMINF)
    AAA = log10(XXYN / YYYN)
    ZZYN = 1. / (1. + AAA * AAA)
    RARB = (XXYN / (1. + (XXYN / YYYN))) * (XF**ZZYN)
    FYRNO3 = RARB / (1. + RARB)
    if (B.strip() == 'A'):
        return R0 * FYRNO3
    elif (B.strip() == 'B'):
        return R0 * (1.E+0 - FYRNO3)


def GCJPL3(k0_300, n, ki_300, m):
    #  Functions given in JPL Booklet
    # REAL k0_300, n, ki_300,m
    # REAL k0, ki
    k0 = k0_300 * ((TEMP / 300.e0)**(-n))
    ki = ki_300 * ((TEMP / 300.e0)**(-m))
    #      GCJPL3=(k0*NUMDEN)/(1+k0*NUMDEN/ki)*0.6** &
    #    ((1+((log10(k0*NUMDEN/ki))**2d0)**-1.0d0))
    GCJPL3 = (
        k0 / (1.e0 + k0 / (ki / NUMDEN))
    ) * 0.6**((1 + ((log10(k0 / (ki / NUMDEN)))**2e0)**1.0e0))
    return GCJPL3 * NUMDEN


def GCJPLEQ(A0, B0, C0, A1, B1, C1, A2, B2, C2, FV, FCT1, FCT2):
    # Function calculates the rate constant of the forward reaction
    # calculates the equilibrium constant
    # Find the backwards reaction by K=kforward/kbackwards
    # REAL A0,B0,C0,A1,B1,C1
    # REAL(kind=dp) :: R0,R1
    # REAL, OPTIONAL :: A2,B2,C2,FV,FCT1,FCT2 !If a P-dependent rxn
    # Calculate Backwards reaction
    R0 = GCARR(A0, B0, C0)
    # Calculate forwards reaction
    if (A2 is not None):  # P-dependent
        if (
            B2 is None or C2 is None or FV is None or FCT1 is None
            or FCT2 is None
        ):
            print('GCJPLEQ: Missing parameters for P-dependent reaction.')
            # Missing params!
            print('GCJPLEQ: Returning zero')
            return 0.E0
        else:
            R1 = GCJPLPR(A1, B1, C1, A2, B2, C2, FV, FCT1, FCT2)
    else:
        R1 = GCARR(A1, B1, C1)  # Std. Arrhenius eqn.

    return R1 / R0


def GCJPLPR(A0, B0, C0, A1, B1, C1, FV, FCT1, FCT2):
    # * PRESSURE-DEPENDENT EFFECTS
    # * ADD THE THIRD BODY EFFECT FOR PRESSURE DEPENDENCE OF RATE
    # * COEFFICIENTS.
    # A0 B0, & C0 are the Arrhenius parameters for the lower-limit
    # rate. A1, B1 & C1 are the upper-limit parameters.
    # FV is the falloff curve paramter, (SEE ATKINSON ET. AL (1992)
    # J. PHYS. CHEM. REF. DATA 21, P. 1145). USUALLY = 0.6
    #
    # REAL A0,B0,C0,A1,B1,C1,FV,FCT1,FCT2
    # REAL FCT,XYRAT,BLOG,RLOW,RHIGH,FEXP

    RLOW = GCARR(A0, B0, C0) * NUMDEN
    RHIGH = GCARR(A1, B1, C1)

    if (FCT2 != 0.):
        FCT = exp(-TEMP / FCT1) + exp(-FCT2 / TEMP)
        XYRAT = RLOW/RHIGH
        BLOG = log10(XYRAT)
        FEXP = 1.e+0 / (1.e+0 + BLOG * BLOG)
        return RLOW * FCT**FEXP / (1e+0 + XYRAT)
    elif (FCT1 != 0.):
        FCT = exp(-TEMP / FCT1)
        XYRAT = RLOW / RHIGH
        BLOG = log10(XYRAT)
        FEXP = 1.e+0 / (1.e+0 + BLOG * BLOG)
        return RLOW * FCT**FEXP / (1e+0 + XYRAT)
    else:
        XYRAT = RLOW / RHIGH
        BLOG = log10(XYRAT)
        FEXP = 1.e+0 / (1.e+0 + BLOG * BLOG)
        return RLOW * FV**FEXP / (1e+0 + XYRAT)


def GCIUPAC3(ko_300, n, ki_300, m, Fc):
    # Function calcualtes the rate constant of 3 body reaction using IUPAC
    # methology
    # REAL ko_300,n,ki_300,m,Fc
    # REAL ko, ki, F, NN

    ko = ko_300 * ((TEMP / 300.e0)**n) * NUMDEN
    ki = ki_300 * ((TEMP / 300.e0)**m)

    NN = 0.75 - 1.27 * log10(Fc)
    F = 10.0**(log10(Fc) / (1.0e0 + (log10(ko / ki) / NN)**2.0))

    return ko / (1 + ko / ki) * F
