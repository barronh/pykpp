from numpy import *
from scipy.constants import *

def ARR( A0, B0, C0 ):
    out = A0 * exp(-B0 / ARR.TEMP) * (ARR.TEMP / 300.0)**(C0)
    return out

def ARR2( A0, B0):
    out = A0 * exp( B0 / ARR2.TEMP )
    return out

def EP2(A0, C0, A2, C2, A3, C3):
    K0 = A0 * exp(-C0 / EP2.TEMP)
    K2 = A2 * exp(-C2 / EP2.TEMP)
    K3 = A3 * exp(-C3 / EP2.TEMP)
    K3 = K3 * EP2.M * 1.0E6
    return K0 + K3 / (1.0 + K3 / K2 )

def EP3(A1, C1, A2, C2):
    K1 = A1 * exp(-C1 / EP3.TEMP)
    K2 = A2 * exp(-C2 / EP3.TEMP)
    return K1 + K2 * (1.0E6 * EP3.M)

def FALL ( A0, B0, C0, A1, B1, C1, CF):
    K0 = ARR(A0, B0, C0)
    K1 = ARR(A1, B1, C1)
    K0 = K0 * FALL.M * 1.0E6
    K1 = K0 / K1
    return (K0 / (1.0 + K1))* CF**(1.0 / (1.0 + (log10(K1))**2))

def k_3rd(temp, cair, k0_300K, n, kinf_300K, m, fc):
    zt_help = 300. / temp
    k0_T    = k0_300K   * zt_help**(n) * cair # k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        # k_inf at current T
    k_ratio = k0_T / kinf_T
    return k0_T / (1. + k_ratio) * fc**(1. / (1. + log10(k_ratio)**2))

tuvdata = """   1 O2 + hv -> O + O                          1.780E-28
   2 O3 -> O2 + O(1D)                          3.756E-05
   3 O3 -> O2 + O(3P)                          4.523E-04
   4 NO2 -> NO + O(3P)                         9.330E-03
   5 NO3 -> NO + O2                            2.288E-02
   6 NO3 -> NO2 + O(3P)                        1.798E-01
   7 N2O5 -> NO3 + NO + O(3P)                  2.155E-09
   8 N2O5 -> NO3 + NO2                         4.732E-05
   9 N2O + hv -> N2 + O(1D)                    9.188E-26
  10 HO2 + hv -> OH + O                        6.019E-23
  11 H2O2 -> 2 OH                              7.614E-06
  12 HNO2 -> OH + NO                           2.059E-03
  13 HNO3 -> OH + NO2                          7.139E-07
  14 HNO4 -> HO2 + NO2                         5.133E-06
  15 CH2O -> H + HCO                           3.293E-05
  16 CH2O -> H2 + CO                           4.801E-05
  17 CH3CHO -> CH3 + HCO                       5.724E-06
  18 CH3CHO -> CH4 + CO                        1.439E-11
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     2.012E-05
  21 CHOCHO -> products                        7.764E-05
  22 CH3COCHO -> products                      1.109E-04
  23 CH3COCH3                                  7.099E-07
  24 CH3OOH -> CH3O + OH                       5.773E-06
  25 CH3ONO2 -> CH3O+NO2                       1.018E-06
  26 PAN + hv -> products                      8.056E-07
  27 ClOO + hv -> Products                     1.904E-20
  28 ClONO2 + hv -> Cl + NO3                   4.388E-05
  29 ClONO2 + hv -> ClO + NO2                  8.948E-06
  30 CH3Cl + hv -> Products                    1.581E-26
  31 CCl2O + hv -> Products                    4.173E-24
  32 CCl4 + hv -> Products                     9.663E-24
  33 CClFO + hv -> Products                    2.367E-24
  34 CF2O + hv -> Products                     5.469E-26
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     2.220E-25
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     1.256E-26
  37 CF3CF2Cl (CFC-115) + hv -> Products       0.000E+00
  38 CCl3F (CFC-11) + hv -> Products           2.160E-24
  39 CCl2F2 (CFC-12) + hv -> Products          9.045E-26
  40 CH3CCl3 + hv -> Products                  3.680E-24
  41 CF3CHCl2 (HCFC-123) + hv -> Products      1.887E-25
  42 CF3CHFCl (HCFC-124) + hv -> Products      2.024E-27
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     2.651E-25
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     1.614E-27
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  4.454E-25
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  2.422E-27
  47 CHClF2 (HCFC-22) + hv -> Products         0.000E+00
  48 BrONO2 + hv -> Br + NO3                   4.037E-04
  49 BrONO2 + hv -> BrO + NO2                  9.885E-04
  50 CH3Br + hv -> Products                    1.682E-23
  51 CHBr3                                     1.444E-06
  52 CF3Br (Halon-1301) + hv -> Products       3.254E-24
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  3.333E-23
  54 CF2Br2 (Halon-1202) + hv -> Products      2.656E-09
  55 CF2BrCl (Halon-1211) + hv -> Products     4.059E-15
  56 Cl2 + hv -> Cl + Cl                       2.526E-03"""
  
tuvbyidx = {}
tuvbykey = {}
for line in tuvdata.split('\n'):
    val = float(line[47:])
    idx = int(line[:4])
    key = line[5:46].strip()
    tuvbyidx[idx] = val
    tuvbykey[key] = val

def TUV_J(idx, scaling = 1.):
    """
        Photolysis rates based on on-line TUV v4.1
        idx - either the integer index in TUV v4.1 or the reaction string (i.e., jlabel)
        scaling - a number to scale result by
    """
    if idx in tuvbyidx:
        return tuvbyidx[idx] * scaling
    elif idx in tuvbykey:
        return tuvbykey[idx] * scaling
    else:
        raise KeyError('Not in tuv data: \n%s' % tuvdata)
  
DP3 = EP3

def k_arr (k_298, tdep, temp):
    return k_298 * exp(tdep * (1. / temp - 3.3540E-3)) # 1/298.15=3.3540e-3

def update_func_world(world):
    try:
        EP3.M = EP2.M = FALL.M = world['M']
    except:
        try:
            EP3.M = EP2.M = FALL.M = eval('P / (R / centi**3) / TEMP * Avogadro', None, world)
        except:
            pass
    ARR2.TEMP = ARR.TEMP = EP2.TEMP = EP3.TEMP = world['TEMP']
