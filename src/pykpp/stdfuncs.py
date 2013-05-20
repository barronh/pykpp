from numpy import *

def ARR( A0, B0, C0, TEMP ):
    out = A0 * exp(-B0 / TEMP) * (TEMP / 300.0)**(C0)
    return out

def ARR2( A0, B0, TEMP):
    out = A0 * exp( B0 / TEMP )
    return out

def EP2(A0, C0, A2, C2, A3, C3, TEMP):
    K0 = A0 * exp(-C0 / TEMP)
    K2 = A2 * exp(-C2 / TEMP)
    K3 = A3 * exp(-C3 / TEMP)
    K3 = K3 * M * 1.0E6
    return K0 + K3 / (1.0 + K3 / K2 )

def EP3(A1, C1, A2, C2, TEMP):
    K1 = A1 * exp(-C1 / TEMP)
    K2 = A2 * exp(-C2 / TEMP)
    return K1 + K2 * (1.0E6 * M)

def FALL ( A0, B0, C0, A1, B1, C1, CF, TEMP, M):
    K0 = ARR(A0, B0, C0)
    K1 = ARR(A1, B1, C1)
    K0 = K0 * M * 1.0E6
    K1 = K0 / K1
    return (K0 / (1.0 + K1))* CF**(1.0 / (1.0 + (log10(K1))**2))

def k_3rd(temp, cair, k0_300K, n, kinf_300K, m, fc):
    zt_help = 300. / temp
    k0_T    = k0_300K   * zt_help**(n) * cair # k_0   at current T
    kinf_T  = kinf_300K * zt_help**(m)        # k_inf at current T
    k_ratio = k0_T / kinf_T
    return k0_T / (1. + k_ratio) * fc**(1. / (1. + log10(k_ratio)**2))


def k_arr (k_298, tdep, temp):
    return k_298 * exp(tdep * (1. / temp - 3.3540E-3)) # 1/298.15=3.3540e-3
