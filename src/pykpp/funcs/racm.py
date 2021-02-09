__all__ = ['RACM_TROE', 'RACM_TROE_EQUIL', 'RACM_THERMAL', 'RACM_THERMAL_T2']
from numpy import exp, log10

TEMP = M = 0


def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    globals().update(world)


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
    # REAL A0, B0, A1, B1, C1
    # REAL(kind=dp) K0, K1
    # REAL(kind=dp), PARAMETER :: CF = 0.6_dp, N = 1._dp
    CF = 0.6
    N = 1.
    K0 = (A0 * (TEMP/300.0)**(-B0))
    K1 = (A1 * (TEMP/300.0)**(-B1))
    K0 = K0 * M
    K1 = K0 / K1
    return (K0 / (1.0 + K1)) *   \
        CF**(1.0 / (1.0 / N + (log10(K1))**2))


def RACM_TROE_EQUIL(A0, B0, A1, B1, A2, C2):
    """
    RACM Troe equilibrium equation as defined in the RACM SBOX model

    #REAL A0, B0, A1, B1, A2, C2
    return RACM_TROE( A0, B0, A1, B1) *  (1./A2 * exp(-C2 / TEMP))
    """
    # REAL A0, B0, A1, B1, A2, C2
    return RACM_TROE(A0, B0, A1, B1) * (1./A2 * exp(-C2 / TEMP))


def RACM_THERMAL(A0, B0):
    """
    RACM Thermal equation as defined in the RACM SBOX model

    #REAL A0, B0
    #   RACM2 reaction rates have the form K = A * EXP(-B / T)
    #
    #   Translation adds a 0 C
    return (A0 * exp(-B0 / TEMP))
    """
    # REAL A0, B0
    #   RACM2 reaction rates have the form K = A * EXP(-B / T)
    #
    #   Translation adds a 0 C
    return (A0 * exp(-B0 / TEMP))


def RACM_THERMAL_T2(A0, B0):
    """
    RACM Thermal T2 equation as defined in the RACM SBOX model

    # REAL A0, B0
    # REAL, PARAMETER :: C0 = 0.
    #
    #   Translation adds a 0 C
    return (A0)*TEMP**2*exp(-(B0)/TEMP)
    """
    # REAL A0, B0
    # REAL, PARAMETER :: C0 = 0.
    #
    #   Translation adds a 0 C
    return (A0)*TEMP**2*exp(-(B0)/TEMP)
