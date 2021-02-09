__all__ = ['CMAQ_1to4', 'CMAQ_5', 'CMAQ_6', 'CMAQ_7',
           'CMAQ_8', 'CMAQ_9', 'CMAQ_10', 'CMAQ_10D', 'OH_CO']
from numpy import exp, log10

# Globals must be udpated by update_func_world for this
# to work
P = M = TEMP = 0


def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    globals().update(world)


def CMAQ_1to4(A0, B0, C0):
    """
    CMAQ reaction rates form 1-4 have the form K = A * (T/300.0)**B * EXP(-C/T)
    """
    # REAL A0, B0, C0

    return (A0 * (TEMP/300.0)**B0 * exp(-C0/TEMP))


def CMAQ_5(A0, B0, C0, Kf):
    """
    CMAQ reaction form 5

    K1 = CMAQ_1to4(A0, B0, C0)

    Returns Kf / K1
    """
    # REAL A0, B0, C0
    # REAL(kind=dp) K1, Kf
    K1 = CMAQ_1to4(A0, B0, C0)
    return Kf / K1


def CMAQ_6(A0, B0, C0, Kf):
    """
    CMAQ reaction form 6

    K1 = CMAQ_1to4(A0, B0, C0)

    Returns Kf * K1
    """
    # REAL A0, B0, C0
    # REAL(kind=dp) K1, Kf
    K1 = CMAQ_1to4(A0, B0, C0)
    return Kf * K1


def CMAQ_7(A0, B0, C0):
    """
    CMAQ reaction form 6

    K0 = CMAQ_1to4(A0, B0, C0)

    Returns K0 * (1 + .6 * PRESS / 101325.) # Pressure is in Pascals
    """
    # REAL A0, B0, C0
    # REAL(kind=dp) K0
    K0 = CMAQ_1to4(A0, B0, C0)
    return K0 * (1 + .6 * P / 101325.)  # Pressure is in Pascals


def CMAQ_8(A0, C0, A2, C2, A3, C3):
    """
    CMAQ reaction form 8

    K0 = (A0) * exp(-(C0) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)
    K3 = (A3) * exp(-(C3) / TEMP)
    K3 = K3 * M

    Returns K0 + K3 / (1.0 + K3 / K2 )
    """
    # REAL A0, C0, A2, C2, A3, C3
    # REAL(kind=dp) K0, K2, K3
    K0 = (A0) * exp(-(C0) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)
    K3 = (A3) * exp(-(C3) / TEMP)
    K3 = K3 * M
    return K0 + K3 / (1.0 + K3 / K2)


def CMAQ_9(A1, C1, A2, C2):
    """
    CMAQ reaction rate form 9

    K1 = (A1) * exp(-(C1) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)

    M = third body concentration (molecules/cm3) and must be
        defined in the stdfuncs namespace

    Returns K1 + K2 * M
    """
    # REAL(kind=dp) A1, C1, A2, C2
    # REAL(kind=dp) K1, K2
    K1 = (A1) * exp(-(C1) / TEMP)
    K2 = (A2) * exp(-(C2) / TEMP)
    return K1 + K2 * M


def CMAQ_10(A0, B0, C0, A1, B1, C1, CF, N):
    """
    CMAQ reaction rate form 10

    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0 * M
    K1 = K0 / K1

    M = third body concentration (molecules/cm3) and must be
        defined in the stdfuncs namespace

    Returns (K0 / (1.0 + K1))*   \
         (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))
    """
    # REAL A0, B0, C0, A1, B1, C1, CF, N
    # REAL(kind=dp) K0, K1
    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0 * M
    K1 = K0 / K1
    return (K0 / (1.0 + K1)) *   \
        (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))


def CMAQ_10D(A0, B0, C0, A1, B1, C1, CF, N):
    """
    Same as reaction rate form 10, but implemented
    to provide compatibility for fortran code
    that need a DOUBLE form
    """
    return CMAQ_10(A0, B0, C0, A1, B1, C1, CF, N)


def OH_CO(A0, B0, C0, A1, B1, C1, CF, N):
    """
    OH + CO reaction rate

    *Note: Mostly like CMAQ_10, but slight difference in K1

    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0
    K1 = K0 / (K1 / M)

    M = third body concentration (molecules/cm3) and must be
        defined in the stdfuncs namespace

    return (K0 / (1.0 + K1))*   \
         (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))

    """
    # REAL A0, B0, C0, A1, B1, C1, CF, N
    # REAL(kind=dp) K0, K1
    K0 = CMAQ_1to4(A0, B0, C0)
    K1 = CMAQ_1to4(A1, B1, C1)
    K0 = K0
    K1 = K0 / (K1 / M)
    return (K0 / (1.0 + K1)) *   \
        (CF)**(1.0 / (1.0 / (N) + (log10(K1))**2))
