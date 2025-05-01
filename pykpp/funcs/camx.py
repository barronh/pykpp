__all__ = ['CAMX_4', 'CAMX_6']

from numpy import exp, log

TEMP = M = 0


def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    global M, TEMP
    globals().update(world)


def CAMX_4(A0, Ea0, B0, Tr0, A1, Ea1, B1, Tr1, F, n):
    k0 = A0 * (TEMP/Tr0)**(B0) * exp(-Ea0/TEMP)
    kinf = A1 * (TEMP/Tr1)**(B1) * exp(-Ea1/TEMP)
    k0M = k0*M
    G = (1+(log(k0M/kinf)/n)**2)**(-1)
    k = (k0M / (1 + k0M/kinf)) * F**(G)
    return k


def CAMX_6(A0, Ea0, B0, Tr0, A2, Ea2, B2, Tr2, A3, Ea3, B3, Tr3):
    k0 = A0 * (TEMP/Tr0)**(B0) * exp(-Ea0/TEMP)
    k2 = A2 * (TEMP/Tr2)**(B2) * exp(-Ea2/TEMP)
    k3 = A3 * (TEMP/Tr3)**(B3) * exp(-Ea3/TEMP)
    k = k0 + k3*M/(1+k3*M/k2)
    return k
