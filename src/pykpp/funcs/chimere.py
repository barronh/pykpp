__all__ = ['CHIMERE_MTROE', 'CHIMERE_TROE', 'CHIMERE_JO3', 'CHIMERE_SPECIAL_1', 'CHIMERE_SPECIAL_2', 'CHIMERE_SPECIAL_3', 'CHIMERE_SPECIAL_4']

from numpy import *

def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    globals().update(world)

def CHIMERE_TROE(A0, B0, C0, A1, B1, C1, N):
    """
    Mapping:
        A0 = tabrate(1,nr)
        B0 = tabrate(2,nr)
        C0 = tabrate(3,nr)
        A1 = tabrate(4,nr)
        B1 = tabrate(5,nr)
        C1 = tabrate(6,nr)
        N  = tabrate(7,nr)

        M  = ai; M = third body concentration (molecules/cm3) and must be 
        defined in the stdfuncs namespace
    
        TEMP = te = bulk air temperature
        
        1. = dun
    Original Code:
        c1 = tabrate(1,nr)*exp(-tabrate(2,nr)/te)                 &
          *(300d0/te)**tabrate(3,nr)
        c2 = tabrate(4,nr)*exp(-tabrate(5,nr)/te)                 &
          *(300d0/te)**tabrate(6,nr)
        c3 = ai*c1
        c4 = c3/c2
        ex = dun/(dun + log10(c4)**2)
        rate(nr,izo,ime,ivert) = c1*tabrate(7,nr)**ex/(dun + c4)
    """

    c1 = A0*exp(-B0/TEMP)*(300./TEMP)**C0
    c2 = A1*exp(-B1/TEMP)*(300./TEMP)**C1
    c3 = M*c1
    c4 = c3/c2
    ex = 1./(1. + log10(c4)**2)
    out = c1*N**ex/(1. + c4)
    return out

def CHIMERE_MTROE(A0, B0, C0, A1, B1, C1, N):
    """
    Mapping:
        A0 = tabrate(1,nr)
        B0 = tabrate(2,nr)
        C0 = tabrate(3,nr)
        A1 = tabrate(4,nr)
        B1 = tabrate(5,nr)
        C1 = tabrate(6,nr)
        N  = tabrate(7,nr)
        M  = ai
        TEMP = te
        1. = dun
    Original Code:
        c1 = tabrate(1,nr)*exp(-tabrate(2,nr)/te)                 &
          *(300d0/te)**tabrate(3,nr)
        c2 = tabrate(4,nr)*exp(-tabrate(5,nr)/te)                 &
          *(300d0/te)**tabrate(6,nr)
        c3 = ai*c1
        c4 = c3/c2
        ex = dun/(dun + ((log10(c4) - 0.12d0)/1.2d0)**2)
        rate(nr,izo,ime,ivert) = c1*tabrate(7,nr)**ex/(dun + c4)
    """
    c1 = A0*exp(-B0/TEMP)*(300./TEMP)**C0
    c2 = A1*exp(-B1/TEMP)*(300./TEMP)**C1
    c3 = M*c1
    c4 = c3/c2
    ex = 1./(1. + ((log10(c4) - 0.12)/1.2)**2)
    out = c1*N**ex/(1. + c4)
    return out

def CHIMERE_JO3(rate):
    ai = M
    te = TEMP
    # Chimere does not use specific humidity, so this
    # is commented out in favor of H2O concentration
    # the units of sphu in the rates.F90 are #/cm**-3
    #
    #Ma = (O2 * 32. + N2 * 28.) / Avogadro
    #Mh2o = H2O * 18. / Avogadro
    #hu = Mh2o / Ma
    #factor = hu/(hu + ai*(0.02909*exp(70./te) + 0.06545*exp(110./te)))
    factor = H2O/(H2O + ai*(0.02909*exp(70./te) + 0.06545*exp(110./te)))
    return rate*factor

def CHIMERE_SPECIAL_1(A1, C1, A2, C2):
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    rate = f1 * f2/(1. + f2)
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    rate = f1 * f2/(1. + f2)
    return rate

def CHIMERE_SPECIAL_2(A1, C1, A2, C2):
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    rate = f1/(1. + f2)
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    rate = f1/(1. + f2)
    return rate


def CHIMERE_SPECIAL_3(A1, C1, A2, C2, A3, C3, A4, C4):
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    f3 = A3*exp(-C3/TEMP)
    f4 = A4*exp(-C4/TEMP)
    rate = 2.*(f1 * f2 * f3 * f4/((1.+f3)*(1.+f4)))**(0.5)
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    f3 = A3*exp(-C3/TEMP)
    f4 = A4*exp(-C4/TEMP)
    rate = 2.*(f1 * f2 * f3 * f4/((1.+f3)*(1.+f4)))**(0.5)
    return rate
    
def CHIMERE_SPECIAL_4(A1, C1, A2, C2, A3, C3, A4, C4):
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    f3 = A3*exp(-C3/TEMP)
    f4 = A4*exp(-C4/TEMP)
    f3 = f3 / (1. + f3)
    f4 = f4 / (1. + f4)
    rate = 2.0*(f1*f2)**(0.5)*(1.-(f3*f4)**(0.5))*(1.-f4)/(2.-f3-f4)
    """
    f1 = A1*exp(-C1/TEMP)
    f2 = A2*exp(-C2/TEMP)
    f3 = A3*exp(-C3/TEMP)
    f4 = A4*exp(-C4/TEMP)
    f3 = f3 / (1. + f3)
    f4 = f4 / (1. + f4)
    rate = 2.0*(f1*f2)**(0.5)*(1.-(f3*f4)**(0.5))*(1.-f4)/(2.-f3-f4)
    return rate
