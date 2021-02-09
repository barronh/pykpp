__all__ = [
    'MCM_KBPAN', 'MCM_KFPAN', 'AM3_STD', 'AM3_TROE', 'AM3_USR1', 'AM3_USR2',
    'AM3_USR3', 'AM3_USR4', 'AM3_USR5', 'AM3_USR6', 'AM3_USR7', 'AM3_USR8',
    'AM3_USR8a', 'AM3_USR9', 'AM3_USR10', 'AM3_USR11', 'AM3_USR12',
    'AM3_USR13', 'AM3_USR14', 'AM3_USR15', 'AM3_USR16', 'AM3_USR17',
    'AM3_USR18', 'AM3_USR19', 'AM3_USR21', 'AM3_USR22', 'AM3_USR24',
    'AM3_USR25', 'AM3_USR51', 'AM3_USR52', 'AM3_USR53', 'AM3_USR54',
    'AM3_USR57', 'AM3_USR58', 'AM3_USR59', 'AM3_USR60', 'AM3_USR61'
]
from numpy import log10, exp

O2 = N2 = M = TEMP = INVTEMP = 0
H2O = T300I = 0


def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    global O2, N2, M, TEMP, INVTEMP
    globals().update(world)
    globals()['INVTEMP'] = 1 / world['TEMP']
    globals()['T300I'] = 300. / world['TEMP']


def MCM_KBPAN():
    KD0 = 1.10e-05 * M * exp(-10100 * INVTEMP)
    KDI = 1.90e17 * exp(-14100 * INVTEMP)
    KRD = KD0 / KDI
    FCD = 0.30
    NCD = 0.75-1.27*(log10(FCD))
    FD = 10**(log10(FCD)/(1+(log10(KRD)/NCD)**2))
    KBPAN = (KD0*KDI)*FD/(KD0+KDI)
    return KBPAN


def MCM_KFPAN():
    KC0 = 3.28e-28*M*(TEMP/300)**-6.87
    KCI = 1.125e-11*(TEMP/300)**-1.105
    KRC = KC0/KCI
    FCC = 0.30
    NC = 0.75-1.27*(log10(FCC))
    FC = 10**(log10(FCC)/(1+(log10(KRC)/NC)**2))
    KFPAN = (KC0*KCI)*FC/(KC0+KCI)
    return KFPAN


def AM3_STD(A0, B0, C0):
    return A0 * exp(B0*INVTEMP)


def AM3_TROE(A0, B0, A1, B1, factor):
    ko = A0 * (T300I)**B0
    kinf = A1 * (T300I)**B1
    xpo = ko * M / kinf
    AM3_TROE = ko * M / (1.0 + xpo)
    xpo = log10(xpo)
    xpo = 1.0 / (1.0 + xpo*xpo)
    return AM3_TROE * factor**xpo


def AM3_USR1():
    return 6.e-34 * (T300I)**2.4


def AM3_USR2():
    return AM3_TROE(2.0e-30, 4.4, 1.4e-12, .7, .6)


def AM3_USR3():
    return AM3_USR2() * 3.704e26 * exp(-11000.0*INVTEMP)


def AM3_USR4():
    return AM3_TROE(1.8e-30, 3.0, 2.8e-11, 0.0, .6)


def AM3_USR5():
    TINV = 1*INVTEMP
    ko = M * 6.5e-34 * exp(1335.*TINV)
    ko = ko / (1. + ko/(2.7e-17*exp(2199.*TINV)))
    return ko + 2.4e-14*exp(460.*TINV)


def AM3_USR6():
    return AM3_TROE(2.0e-31, 3.4, 2.9e-12, 1.1, .6)


def AM3_USR7():
    return AM3_USR6() * exp(-10900.*INVTEMP) / 2.1e-27


def AM3_USR8():
    return AM3_TROE(5.9e-33, 1.4, 1.1e-12, -1.3, .6)


def AM3_USR8a():
    return AM3_TROE(1.5e-13, -0.6, 2.1e9, -6.1, .6) / M


def AM3_USR9():
    tinv = 1.0*INVTEMP
    ko = 2.3e-13 * exp(600.0*tinv)
    kinf = 1.7e-33 * M * exp(1000.0*tinv)
    fc = 1.0 + 1.4e-21 * H2O * exp(2200.0*tinv)

    return (ko + kinf) * fc


def AM3_USR10():
    return AM3_TROE(8.e-27, 3.5, 3.e-11, 0, .5)


def AM3_USR11():
    return AM3_TROE(9.7e-29, 5.6, 9.3e-12, 1.5, .6)


def AM3_USR12():
    return AM3_USR11() * 1.111e28 * exp(-14000.0 * INVTEMP)


def AM3_USR13():
    return AM3_TROE(1.e-28, 4.5, 7.5e-12, 0.85, .6)


def AM3_USR14():
    return 9.3e-12 * T300I / M


def AM3_USR15():
    return AM3_USR14() * 1.111e28 * exp(-14000.0 * INVTEMP)


def AM3_USR16():
    return 0


def AM3_USR17():
    return 0


def AM3_USR18():
    return 0


def AM3_USR19():
    return 0


def AM3_USR21():
    return TEMP**2 * 7.69e-17 * exp(253.0*INVTEMP)


def AM3_USR22():
    return 3.82e-11 * exp(-2000.0*INVTEMP) + 1.33e-13


def AM3_USR24():
    ko = 1.0 + 5.5e-31 * exp(7460.0*INVTEMP) * M * 0.21
    return 1.7e-42 * exp(7810.0*INVTEMP) * M * 0.21 / ko


def AM3_USR25():
    return 0.0


def AM3_USR51():
    return AM3_TROE(9.0e-28, 8.9, 7.7e-12, .2, .6)


def AM3_USR52():
    return AM3_USR51()*1.111e28 * exp(-14000.0 * INVTEMP)


def AM3_USR53():
    return AM3_TROE(9.0e-28, 8.9, 7.7e-12, .2, .6)


def AM3_USR54():
    return AM3_USR53()*1.111e28 * exp(-14000.0 * INVTEMP)


def AM3_USR57():
    FRAC = 1-11.0729*exp(-(1./73.)*TEMP)
    if (FRAC < 0):
        FRAC = 0.01
    return AM3_STD(8.00e-12, 0.0, 0.0)*FRAC


def AM3_USR58():
    FRAC = 1-11.0729*exp(-(1./73.)*TEMP)
    if (FRAC < 0):
        FRAC = 0.01
    return AM3_STD(8.00e-12, 0.0, 0.0)*(1.0-FRAC)


def AM3_USR59():
    K1 = AM3_STD(1.40e-12, -1860.0, 0.0)
    return K1 * (O2 + 3.5e18) / (2.0 * O2 + 3.5e18)


def AM3_USR60():
    FRAC = 1-23.7*exp(-(1./60)*TEMP)
    if (FRAC < 0):
        FRAC = 0.01
    return AM3_STD(2.15e-12, 305, 0.0)*FRAC


def AM3_USR61():
    FRAC = 1-23.7*exp(-(1./60)*TEMP)
    if (FRAC < 0):
        FRAC = 0.01
    return AM3_STD(2.15e-12, 305, 0.0)*(1.0-FRAC)
