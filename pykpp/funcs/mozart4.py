__all__ = [
    'MZ4_TROE', 'MZ4_USR1', 'MZ4_USR2', 'MZ4_USR3', 'MZ4_USR4', 'MZ4_USR5',
    'MZ4_USR6', 'MZ4_USR7', 'MZ4_USR8', 'MZ4_USR9', 'MZ4_USR10', 'MZ4_USR11',
    'MZ4_USR12', 'MZ4_USR14', 'MZ4_USR21', 'MZ4_USR22', 'MZ4_USR23',
    'MZ4_USR24'
]
from numpy import exp, log10

H2O = M = TEMP = 0


def update_func_world(mech, world):
    """
    Function to update globals for user defined functions
    """
    globals().update(world)


def MZ4_TROE(A0, B0, A1, B1, factor):
    """
    Troe fall off equation as calculated in
    MOZART4
    """
    # REAL(kind=dp) A0, B0, factor, A1, B1
    # REAL(kind=dp) ko, kinf, xpo
    ko = A0 * (300.e0 / TEMP)**B0
    kinf = A1 * (300.e0 / TEMP)**B1
    xpo = ko * M / kinf
    MZ4_TROE_TMP = ko / (1. + xpo)
    xpo = log10(xpo)
    xpo = 1. / (1. + xpo*xpo)
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
    return MZ4_USR2() * 3.333e26 * exp(-10990.e0/TEMP)


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
    # REAL(kind=dp) KO, TINV

    tinv = 1/TEMP
    ko = M * 6.5e-34 * exp(1335.*tinv)
    ko = ko / (1. + ko/(2.7e-17*exp(2199.*tinv)))
    return ko + 2.4e-14*exp(460.*tinv)


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
    return MZ4_USR6() * exp(-10900./TEMP) / 2.1e-27


def MZ4_USR8():
    """
    USR8 reaction rate as defined in MOZART4

    Returns 1.5e-13 * (1. + 6.e-7 * boltz * M * TEMP)
    """
    # real, parameter ::  boltz    = 1.38044e-16      ! erg/k
    return 1.5e-13 * (1. + 6.e-7 * 1.38044e-16 * M * TEMP)


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
    # REAL(kind = dp) ko, kinf, fc, tinv
    tinv = 1.0/TEMP
    ko = 2.3e-13 * exp(600.0*tinv)
    kinf = 1.7e-33 * M * exp(1000.0*tinv)
    fc = 1.0 + 1.4e-21 * H2O * exp(2200.0*tinv)

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
    return MZ4_USR11() * 1.111e28 * exp(-14000.0 / TEMP)


def MZ4_USR14():
    """
    USR14 reaction rate as defined in MOZART4

    Return 1.1e-11 * 300.e0/ TEMP / M
    """
    return 1.1e-11 * 300.e0 / TEMP / M


def MZ4_USR15():
    """
    USR15 reaction rate as defined in MOZART4

    Returns MZ4_USR14() * 1.111e28 *  exp( -14000.e0 / TEMP )
    """
    return MZ4_USR14() * 1.111e28 * exp(-14000.e0 / TEMP)


def MZ4_USR21():
    """
    USR21 reaction rate as defined in MOZART4

    Returns TEMP**2 * 7.69e-17 * exp( 253.e0/TEMP )
    """
    return TEMP**2 * 7.69e-17 * exp(253.e0/TEMP)


def MZ4_USR22():
    """
    USR22 reaction rate as defined in MOZART4

    Returns 3.82e-11 * exp( -2000.0/TEMP ) + 1.33e-13
    """
    return 3.82e-11 * exp(-2000.0/TEMP) + 1.33e-13


def MZ4_USR23():
    """
    USR23 reaction rate as defined in MOZART4

    fc = 3.0e-31 *(300.0/TEMP)**3.3e0
    ko = fc * M / (1.0 + fc * M / 1.5e-12)
    Returns ko * .6e0**(1. + (log10(fc * M / 1.5e-12))**2.0)**(-1.0)
    """
    # REAL(kind=dp) ko, fc

    fc = 3.0e-31 * (300.0/TEMP)**3.3e0
    ko = fc * M / (1.0 + fc * M / 1.5e-12)
    return ko * .6e0**(1. + (log10(fc * M / 1.5e-12))**2.0)**(-1.0)


def MZ4_USR24():
    """
    USR24 reaction rate as defined in MOZART4

    #REAL(kind=dp) ko
    ko = 1.0 + 5.5e-31 * exp( 7460.0/TEMP ) * M * 0.21e0
    return 1.7e-42 * exp( 7810.0/TEMP ) * M * 0.21e0 / ko
    """

    # REAL(kind=dp) ko
    ko = 1.0 + 5.5e-31 * exp(7460.0/TEMP) * M * 0.21e0
    return 1.7e-42 * exp(7810.0/TEMP) * M * 0.21e0 / ko
