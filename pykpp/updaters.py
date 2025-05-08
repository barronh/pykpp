from __future__ import print_function
__all__ = [
    'Monitor', 'interp_updater', 'spline_updater', 'Update_RCONST',
    'Update_SUN', 'Update_THETA', 'Update_M', 'interpolated_from_csv',
    'splined_from_csv', 'add_time_interpolated',
    'add_time_interpolated_from_csv', 'code_updater', 'add_code_updater',
    'func_updater', 'solar_declination'
]

import numpy as np
from scipy import constants
from scipy.constants import *
from scipy.constants import centi, R, Avogadro, kilo, Boltzmann, nano
from numpy import cos, sin, degrees, radians, arccos, pi, interp, inf, array
_np = np  # to avoid flake8 F401
_constants = constants  # to avoid flake8 F401


class updater:
    def __init__(self, incr=600, allowforce=True, verbose=0, **props):
        """
        Initialize updater. updatenow method is used to determine if updater
        should be called. If True, __call__ updates the world dictionary
        based on properties.

        Arguments
        ---------
        incr : float
            Elapsed time after which updater should be called. Defaults to 600s
        allowforce : bool
            If True, return True even if elapsed time is too short unless
            the class does not allowforce.
        verbose : int
            Level of verbosity
        props : mappable
            Addition properties provided.
        """
        self.reset()
        self.verbose = verbose
        self.incr = incr
        self.props = props
        self.allowforce = allowforce

    def _call(self, mech, world):
        if self.verbose > 0:
            print('calling updater')
        world.update(**self.props)

    def __call__(self, mech, world, force=False):
        """
        Arguments
        ---------
        mech : pykpp.mech.Mech
            Mechanism to be used in the context of world
        world : mappable
            Dictionary of system state
        force : bool
            If True, call updater even if unecessary.

        Returns
        -------
        update : bool
            True if updated, False otherwise
        """
        t = world['t']
        update = self.updatenow(t, force=force)
        if update:
            self._call(mech, world)

        return update

    def reset(self):
        """
        Set last to -inf so that next updatenow will return True.

        Returns
        -------
        None
        """
        self.last = -inf

    def updatenow(self, t, force=False):
        """
        test if incr has passed since last update

        Arguments
        ---------
        t : float
            Time (typically s) in system.
        force : bool
            If True, return True even if elapsed time is too short unless
            the class does not allowforce.

        Returns
        -------
        update : bool
            If True, the time since last called is large enough that the
            instance should be called.
        """
        tsince = abs(t - self.last)
        if tsince > self.incr or (force and self.allowforce):
            self.last = t
            return True
        else:
            return False


class interp_updater(updater):
    def __init__(self, incr=600, allowforce=True, verbose=False, time=None, **props):
        """
        Initialize updater. updatenow method is used to determine if updater
        should be called. If True, __call__ updates the world dictionary
        based on properties.

        Arguments
        ---------
        incr : float
            Elapsed time after which updater should be called.
        allowforce : bool
            If True, return True even if elapsed time is too short unless
            the class does not allowforce.
        verbose : int
            Level of verbosity
        time : array
            Time (typically s) in system.
        props : mappable
            Addition property arrays with the same number of elements as time.
        """
        self.reset()
        self.verbose = verbose
        if time is None:
            for tkey in ['t', 'time', 'Time']:
                if tkey in props:
                    time = props.pop(tkey)
        self.time = time
        self.incr = incr
        self.props = props
        self.allowforce = allowforce

    def _call(self, mech, world):
        """
        update world dictionary with keywords interpolated with respect to
        time from props supplied during initialization.

        Arguments
        ---------
        mech : pykpp.mech.Mech
            Mechanism to be used in the context of world
        world : mappable
            Dictionary of system state

        Returns
        -------
        None
        """
        if self.verbose:
            print(
                "Updating %s: %s"
                % (', '.join(self.props.keys()), self.last)
            )
        t = world['t']
        for k, vs in self.props.items():
            world[k] = interp(t, self.time, vs)


class spline_updater(updater):
    def __init__(self, incr=600, allowforce=True, verbose=False, time=None, **props):
        """
        Initialize updater. updatenow method is used to determine if updater
        should be called. If True, __call__ updates the world dictionary
        based on properties.

        Arguments
        ---------
        incr : float
            Elapsed time after which updater should be called.
        allowforce : bool
            If True, return True even if elapsed time is too short unless
            the class does not allowforce.
        verbose : int
            Level of verbosity
        time : array
            Time (typically s) in system.
        props : mappable
            Addition property arrays with the same number of elements as time.
        """
        from scipy.interpolate import splrep
        self.reset()
        self.verbose = verbose
        if time is None:
            for tkey in ['t', 'time', 'Time']:
                if tkey in props:
                    time = props.pop(tkey)

        self.time = time
        self.incr = incr
        self.splreps = dict([
            (k, splrep(self.time, v, s=0)) for k, v in props.items()
        ])
        self.allowforce = allowforce

    def _call(self, mech, world):
        """
        update world dictionary with keywords interpolated using a spline
        with respect to time from props supplied during initialization.

        Arguments
        ---------
        mech : pykpp.mech.Mech
            Mechanism to be used in the context of world
        world : mappable
            Dictionary of system state

        Returns
        -------
        None
        """
        from scipy.interpolate import splev
        if self.verbose:
            print(
                "Updating %s: %s"
                % (', '.join(self.props.keys()), self.last)
            )
        t = world['t']
        for k, vs in self.splreps.items():
            world[k] = splev(t, vs, der=0)


class func_updater(updater):
    def __init__(self, func, incr=600, allowforce=True, verbose=False):
        """
        Initialize updater. updatenow method is used to determine if updater
        should be called. If True, __call__ updates the world dictionary
        based on properties.

        Arguments
        ---------
        func : callable
            Function to all after time elaped.
        allowforce : bool
            If True, return True even if elapsed time is too short unless
            the class does not allowforce.
        verbose : int
            Level of verbosity
        time : array
            Time (typically s) in system.
        props : mappable
            Addition property arrays with the same number of elements as time.
        """
        self.reset()
        self.verbose = verbose
        self.incr = incr
        self._call = func
        self.allowforce = allowforce


class code_updater(updater):
    def __init__(
        self, code, incr=600, allowforce=True, verbose=False, message='code'
    ):
        """
        code - string that can be compiled as exec
        incr - frequency to re-execute code
        verbose - show update status
        message - indentify this updater as message
        """
        self.verbose = verbose
        self.block = compile(code, '<user>', 'exec')
        self.incr = incr
        self.message = message
        self.reset()
        self.allowforce = allowforce

    def _call(self, mech, world):
        """
        if time increment incr has passed or force, exec code from
        __init__ with world as locals dictionary

        Arguments
        ---------
        mech : pykpp.mech.Mech
            Mechanism to be used in the context of world
        world : mappable
            Dictionary of system state

        Returns
        -------
        None
        """
        if self.verbose:
            print("Updating %s: %s" % (self.message, self.last))

        exec(self.block, globals(), world)


def add_time_interpolated(time, incr=0, verbose=False, **props):
    """
    Shortcut to
        add_world_updater(
            interp_updater(time=time, incr=incr, verbose=verbose, **props)
        )
    """
    global add_world_updater
    add_world_updater(
        interp_updater(time=time, incr=incr, verbose=verbose, **props)
    )


def interpolated_from_csv(path, timekey, incr=0, delimiter=',', verbose=False):
    import pandas as pd

    data = pd.read_csv(path, delimiter=delimiter)
    datadict = dict([(k, data[k]) for k in data.keys()])
    time = datadict.pop(timekey)
    return interp_updater(time=time, incr=incr, verbose=verbose, **datadict)


def splined_from_csv(path, timekey, incr=0, delimiter=',', verbose=False):
    import pandas as pd

    data = pd.read_csv(path, delimiter=delimiter)
    datadict = dict([(k, data[k]) for k in data.keys()])
    time = datadict.pop(timekey)
    return spline_updater(time=time, incr=incr, verbose=verbose, **datadict)


def add_time_interpolated_from_csv(path, timekey, incr=0):
    """
    Shortcut to add_time_interpolated from data in a csv file
    """
    global add_world_updater
    add_world_updater(interpolated_from_csv(path, timekey, incr))


def add_code_updater(code, incr=0, verbose=False, message='code'):
    """
    Shortcut to:
    add_world_updater(
        code_updater(code=code, incr=incr, verbose=verbose, message=message)
    )
    """
    global add_world_updater
    add_world_updater(
        code_updater(code=code, incr=incr, verbose=verbose, message=message)
    )


# def add_world_updater(func, incr=0, verbose=False):
#     """
#     Add func to be called with mech and world
#     to update the world environment
#     """
#     Update_World.add(func_updater(func=func, incr=incr, verbose=verbose))

# add_func_updater = add_world_updater

def Monitor(mech, world=None):
    """
    Call the print_monitor(y, t) where t is in world and y = mech.get_y()

    Arguments
    ---------
    mech : pykpp.mech.Mech
        Mechanism with properties relevant to the world
    world : dict
        Dictionary with all properties of the world. Concentrations are
        typically in molecules/cm3, but the only requirement is that they
        match the units of rates.

    Returns
    -------
    None
    """
    try:
        t = eval('t', None, world)
        y = mech.get_y()
        mech.print_monitor(y, t)
    except Exception as e:
        print(str(e))
        pass


def Update_RCONST(mech, world=None):
    """
    Evaluate rconst in the context of world.

    Arguments
    ---------
    mech : pykpp.mech.Mech
        Mechanism with properties relevant to the world
    world : dict
        Dictionary with all properties of the world. Concentrations are
        typically in molecules/cm3, but the only requirement is that they
        match the units of rates.

    Returns
    -------
    None
    """
    for rconst in mech._parsed['RCONST']:
        exec(rconst, None, world)


def Update_SUN(mech, world):
    """
    Updates world dectionary to contain
      SUN - scaling variable between 0 and 1 (following Sandu et al.)

    if t < SunRise or t > SunSet: SUN = 0.

    hour = time since noon
    squared = abs(hour) * hour

    Arguments
    ---------
    mech : pykpp.mech.Mech
        Mechanism with properties relevant to the world
    world : dict
        Dictionary with all properties of the world. Concentrations are
        typically in molecules/cm3, but the only requirement is that they
        match the units of rates.

    Returns
    -------
    None
    """
    t = world['t']
    SunRise = world['SunRise']
    SunSet = world['SunSet']
    Thour = t / 3600.0
    Tlocal = Thour % 24.

    if (Tlocal >= SunRise) and (Tlocal <= SunSet):
        Ttmp = (2.0 * Tlocal - SunRise - SunSet) / (SunSet - SunRise)
        if (Ttmp > 0.):
            Ttmp = Ttmp * Ttmp
        else:
            Ttmp = -Ttmp * Ttmp

        SUN = (1.0 + cos(pi * Ttmp)) / 2.0
    else:
        SUN = 0.0
    world['SUN'] = SUN


def solar_declination(N):
    """
    N - julian day 1-365 (1 = Jan 1; 365 = Dec 31)
    Returns solar declination in radians

    wikipedia.org/wiki/Declination_of_the_Sun
    dec_deg = -23.44 * cos_deg(360./365 * (N + 10))
    dec_rad = (pi / 180. * -23.44) * cos_rad(pi / 180. * 360./365 * (N + 10))
    """
    return -0.40910517666747087 * cos(0.017214206321039961 * (N + 10.))


def solar_zenith_angle(t, dec=None, jday=1, phi=0.6981317007977318):
    """
    Arguments
    ---------
    t : float
        Time in seconds in local solar time (time_utc + lon / 15)
    dec : float
        Declination angle in radians
    jday : int
        Julian date (days since Jan 1)
    phi : float
        Latitude in radians (defautls to radians(40 degrees))

    Returns
    -------
    sza : float
        Solar Zenith Angle (angle from zenith) in radians
    """
    if dec is None:
        dec = solar_declination(jday)

    hlocal = (t / 3600.) % 24.
    houra = radians((hlocal - 12.) * 15.)
    sza = arccos(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(houra))
    return sza


def Update_THETA(mech, world):
    """
    Adds solar zenith angle (THETA; angle from solar noon) in degrees
    to the world dictionary based on time (t), and latitude
    (Latitude_Radians) and StartJday or SolarDeclination_Radians

    THETA = arccos(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(houra))

    Arguments
    ---------
    mech : pykpp.mech.Mech
        Mechanism with properties relevant to the world
    world : dict
        Dictionary with all properties of the world. Concentrations are
        typically in molecules/cm3, but the only requirement is that they
        match the units of rates.

    Returns
    -------
    None
    """
    phi = world['Latitude_Radians']
    t = world['t']
    if 'SolarDeclination_Radians' in world:
        dec = world['SolarDeclination_Radians']
        THETA = solar_zenith_angle(t, dec=dec, phi=phi)
    else:
        StartJday = world['StartJday']
        N = StartJday + (t / 3600.) // 24
        THETA = solar_zenith_angle(t, jday=N, phi=phi)

    world['THETA'] = degrees(THETA)


def Update_M(mech, world):
    """
    Adds concentrations (molecules/cm3) to world namespace for:
        M (air) = P / R / TEMP * Avogadro * 1e-6,
        O2 (0.20946*M),
        N2 (0.78084*M), and
        H2 (550 ppb)
    based on:
        Pressure (P [=] Pascals),
        Temperature (TEMP [=] Kelvin), and
        R is provided [=] m**3 * Pascals/K/mol

    TEMP and P must be defined either in world or in stdfuncs

    Arguments
    ---------
    mech : pykpp.mech.Mech
        Mechanism with properties relevant to the world
    world : dict
        Dictionary with all properties of the world. Concentrations are
        typically in molecules/cm3, but the only requirement is that they
        match the units of rates.

    Returns
    -------
    None
    """
    from warnings import warn
    try:
        M = eval('P / R / TEMP * Avogadro * 1e-6', None, world)
    except Exception as e:
        warn('Cannot eval M; requires P and TEMP;' + str(e))
        M = float(world['M'])
    world['M'] = M
    world['O2'] = 0.20946 * M
    world['N2'] = 0.78084 * M
    world['H2'] = 0.00000055 * M
