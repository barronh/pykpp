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
    def __init__(self, *args, **kwds):
        self.reset()

    def __call__(self, mech, world, force=False):
        return False

    def reset(self):
        """
        start increment testing over
        """
        self.last = -inf

    def updatenow(self, t, force=False):
        """
        test if incr has passed since last update
        """
        tsince = abs(t - self.last)
        if tsince > self.incr or (force and self.allowforce):
            self.last = t
            return True
        else:
            return False


class interp_updater(updater):
    def __init__(self, time, incr, allowforce=True, verbose=False, **props):
        """
        time - time array
        incr - frequency to re-execute code
        verbose - show update status
        props - keyword variables to update
        """
        self.reset()
        self.verbose = verbose
        self.time = time
        self.incr = incr
        self.props = props
        self.allowforce = allowforce

    def __call__(self, mech, world, force=False):
        """
        mech - mechanism object
        world - dictionary representing state
        force - update even if frequency indicates not necessary
        update world dictionary with keywords from props in __init__
        """
        t = world['t']
        update = self.updatenow(t, force=force)
        if update:
            if self.verbose:
                print(
                    "Updating %s: %s"
                    % (', '.join(self.props.keys()), self.last)
                )
            for k, vs in self.props.items():
                world[k] = interp(t, self.time, vs)

        return update


class spline_updater(updater):
    def __init__(self, time, incr, allowforce=True, verbose=False, **props):
        """
        time - time array
        incr - frequency to re-execute code
        verbose - show update status
        props - keyword variables to update
        """
        from scipy.interpolate import splrep
        self.reset()
        self.verbose = verbose
        self.time = time
        self.incr = incr
        self.splreps = dict([
            (k, splrep(self.time, v, s=0)) for k, v in props.items()
        ])
        self.allowforce = allowforce

    def __call__(self, mech, world, force=False):
        """
        mech - mechanism object
        world - dictionary representing state
        force - update even if frequency indicates not necessary
        update world dictionary with keywords from props in __init__
        """
        from scipy.interpolate import splev
        t = world['t']
        update = self.updatenow(t, force=force)
        if update:
            if self.verbose:
                print(
                    "Updating %s: %s"
                    % (', '.join(self.props.keys()), self.last)
                )
            for k, vs in self.splreps.items():
                world[k] = splev(t, vs, der=0)

        return update


class func_updater(updater):
    def __init__(self, func, incr, allowforce=True, verbose=False):
        """
        func - function that takes mech, and world
        incr - frequency to re-execute code
        verbose - show update status
        message - indentify this updater as message
        """
        self.reset()
        self.verbose = verbose
        self.incr = incr
        self.func = func
        self.allowforce = allowforce

    def __call__(self, mech, world, force=False):
        """
        mech - mechanism object
        world - dictionary representing state
        force - update even if frequency indicates not necessary

        if time increment incr has passed or force, call func from
        __init__ with mech and world as arguments
        """
        t = world['t']
        update = self.updatenow(t, force=force)
        if update:
            if self.verbose:
                print(
                    "Updating %s: %s"
                    % (str(self.func).split(' ')[1], self.last)
                )
            self.func(mech, world)

        return update


class code_updater(updater):
    def __init__(
        self, code, incr, allowforce=True, verbose=False, message='code'
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

    def __call__(self, mech, world, force=False):
        """
        mech - mechanism object
        world - dictionary representing state
        force - update even if frequency indicates not necessary

        if time increment incr has passed or force, exec code from
        __init__ with world as locals dictionary
        """
        last = self.last
        t = world['t']
        update = self.updatenow(t, force=force)
        if update:
            if self.verbose:
                print("Updating %s: %s" % (self.message, last))

            exec(self.block, globals(), world)

        return update


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
    try:
        t = eval('t', None, world)
        y = mech.get_y()
        mech.print_monitor(y, t)
    except Exception as e:
        print(str(e))
        pass


def Update_RCONST(mech, world=None):
    for rconst in mech._parsed['RCONST']:
        exec(rconst, None, world)


def Update_SUN(mech, world):
    """
    Updates world dectionary to contain
      SUN - scaling variable between 0 and 1 (following Sandu et al.)

    if t < SunRise or t > SunSet: SUN = 0.

    hour = time since noon
    squared = abs(hour) * hour

    """
    t = world['t']
    SunRise = world['SunRise']
    SunSet = world['SunSet']
    Thour = t/3600.0
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


def Update_THETA(mech, world):
    """
    Adds solar zenith angle (THETA; angle from solar noon) in degrees
    to the world dictionary based on time

    THETA = arccos(sin(lat) * sin(dec) + cos(lat) * cos(dec) * cos(houra))
    """
    phi = world['Latitude_Radians']
    t = world['t']
    if 'SolarDeclination_Radians' in world:
        dec = world['SolarDeclination_Radians']
    else:
        StartJday = world['StartJday']
        N = StartJday + (t / 3600.) // 24
        dec = solar_declination(N)
    Tlocal = (t / 3600.) % 24.
    houra = radians((Tlocal - 12.) * 15.)
    THETA = arccos(sin(phi) * sin(dec) + cos(phi) * cos(dec) * cos(houra))
    world['THETA'] = degrees(THETA)


def Update_M(mech, world):
    """
    Adds concentrations (molecules/cm3) to world namespace for:
        M (air),
        O2 (0.20946 M),
        N2 (0.78084 M), and
        H2 (500 ppb)
    based on:
        Pressure (P in Pascals),
        Temperature (TEMP in Kelvin), and
        R is provided in m**3 * Pascals/K/mol

    TEMP and P must be defined either in world or in stdfuncs
    """
    from warnings import warn
    try:
        M = eval('P / (R / centi**3) / TEMP * Avogadro', None, world)
    except Exception as e:
        warn('Cannot eval M; requires P and TEMP')
        M = float(world['M'])
    world['M'] = M
    world['O2'] = 0.20946 * M
    world['N2'] = 0.78084 * M
    world['H2'] = 0.00000055 * M
