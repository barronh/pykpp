def test_sza():
    """
    Testing the actual values returned for solar zenith angle
    based on a reference method maintained by NOAA GML.
    Results bbased on May 9 2025 at 40N 0E

    https://gml.noaa.gov/grad/solcalc/azel.html
    """
    import numpy as np
    from ..updaters import solar_zenith_angle, Update_THETA
    decs = np.radians(17 + np.array([.45, .48, .5, .52, .54, .56]))
    ts = np.array([8, 10, 12, 14, 16, 18]) * 3600.
    szachks = np.radians(
        90 - np.array([34.63, 56.25, 67.49, 55.13, 33.33, 10.6])
    )
    world = {'StartJday': 129, 'Latitude_Radians': np.radians(40.)}
    for t, dec, chk in zip(ts, decs, szachks):
        world['t'] = t
        world['SolarDeclination_Radians'] = dec
        Update_THETA(None, world)
        sza = solar_zenith_angle(t, dec=dec, jday=129)
        fracdiff = abs(1 - sza / chk)
        # print(t / 3600, np.degrees(sza), np.degrees(chk), fracdiff)
        assert fracdiff < 0.02  # error is less than 2%
        assert world['THETA'] == np.degrees(sza)
    

def test_funcupdater():
    """
    Mostly testing that this works. The function being tested is super simple
    and arbitrary.
    """
    import numpy as np
    from ..updaters import func_updater
    def sqt(mech, world):
        world['t2'] = world['t']**2

    fu = func_updater(sqt, incr=600)
    ts = [0, 1200, 1500, 2400]
    chks = [0, 1200**2, 1200**2, 2400**2]
    world = {}
    for t, chk in zip(ts, chks):
        print(t, chk)
        world['t'] = t
        fu(None, world)
        assert np.allclose(world['t2'], chk)


def test_codeupdater():
    """
    Mostly testing that this works. The code being tested is super simple and
    arbitrary.
    """
    import numpy as np
    from ..updaters import code_updater
    code = 't2 = t**2'

    fu = code_updater(code, incr=600)
    ts = [0, 1200, 1500, 2400]
    chks = [0, 1200**2, 1200**2, 2400**2]
    world = {}
    for t, chk in zip(ts, chks):
        print(t, chk)
        world['t'] = t
        fu(None, world)
        assert np.allclose(world['t2'], chk)


def test_interps():
    """
    Testing interpolation and spline methods. Tests exercise the functionality
    of both interpolation and incr time checking.
    """
    from tempfile import TemporaryDirectory as TmpDir
    import numpy as np
    from ..updaters import interpolated_from_csv, splined_from_csv
    with TmpDir() as td:
        tfpath = f'{td}/test.csv'
        with open(tfpath, 'w') as tf:
            tf.write("""time	PBLH	TEMP
0.00	250.	288.15
14400	250.	288.15
18000	300.	288.16
21600	550.	289
""")
        liniu = interpolated_from_csv(tfpath, 'time', incr=600, delimiter='\t')
        spiu = splined_from_csv(tfpath, 'time', incr=600, delimiter='\t')
        linchecks = [
            {'t': 0., 'PBLH': 250, 'TEMP': 288.15},        # initial value
            {'t': 14300., 'PBLH': 250, 'TEMP': 288.15},    # changed, but sam value
            {'t': 14850., 'PBLH': 250, 'TEMP': 288.15},    # no change 500s
            {'t': 16200., 'PBLH': 275., 'TEMP': 288.155},  # intermediate value
            {'t': 19800., 'PBLH': 425, 'TEMP': 288.58},    # intermediate value
            {'t': 21600., 'PBLH': 550, 'TEMP': 289.},     # intermediate value
        ]
        spchecks = [
            {'t': 0.0, 'PBLH': 250., 'TEMP': 288.15},
            {'t': 14300.0, 'PBLH': 250.59767233, 'TEMP': 288.15758533},
            {'t': 14850.0, 'PBLH': 250.59767233, 'TEMP': 288.15758533},
            {'t': 16200.0, 'PBLH': 255.625, 'TEMP': 288.0770625},
            {'t': 19800.0, 'PBLH': 394.375, 'TEMP': 288.4504375},
            {'t': 21600.0, 'PBLH': 550., 'TEMP': 289.},
        ]
        world = {}
        chkkeys = ['PBLH', 'TEMP']
        iuchecks = [('lin', liniu, linchecks), ('spline', spiu, spchecks)]
        for ikey, iu, checks in iuchecks:
            for chk in checks:
                world['t'] = chk['t']
                iu(None, world)
                assert np.allclose(
                    [world[k] for k in chkkeys],
                    [chk[k] for k in chkkeys]
                )


def test_m():
    """
    Basically, reimplemented code and checking that it works. If something
    changes, we should know about it.
    """
    import numpy as np
    from ..updaters import Update_M
    chkkeys = ['M', 'O2', 'N2', 'H2']
    world = {'t': 0., 'P': 101325., 'TEMP': 298.}
    chk = {}
    Update_M(None, world)
    chk['M'] = eval('P / TEMP', None, world) / 8.314462618 * 6.02214076e+17
    chk['O2'] = 0.20946 * chk['M']
    chk['N2'] = 0.78084 * chk['M']
    chk['H2'] = 550e-9 * chk['M']
    assert np.allclose(
        [world[k] for k in chkkeys],
        [chk[k] for k in chkkeys],
    )
    