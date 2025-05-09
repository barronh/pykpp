"""
Advanced Trajectory Simulation
------------------------------

The previous example used a reference set of inputs to create a reasonable
summer day. This can then be adapted to make a specific summer day. This
can be done using a single trajectory or an ensemble of trajectories.

Prerequisites
^^^^^^^^^^^^^

    python -m pip install pyrsig netcdf4 pycno pyproj

Trajectory Environment
^^^^^^^^^^^^^^^^^^^^^^

You can replace teh "Summer" PBL and temperature using HYSPLIT. HYSPLIT runs
back or forward trajectories and allows you to save out mixing depth and temp.
The output is a text dump (tdump) file that you can eaisly turn into a csv
like the one used here. Similarl

Trajectory Emissions
^^^^^^^^^^^^^^^^^^^^

The emissions can be replaced with new values several ways. For example, EPA
routinely puts gridded annual emissions on their ftp site. And, the 2022
modeling platform is on AWS. You can sample the emission inventory along the
trajectory environment (lon/lat). For the full modeling platform, you will
have to make some decisions about what emissions to include and how. The
gridded emissions are easy because they can be assumed to be emitted within
the mixed depth. It is more complicated to decide how to treat point emissions.
The files do not include plume rise, nor does it tell you how to treat
horizontal dillution. This tutorial will not answer those questions, because
it will use the annual average reports that are at 12-km resolution with no
vertical resolution.
"""

# %%
# Part 1: Trajectory Environment Offline
# --------------------------------------
#
# 1. Go to https://www.ready.noaa.gov/HYSPLIT_traj.php
# 2. Choose "compute archive Trajectories"
# 3. Choose 1 "Starting Locations" in "Normal" mode; next
# 4. Choose HRRR 3km (sigma, U.S., 06/2019-present); next
# 5. Choose the JLG Super Site coordiantes 33.503833 N and 112.095767 W; next
# 6. Choose 20240722_18-23_hrrr; next
# 7. Choose "Backward"
# 8. Choose hour 21Z (peak ozone 86 ppb)
# 9. Enable all "dump meteorology" options (Terrain Height, Ambient Temperature, Mixed Depth,e t)
# 10. Run trajectory and wait
# 11. When complete, download the tdump file save here as tdump.inp

# %%
# Read Environment
# ^^^^^^^^^^^^^^^^
#
import pandas as pd

std = ['tid', 'gid', 'yy', 'mm', 'dd', 'HH', 'MM', 'fhr', 'age', 'lat', 'lon', 'alt']
trajpath = 'trajdump.txt'
with open(trajpath, 'r') as tf:
    for _l in tf:
        if 'PRESSURE' in _l:
            break

    keys = std + _l.split()[1:]
    trajdf = pd.read_csv(tf, names=keys, sep=r'\s+')

trajdf = trajdf[[k for k in trajdf.columns if k != 'THETA']]
ttime = trajdf.eval('yy * 1e8 + mm * 1e6 + dd * 1e4 + HH * 1e2 + MM')  # calc time
trajdf['TSTEP'] = pd.to_datetime(ttime, format='%y%m%d%H%M') # add time coord as TSTEP
trajdf['time_lst'] = trajdf['TSTEP'] + pd.to_timedelta(-112 / 15., unit='h')
trajdf['t'] = (trajdf['time_lst'] - trajdf['time_lst'].min().floor('1d')).dt.total_seconds()


# %%
# Download Emissions
# ^^^^^^^^^^^^^^^^^^
#

import os
import requests
import pyrsig

# files are Mg/yr and need to be mole/s...
root = 'https://gaftp.epa.gov/Air/emismod/2022/v1/gridded'
# if you want sectors
# anthpath = '2022hc_cb6_22m_total_nobeis_nofires_withfert_12US1_annual.ncf'
allpath = '2022hc_cb6_22m_total_withbeis_withfires_withfert_12US1_annual.ncf'

if not os.path.exists(allpath):
    with requests.get(f'{root}/{allpath}', verify=False) as r, open(allpath, 'wb') as out:
        out.write(r.content)

ef = pyrsig.open_ioapi(allpath).drop_vars('TFLAG')

# %%
# Get Emissions Along Path
# ^^^^^^^^^^^^^^^^^^^^^^^^
import pyproj
import copy

proj = pyproj.Proj(ef.crs_proj4)
trajdf['COL'], trajdf['ROW'] = proj(trajdf.lon, trajdf.lat)
trajdf['LAY'] = 0
t = trajdf['TSTEP'].to_xarray()  # input is time invariant
k = trajdf['LAY'].to_xarray()
r = trajdf['ROW'].to_xarray()
c = trajdf['COL'].to_xarray()

sopt = dict(TSTEP=t, LAY=k, ROW=r, COL=c, method='nearest')
erdf = ef.sel(**sopt).to_dataframe().drop(ef.dims, axis=1)

# Convert gases from mole/s to molecules/cm3/s using mixdepth
V = trajdf['MIXDEPTH'] * ef.XCELL * ef.YCELL * 1e6  # cm**3 [=] m * m * m * cm**3 / m**3
# molec/s [=] g/ton * y/d * d/h * h/s * mol/g * molec / mol
unitfactor = 6.022e23 * 907185 / 365 / 24 / 3600 / pd.Series({
    'NO2': 46.0, 'NO': 30.0, 'HONO': 47.0,
    'SO2': 64.0, 'PACD': 76.0, 'AACD': 60.0, 'ALD2': 44.0, 'FORM': 30.0,
    'MEOH': 32.0, 'FACD': 46.0, 'CO': 28.0, 'ALDX': 58.1, 'GLYD': 60.0,
    'GLY': 58.0, 'MGLY': 72.0, 'ETHA': 30.1, 'ETOH': 46.1, 'KET': 72.1,
    'PAR': 14.0, 'ACET': 58.1, 'PRPA': 44.1, 'ETHY': 26.0, 'ETH': 28.0,
    'OLE': 42.1, 'IOLE': 56.1, 'ISOP': 68.1, 'ISPD': 70.1, 'TERP': 136.2,
    'APIN': 136.2, 'TOL': 92.1, 'XYLMN': 106.2, 'NAPH': 128.2, 'CL2': 71.0,
    'HCL': 36.5, 'SESQ': 204.0, 'BUTADIENE13': 54.0, 'ACROLEIN': 56.1
})

edf = erdf[list(unitfactor.index)].multiply(unitfactor).divide(V, axis=0)

emis = '\n'.join([
    f'<E{k}> EMISSION = {k} : emis_{k} ;'
    for k in edf.columns
])

# deposition velocities
vd = dict(
    CO=0.03,
    SO2=1.8,
    NO=0.016,
    NO2=0.1,
    NO3=0.1,
    N2O5=4,
    HONO=2.2,
    HNO3=4,
    HNO4=4,
    O3=0.4,
    H2O2=0.5,
    CH3OOH=0.1,
)
vdf = pd.DataFrame({
    k: v / (trajdf['MIXDEPTH'] * 100)  # convert from cm/s to 1/s
    for k, v in vd.items()
})

depn = """
<DDEPCO> CO = DUMMY : vd_CO ;
<DDEPSO2> SO2 = DUMMY : vd_SO2 ;
<DDEPNO> NO = DUMMY : vd_NO ;
<DDEPNO2> NO2 = DUMMY : vd_NO2 ;
<DDEPNO3> NO3 = DUMMY : vd_NO3 ;
<DDEPN2O5> N2O5 = DUMMY : vd_N2O5 ;
<DDEPHONO> HONO = DUMMY : vd_HONO ;
<DDEPHNO3> HNO3 = DUMMY : vd_HNO3 ;
<DDEPHNO4> PNA = DUMMY : vd_HNO4 ;
<DDEPO3> O3 = DUMMY : vd_O3 ;
<DDEPH2O2> H2O2 = DUMMY : vd_H2O2 ;
<DDEPCH3OOH> MEPX = DUMMY : vd_CH3OOH ;
"""

trajedf = trajdf.join(
    edf.rename(columns=lambda k: 'emis_' + k)
).join(vdf.rename(columns=lambda k: 'vd_' + k))
trajedf['P'] = trajedf['PRESSURE'] * 100
trajedf['TEMP'] = trajedf['AIR_TEMP']


# %%
# Create Mechanism
# 
from pykpp.mech import Mech
from pykpp.updaters import Update_M, Update_THETA
from pykpp.stdfuncs import update_func_world
from pykpp.tuv.tuv5pt0 import TUV_J5pt0
import numpy as np

mechdef = """
#INLINE PY_INIT
TEMP = 298.15
P = 99600.
t = TSTART = 0 * 3600.
TEND = TSTART + 3600. * 24
DT = 600.
MONITOR_DT = 3600000
StartDate = 'datetime(2024, 7, 22)'
Latitude_Degrees = 33.503833
Longitude_Degrees = -112.095767
#ENDINLINE

#INITVALUES
CFACTOR = P * Avogadro / R / TEMP * centi **3 * nano {ppb-to-molecules/cm3}
ALL_SPEC=1e-32*CFACTOR;
TOTALNOx=1.
M=1e9
O2=.21*M
N2=.79*M
H2O=0.01*M
CH4=1750
CO=100
O3=60.
NO = 0.1 * TOTALNOx
NO2 = 0.9 * TOTALNOx
SO2 = 1
N2O = 320
{B = 210.; what to do about B?}

#LOOKAT ALL;

#include cb6r5m_ae7_aq.eqn
""" + f"""
{depn}

{emis}
"""
mech = Mech(
    mechdef, mechname='test', timeunit='utc', monitor_incr=None,
    add_default_funcs=True
)

dt = mech.world['DT']
ts = np.arange(trajedf['t'].min(), trajedf['t'].max(), dt)
envdf = trajedf.set_index('t').to_xarray().interp(t=ts, method='linear').to_dataframe()
rows = []
olddepth = envdf.iloc[0]['MIXDEPTH']
fast = ('O', 'O1D', 'HCO3', 'NTR1')
mech.world['atol'] = np.array([10 if k in fast else 1e-3 for k in mech.allspcs])
mech.world['rtol'] = np.array([0.1 if k in fast else 1e-5 for k in mech.allspcs])
ybkg = mech.get_y().copy()
for t in ts:
    mech.world['t'] = t
    mech.world.update(envdf.loc[t].to_dict())
    if olddepth < mech.world["MIXDEPTH"]:
        fold = olddepth / mech.world["MIXDEPTH"]
        print(fold)
        mech.world['O3'] = fold * mech.world['O3'] + (1 - fold) * 60 * mech.world['CFACTOR']
        mech.world['CO'] = fold * mech.world['CO'] + (1 - fold) * 100 * mech.world['CFACTOR']
        olddepth = mech.world["MIXDEPTH"]
    mech.run(tstart=t, tend=t + dt)
    CFACTOR = mech.world['CFACTOR']
    row = {k: mech.world[k] / CFACTOR for k in mech.allspcs}
    for k in ['t', 'CFACTOR', 'TEMP', 'P', 'MIXDEPTH', 'THETA', 'emis_NO', 'SUN_FLUX']:
        row[k] = mech.world[k]
    row['jNO2'] = TUV_J5pt0('NO2 -> NO + O(3P)', mech.world['THETA'])
    rows.append(row)

# %%
# Plot
# ^^^^
#
import matplotlib.pyplot as plt

outdf = pd.DataFrame.from_dict(rows)

outdf['h_LST'] = outdf['t'] / 3600
noxspcs = ['NO', 'NO2', 'N2O5', 'HONO', 'NO3']
outdf['NOx'] = sum([outdf[ns] for ns in noxspcs if ns in outdf.columns])
outdf['TEMP/25 [C]'] = (outdf['TEMP'] - 273.15) / 25.
outdf['OH/10 [ppqv]'] = outdf['OH'] * 1e5
outdf['HO2 [pptv]'] = outdf['HO2'] * 1e3
outdf['PBLH [km]'] = outdf['MIXDEPTH'] * 1e-3
outdf['NO_molps'] = outdf.eval('emis_NO / 6.022e23 * MIXDEPTH * 1e6 ') * ef.XCELL * ef.YCELL
outdf['jNO2*100'] = outdf['jNO2'] * 100

fig, ax = plt.subplots()
outdf.set_index('h_LST')[['NO_molps', 'jNO2*100', 'PBLH [km]', 'TEMP/25 [C]']].plot(ax=ax)
fig.savefig('hysplit_physical.png')

fig, ax = plt.subplots()
outdf.set_index('h_LST')[['NOx'] + noxspcs].plot(ax=ax)
ax.set(xlabel='hour [LST]', ylabel='ppb', ylim=(.1, 5), yscale='log')
fig.savefig('hysplit_nox.png')

fig, ax = plt.subplots()
outdf.set_index('h_LST')[['O3', 'OH/10 [ppqv]', 'HO2 [pptv]']].plot(ax=ax)
ax.set(xlabel='hour (LST)')
fig.savefig('hysplit_ox.png')

outdf.to_csv('hysplit.csv')