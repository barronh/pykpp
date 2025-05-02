"""
Upper Tropospheric Simulation
=============================

    goals: Predict ozone in the Upper Troposphere
    author: Barron Henderson
    last updated: 2025-05-02

Henderson et al. (2011) used an ensemble of simulation to create a syntheic
upper troposphere. Each simulation was initialized by observed composition from
a set of measurements identified as "fresh convection." The simulation the aged
for 10 days. The ensemble was then sampled to mimic the effect of convection.

This examples shows how to use a single simulation based on the median of
observations from Table 2.

Henderson et al. (2011) doi:10.5194/acp-11-275-2011
"""

import time
import numpy as np
from pykpp.mech import Mech
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
from pykpp.updaters import func_updater

mdnobs_ppb = pd.read_csv(StringIO("""spc,bkg,init
HO,0.5396e-3,0.6101e-3
HO2,13.16e-3,11.24e-3
O3,77.76,70.61
NO2,95.52e-3,153.6e-3
NO,203.3e-3,411.8e-3
HNO3,280.1e-3,125.9e-3
HO2NO2,82.e-3,67.8e-3
H2O2,234.2e-3,195.9e-3
CO,98.36,108.0
CH4,1.789e3,1.784e3
C2H6,790.e-3,800e-3
C3H8,146e-3,153.5e-3
C2H4,1.5e-3,1.5e-3
RNO3,8.63e-3,8.63e-3
CH2O,174.5e-3,437e-3
CH3CHO,83.8e-3,117.5e-3
CH3COCH3,1475e-3,1375e-3
CH3COC2H5,71.25e-3,95.e-3
CH3CONO2,374.9e-3,370.6e-3
CH3COOOH,172.8e-3,226.1e-3
"""))

mdnobs_ppb = mdnobs_ppb.set_index('spc').rename(index={
    'HO': 'OH', 'HO2NO2': 'HNO4', 'C2H6': 'ETHA', 'C3H8': 'PRPA',
    'C2H4': 'ETH', 'RNO3': 'NTR1', 'CH2O': 'FORM', 'CH3CHO': 'ALD2',
    'CH3COCH3': 'ACET', 'CH3COC2H5': 'KET', 'CH3CONO2': 'PAN',
    'CH3COOOH': 'PAA'
}).groupby('spc').sum()
print(list(mdnobs_ppb.index))
inittxt = '\n'.join([f'{k} = {v}' for k, v in mdnobs_ppb.init.items()])

infile = StringIO("""
#INLINE PY_INIT
TEMP = 233.7 { median upper trop temp }
P = 30060. { median upper trop pressure }
t = TSTART = 12 * 3600. { noon start}
TEND = TSTART + 3600. * 24 * 10 { run for 10 days }
DT = 60. { solve at 60s }
MONITOR_DT = 3600 { save at 1h }
StartDate = 'datetime(2010, 7, 14)'
Latitude_Degrees = 40. { latitude is nominally mid-latitudes }
Longitude_Degrees = 0.00E+00 { longitude is nominal so LST = UTC }

add_world_updater(func_updater(Update_M, incr=360., verbose=False))
add_world_updater(func_updater(Update_THETA, incr=360., verbose=False))
add_world_updater(func_updater(
    Monitor, incr=7200., allowforce=False, verbose=False
))
#ENDINLINE

#INITVALUES
CFACTOR = P * Avogadro / R / TEMP * centi**3 * nano {ppb-to-molecules/cm3}
ALL_SPEC=1e-32*CFACTOR;
TOTALNOx=1.
M=1e9
O2=.21*M
N2=.79*M
H2O=0.01*M

""" + inittxt + """

#MONITOR t//3600; THETA; TEMP; O3;
#LOOKAT ALL;

#include cb6r5m_ae7_aq.eqn

""")

keywords = [
    'DUMMY', 'EMISSION', 'BNZHRXN', 'BNZNRXN', 'ISOPRXN', 'SESQRXN', 'SULRXN',
    'TOLHRXN', 'TOLNRXN', 'TRPRXN', 'XYLHRXN', 'XYLNRXN'
]

mech = Mech(
    infile, mechname='test', incr=360, add_default_funcs=False,
    keywords=keywords
)

unused_obs = [spc for spc in list(mdnobs_ppb.index) if spc not in mech.allspcs]
if len(unused_obs) > 0:
    print('**WARN:: unused obs')
    print('**', unused_obs)
    print('**WARN:: unused obs')

# Create time start/stop matrix
nhour = 24. * 10
start_end_ts = np.linspace(
    0, nhour, int(nhour) * 12 + 1
).repeat(2, 0)[1:-1].reshape(-1, 2) * 3600 + 12 * 3600


# Create specific tolerance values for radical species
# great care should be taken with O1D.
# note that CMAQ EBI solver has no absolute tolerance and
# uses 1 rtol for radicals and rxn counters.
#
fastspc = ['O', 'O1D', 'HCO3', 'NTR1']
fastspc += [s for s in mech.allspcs if s.endswith('RXN')]
atol = np.array([10 if spc in fastspc else 1e-3 for spc in mech.allspcs])
rtol = np.array([0.1 if spc in fastspc else 1e-5 for spc in mech.allspcs])

# Start timing the process
runstart = time.time()

# Update world for all interpolated values
mech.world['t'] = start_end_ts.min()
mech.update_y_from_world()
mech.Update_World(forceupdate=True)

# %%
# Run Model
# ---------
#

# Archive initial values
mech.archive()

# For each start/stop combination, update the world and integrate
for t0, t1 in start_end_ts:
    # set current time to t0
    t = t0
    # Get y vector from world state
    y = mech.update_y_from_world()
    # integrate from t to t1 with initial state of y
    ts, Y = mech.integrate(
        t, t1, y0=y, solver='odeint', mxords=3, mxordn=3, atol=atol,
        rtol=rtol, mxstep = 1000, hmax = 300, verbose=False
    )
    # ts, Y = mech.integrate(
    #     t, t1, y0=y, solver='lsoda', atol=atol, rtol=rtol, nsteps=1000,
    #     max_step=300, max_order_ns=3, max_order_s=3, verbose=False
    # )
    mixdt, = np.diff(ts)
    # 5% per day; converted to % per second and % per dt
    mech.world['y'] = Y[-1]
    if False:
        mixf = .05 / 24 / 3600 * mixdt
        Ybkg = np.copy(Y[-1])
        for s, b in mdnobs_ppb.bkg.items():
            if s in mech.allspcs:
                # convert ppb to #/cm3
                Ybkg[mech.allspcs.index(s)] = b * mech.world['CFACTOR']
        # Update world to match new time
        mixed = Y[-1] * (1 - mixf) + mixf * Ybkg
        mech.world['y'][:] = mixed

    mech.update_world_from_y()
    mech.Update_World(forceupdate=True)
    # Update y vector for mixing of PBL
    # Set new time to last time of integration
    t = ts[-1]
    # Update world from y vector
    mech.update_world_from_y()
    mech.archive()


# Optionally save to disk
# mech.output()

runend = time.time()
print((runend - runstart), 'seconds')

# %%
# Make Plots
# ----------

data = mech.get_output()
data['h_LST'] = data['t'] / 3600
noxspcs = ['NO', 'NO2', 'N2O5', 'HONO', 'NO3']
data['NOx'] = sum([data[ns] for ns in noxspcs if ns in data.columns])
data['TEMP/25 [C]'] = (data['TEMP'] - 273.15) / 25.
data['OH/10 [ppqv]'] = data['OH'] * 1e5
data['HO2 [pptv]'] = data['HO2'] * 1e3

fig, ax = plt.subplots()
data.set_index('h_LST')[['NOx'] + noxspcs].plot(ax=ax)
ax.set(xlabel='hour [LST]', ylabel='ppb', ylim=(.1, 5), yscale='log')
fig.savefig('henderson_nox.png')

fig, ax = plt.subplots()
data.set_index('h_LST')[['O3', 'OH/10 [ppqv]', 'HO2 [pptv]']].plot(ax=ax)
ax.set(xlabel='hour (LST)')
fig.savefig('henderson_ox.png')
