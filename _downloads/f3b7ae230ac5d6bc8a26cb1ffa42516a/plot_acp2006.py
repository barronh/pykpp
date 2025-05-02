"""
Simple Chemistry Example
========================

This test is a simple example from Seinfeld and Pandis ACP Second edition
on pg 240.[1].

    we extract a basic chemical mechanism with 1 primary organic oxidation
    (R1). Under high NOx conditions, the primary organic oxidation produces
    2 peroxy radical NO oxidations (R2,R3) to produce more radicals (i.e.,
    HO2 or OH). The produced NO2 can photolyze (R4) to reproduce NO and an
    odd oxygen (O3P) that is assumed to instantaneously produce ozone. Under
    high-NOx conditions, the ozone can be lost by oxidizing NO (R5). Finally,
    radicals can be removed lost by producing nitric acid (R6) or peroxides
    (R7,R8). Lastly, we add an artificial source of HOx (here defined as HO2
    + OH) (R9).

[Seinfeld and Pandis Pg 240](https://books.google.com/books?id=YH2K9eWsZOcC&lpg=PT220&vq=RH%20PHOx&pg=PT220#v=onepage&f=false) 
"""

# %%
# Imports
# -------
# 

import io
import matplotlib.pyplot as plt
from pykpp.mech import Mech
import numpy as np

# %%
# Define Model
# ------------
# A model is made up of:
#
# - equations to be solved
# - inline code to be run before equations are solved (PY_INIT),
# - options for saving data, and
# - state initialization (INITVALUES).
# 

mechstr = """
#EQUATIONS
{R1} RH + OH = RO2 + H2O : 26.3e-12;
{R2} RO2 + NO = NO2 + RCHO + HO2 : 7.7e-12;
{R3} HO2 + NO = NO2 + OH : 8.1e-12;
{R4} NO2 = NO + O3 : jno2 ;
{R5} O3 + NO = NO2 + O2 : 1.9e-14 ;
{R6} OH + NO2 = HNO3 : kohno2 ; 
{R7} HO2 + HO2 = H2O2 + O2 : 2.9e-12 ;
{R8} RO2 + HO2 = ROOH + O2 : 5.2e-12 ;
{R9} EMISSION = OH : PHOx;

#INLINE PY_INIT
t=TSTART=6*3600
TEND=10*3600+TSTART
P = 99600.
TEMP = 298.15
DT = 60.
MONITOR_DT = 3600.
StartDate = 'datetime(2010, 7, 14)'
Latitude_Degrees = 40.
Longitude_Degrees = 0.00E+00

jno2 = .015;
kohno2 = 1.1e-11;
#ENDINLINE

#MONITOR O3; RH; NO; NO2; OH; HO2;
#INTEGRATOR odeint;

#INITVALUES
CFACTOR = P * Avogadro / R / TEMP * centi **3 * nano {ppb-to-molecules/cm3}
ALL_SPEC=1e-32*CFACTOR;
M=1e9
TOTALNOx=10.
RH = 200.
O2=.21*M
N2=.79*M
H2O=0.01*M
O3=30.
NO = TOTALNOx * 2 / 3
NO2 = TOTALNOx * 1 / 3
PHOx = .1e-3 * CFACTOR
{B = 210.; what to do about B?}
"""

# %%
# Run the model and visualize
# ---------------------------
# - load the model in a mechanism,
# - run it using the odeint solver,
# - get the output as a pandas dataframe,
# - add some derived variables for plotting, and
# - make one plot for NOx species, and one for ozone, OH, and HO2

mech = Mech(io.StringIO(mechstr), incr=3600, verbose=0)    
runtime = mech.run(verbose=1, debug=False, solver='odeint')
df = mech.get_output()

df['h_LST'] = df.eval('t / 3600')
df['NOx'] = df.eval('NO + NO2')
df['OH/10 ppqv'] = df.eval('OH * 1e5')
df['HO2 pptv'] = df.eval('HO2 * 1e3')

fig, ax = plt.subplots()
df.set_index('h_LST')[['NO', 'NO2', 'NOx']].plot(ax=ax)
ax.set(yscale='log')
fig.savefig('simple_nox.png')

fig, ax = plt.subplots()
df.set_index('h_LST')[['O3', 'OH/10 ppqv', 'HO2 pptv']].plot(ax=ax)
fig.savefig('simple_ox.png')

# %%
# Create Ozone Isopleths like Figure 6.12
# ---------------------------------------
# - reload the model in a mechanism and disable monitoring,
# - define the NOx and VOC increments to run to find ozone points,
# - run each combniation of NOx and VOC it using the vode solver,
# - save ozone as a 2d array,
# - make a plot of the contours of ozone from all NOx and VOC combos.

mech = Mech(io.StringIO(mechstr), monitor_incr=None)
noxppbs = [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 20, 30 , 40, 50]
vocppbs = [50, 60, 70, 80, 90, 100, 200, 300, 400, 500]

o3vals = np.zeros((len(noxppbs), len(vocppbs)), dtype = 'f')

for ni, n in enumerate(noxppbs):
    print(f'NOx [ppb]: {n}', end=': VOC [ppb]:')
    for vi, v in enumerate(vocppbs):
        print(v, end='.', flush=True)
        # reset world to ensure no carryover from previous simulation
        mech.resetworld()
        CFACTOR = mech.world['CFACTOR']
        # override the default initialization of NO, NO2 and RH (VOC)
        mech.world['NO'] = n * CFACTOR * 2 / 3
        mech.world['NO2'] = n * CFACTOR * 1 / 3
        mech.world['RH'] = v * CFACTOR
        # Run for this combination of NOx and VOC and store result
        mech.run(solver = 'vode')
        o3vals[ni, vi] = eval('O3 / CFACTOR', None, mech.world)
    print()

fig, ax = plt.subplots()
levels = np.linspace(30, 360, 12)
cs = ax.contourf(vocppbs, noxppbs, o3vals, levels=levels)
fig.colorbar(cs, label='ozone ppb')
ax.set_yscale('log', base=10, subs=[2,3,5])
ax.set_xscale('log', base=10, subs=[2,3,5])
ax.set(title='Ozone Isoplets (PHOx=0.1 ppt/s)', xlabel='RH [ppb]', ylabel='NOx [ppb]')
fig.savefig('ekma.png')

