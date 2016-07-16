import time
import numpy as np
from pykpp.mech import Mech
from io import StringIO
import pandas as pd
import matplotlib.pyplot as plt
from pykpp.updaters import func_updater, Monitor

mech = Mech(StringIO("""
#INLINE PY_INIT
TEMP = 298.15
P = 99600.
t = TSTART = 0 * 3600.
TEND = TSTART + 3600. * 24
DT = 60.
MONITOR_DT = 3600
StartDate = 'datetime(2010, 7, 14)'
Latitude_Degrees = 40.
Longitude_Degrees = 0.00E+00
KDIL = 0.
PBLH_OLD = 250.
vd_CO = 0.03
vd_SO2 = 1.8
vd_NO = 0.016
vd_NO2 = 0.1
vd_NO3 = 0.1
vd_N2O5 = 4
vd_HONO = 2.2
vd_HNO3 = 4
vd_HNO4 = 4
vd_O3 = 0.4
vd_H2O2 = 0.5
vd_CH3OOH = 0.1
molps_to_molpcm2ps = Avogadro / 1200000**2

add_world_updater(func_updater(Update_M, incr = 360., verbose = False))
add_world_updater(func_updater(Update_THETA, incr = 360., verbose = False))
add_world_updater(interpolated_from_csv('summerenv.tsv', 'time', incr = 360., delimiter = '\\t', verbose = False))
add_world_updater(interpolated_from_csv('mean_emis.csv', 'time', incr = 360, verbose = False))
add_world_updater(func_updater(Monitor, incr = 7200., allowforce = False, verbose = False))
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
O3=30.
NO = 0.1 * TOTALNOx
NO2 = 0.9 * TOTALNOx
SO2 = 1
N2O = 320
{B = 210.; what to do about B?}

#MONITOR t/3600.; THETA; PBLH; TEMP; O3;
#LOOKAT THETA; t/3600; PBLH; TEMP; O1D; OH; HO2; O3; NO; NO2; eNOx; TUV_J(6, THETA); emis_NO; emis_NO * molps_to_molpcm2ps / PBLH / 100 / CFACTOR * 3600;

#include cb05cl.eqn
<EMISALD2> EMISSION = ALD2 : emis_ALD2 * molps_to_molpcm2ps / PBLH / 100 ;
<EMISALDX> EMISSION = ALDX : emis_ALDX * molps_to_molpcm2ps / PBLH / 100 ;
<EMISCH4> EMISSION = CH4 : emis_CH4 * molps_to_molpcm2ps / PBLH / 100 ;
<EMISCL2> EMISSION = CL2 : emis_CL2 * molps_to_molpcm2ps / PBLH / 100 ;
<EMISCO> EMISSION = CO : emis_CO * molps_to_molpcm2ps / PBLH / 100 ;
<EMISETH> EMISSION = ETH : emis_ETH * molps_to_molpcm2ps / PBLH / 100 ;
<EMISETHA> EMISSION = ETHA : emis_ETHA * molps_to_molpcm2ps / PBLH / 100 ;
<EMISETOH> EMISSION = ETOH : emis_ETOH * molps_to_molpcm2ps / PBLH / 100 ;
<EMISFORM> EMISSION = FORM : emis_FORM * molps_to_molpcm2ps / PBLH / 100 ;
<EMISHONO> EMISSION = HONO : emis_HONO * molps_to_molpcm2ps / PBLH / 100 ;
<EMISIOLE> EMISSION = IOLE : emis_IOLE * molps_to_molpcm2ps / PBLH / 100 ;
<EMISISOP> EMISSION = ISOP : emis_ISOP * molps_to_molpcm2ps / PBLH / 100 ;
<EMISMEOH> EMISSION = MEOH : emis_MEOH * molps_to_molpcm2ps / PBLH / 100 ;
<EMISNH3> EMISSION = NH3 : emis_NH3 * molps_to_molpcm2ps / PBLH / 100 ;
<EMISNO> EMISSION = NO + eNOx : emis_NO * molps_to_molpcm2ps / PBLH / 100 ;
<EMISNO2> EMISSION = NO2 + eNOx : emis_NO2 * molps_to_molpcm2ps / PBLH / 100 ;
<EMISOLE> EMISSION = OLE : emis_OLE * molps_to_molpcm2ps / PBLH / 100 ;
<EMISPAR> EMISSION = PAR : emis_PAR * molps_to_molpcm2ps / PBLH / 100 ;
<EMISSESQ> EMISSION = SESQ : emis_SESQ * molps_to_molpcm2ps / PBLH / 100 ;
<EMISSO2> EMISSION = SO2 : emis_SO2 * molps_to_molpcm2ps / PBLH / 100 ;
<EMISSULF> EMISSION = SULF : emis_SULF * molps_to_molpcm2ps / PBLH / 100 ;
<EMISTERP> EMISSION = TERP : emis_TERP * molps_to_molpcm2ps / PBLH / 100 ;
<EMISTOL> EMISSION = TOL : emis_TOL * molps_to_molpcm2ps / PBLH / 100 ;
<EMISXYL> EMISSION = XYL : emis_XYL * molps_to_molpcm2ps / PBLH / 100 ;
<DDEPCO> CO = DUMMY : vd_CO / PBLH / 100 ;
<DDEPSO2> SO2 = DUMMY : vd_SO2 / PBLH / 100 ;
<DDEPNO> NO = DUMMY : vd_NO / PBLH / 100 ;
<DDEPNO2> NO2 = DUMMY : vd_NO2 / PBLH / 100 ;
<DDEPNO3> NO3 = DUMMY : vd_NO3 / PBLH / 100 ;
<DDEPN2O5> N2O5 = DUMMY : vd_N2O5 / PBLH / 100 ;
<DDEPHONO> HONO = DUMMY : vd_HONO / PBLH / 100 ;
<DDEPHNO3> HNO3 = DUMMY : vd_HNO3 / PBLH / 100 ;
<DDEPHNO4> PNA = DUMMY : vd_HNO4 / PBLH / 100 ;
<DDEPO3> O3 = DUMMY : vd_O3 / PBLH / 100 ;
<DDEPH2O2> H2O2 = DUMMY : vd_H2O2 / PBLH / 100 ;
<DDEPCH3OOH> MEPX = DUMMY : vd_CH3OOH / PBLH / 100 ;

"""), mechname = 'summercb05', incr = 360, add_default_funcs = False, keywords = ['DUMMY', 'EMISSION', 'BNZHRXN', 'BNZNRXN', 'ISOPRXN', 'SESQRXN', 'SULRXN', 'TOLHRXN', 'TOLNRXN', 'TRPRXN', 'XYLHRXN', 'XYLNRXN'])

def UpdatePBL(mech, world):
    """
    Defining updater for PBL rise
    """
    if 'y' not in world: return
    PBLH_OLD = world['PBLH_OLD']
    ybkg = world['ybkg']
    y = world['y']
    PBLH = mech.world['PBLH']
    if PBLH != PBLH_OLD:
        KDIL = np.minimum(1, PBLH_OLD / PBLH)
        mech.world['KDIL'] = KDIL
        world['PBLH_OLD'] = PBLH
        ydil = (KDIL * y) + (1 - KDIL) * ybkg
        nonzero = ybkg != 2.4195878568940516e-22
        y[nonzero] = ydil[nonzero]
    else:
        mech.world['KDIL'] = 1
    mech.update_world_from_y(y)

# Create time start/stop matrix
nhour = 24.
start_end_ts = np.linspace(0, nhour, nhour*12+1).repeat(2, 0)[1:-1].reshape(-1, 2)*3600

# Capture initial values as "background"
# concentrations
mech.world['ybkg'] = mech.get_y(mech.parsed_world)

# Create specific tolerance values for radical species
# great care should be taken with O1D.
# note that CMAQ EBI solver has no absolute tolerance and
# uses 1 rtol for radicals and rxn counters.
#
atol = np.array([10 if (spc in ('O', 'O1D', 'HCO3', 'NTR') or spc[-3:] == 'RXN') else 1e-3 for spc in mech.allspcs])
rtol = np.array([0.1 if (spc in ('O', 'O1D', 'HCO3', 'NTR') or spc[-3:] == 'RXN') else 1e-5 for spc in mech.allspcs])

# Start timing the process
runstart = time.time()

# Update world for all interpolated values
mech.world['t'] = start_end_ts.min()
mech.update_y_from_world()
mech.Update_World(forceupdate = True)

# Add an updater that depends on the y vector
mech.add_world_updater(func_updater(UpdatePBL, incr = 360., verbose = False))
mech.world['PBLH_OLD'] = mech.world['PBLH']

# Archive initial values
mech.archive()
# NEI 2011 has Biogenics already
#def AddQuasiBio(mech, world):
#    # quasi Biogenic
#    world['emis_ISOP'] *= 15.
#
#mech.add_world_updater(func_updater(AddQuasiBio, incr = 360, verbose = False))

# For each start/stop combination, update the world and integrate
for t0, t1 in start_end_ts:
    # set current time to t0
    t = t0
    # Get y vector from world state
    y = mech.update_y_from_world()
    # integrate from t to t1 with initial state of y
    #ts, Y = mech.integrate(t, t1, y0 = y, solver = 'odeint', mxords = 3, mxordn = 3, atol = atol, rtol = rtol, mxstep = 1000, hmax = 300, verbose = False)
    ts, Y = mech.integrate(t, t1, y0 = y, solver = 'lsoda', atol = atol, rtol = rtol, nsteps = 1000, max_step = 300, max_order_ns = 3, max_order_s = 3, verbose = False)
    # Update world to match new time
    mech.world['y'] = Y[-1]
    mech.update_world_from_y()
    mech.Update_World(forceupdate = True)
    # Update y vector for mixing of PBL
    # Set new time to last time of integration
    t = ts[-1]
    # Update world from y vector
    mech.update_world_from_y()
    mech.archive()
    

mech.output()
runend = time.time()
print((runend - runstart), 'seconds')


mech._archive.seek(0, 0)
data = pd.read_csv(mech._archive, sep = '\t')
plt.figure()
plt.plot(data['t/3600'], data['emis_NO'], label = 'eNOx')
plt.plot(data['t/3600'], data['TUV_J(6,THETA)']*100, label = 'jNO2*100')
plt.plot(data['t/3600'], data['PBLH']/1000, label = 'PBLH/1000')
plt.plot(data['t/3600'], (data['TEMP'] - 273.15)/25., label = 'TEMPC/25')
plt.legend()
plt.savefig('physical.pdf')
plt.figure()
plt.plot(data['t/3600'], data['NO'] + data['NO2'], label = 'NOx')
plt.plot(data['t/3600'], data['NO'], label = 'NO')
plt.plot(data['t/3600'], data['NO2'], label = 'NO2')
plt.legend()
plt.ylim(0.1, 5)
plt.yscale('log')
plt.savefig('chemical_nox.pdf')
plt.figure()
plt.plot(data['t/3600'], data['O3'], label = 'O3 ppb')
plt.plot(data['t/3600'], data['OH']/10*1000**2, label = 'OH/10 ppqv')
plt.plot(data['t/3600'], data['HO2']*1000, label = 'HO2 pptv')
plt.legend()
plt.savefig('chemical_ozone.pdf')
