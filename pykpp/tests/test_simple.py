__doc__ = """
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
NO = 0.1 * TOTALNOx
NO2 = 0.9 * TOTALNOx
PHOx = .1e-3 * CFACTOR
{B = 210.; what to do about B?}
"""


def test_simple():
    from warnings import warn
    import numpy as np
    import io
    import pandas as pd
    from ..mech import Mech

    mechdefn = io.StringIO(mechstr) 
    mech1 = Mech(mechdefn, incr=3600, verbose=0)    
    mech1.run(verbose=1, debug=False, solver='odeint')
    ode_df = mech1.get_output()
    mechdefn = io.StringIO(mechstr) 
    mech2 = Mech(mechdefn, incr=3600, verbose=0)    
    mech2.run(
        verbose=1, debug=False, solver='lsoda', max_order_ns=2, max_order_s=2
    )
    lsoda_df = mech2.get_output()
    mechdefn = io.StringIO(mechstr) 
    mech3 = Mech(mechdefn, incr=3600, verbose=0)    
    mech3.run(verbose=1, debug=False, solver='vode')
    vode_df = mech3.get_output()
    for ok in ['t', 'O3']:
        sumdf = pd.DataFrame(dict(
            ode=ode_df[ok], lsoda=lsoda_df[ok], vode=vode_df[ok]
        ))
        sumdf['del'] = sumdf.max(axis=1) - sumdf.min(axis=1)
        sumdf['del%'] = sumdf['del'] / sumdf.mean(axis=1) * 100
        warn(f'\nKey: {ok}\n' + sumdf.iloc[::60].to_markdown())
    assert np.allclose(ode_df['O3'], lsoda_df['O3'], rtol=5)
    assert np.allclose(ode_df['O3'], vode_df['O3'], rtol=5)
