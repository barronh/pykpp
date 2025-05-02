pykpp
=====

pykpp is a KPP-like chemical mechanism parser that produces a box model solvable by SciPy's odeint solver


Run Examples No installation
----------------------------
* Requires a google account.
* Uses Google's Colab
* Click Badge [![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/barronh/pykpp/)


Running Instructions
--------------------

1. Open DOS Prompt on Windows or Terminal on Linux/Mac OSX
2. Type "python -m pykpp model.kpp" where model.kpp is a valid KPP input file.
3. Example below assumes Installation and that mcm_ch4.kpp has been downloaded to the working folder

```
$ python -m pykpp mcm_ch4.kpp --solver="odeint,rtol=1e-5,atol=[1e-5 if s != 'O1D' else 1 for s in out.allspcs]"
Found PYTHON!
Species: 30 
Reactions: 68
/Users/barronh/Development/pykpp/src/pykpp/updaters.py:7: UserWarning: Using RO2 = 1e-32 @ 18000.0
  def __init__(self, *args, **kwds):
 0.0%: {t:18000,O3:60,NO2:4.1E-43,RO2/CFACTOR:4.1E-43}
 1.0%: {t:18301,O3:59,NO2:0.98,RO2/CFACTOR:4.1E-43}
 2.1%: {t:18606,O3:59,NO2:0.98,RO2/CFACTOR:7.7E-05}
...
98.8%: {t:46464,O3:58,NO2:0.041,RO2/CFACTOR:0.046}
100.1%: {t:46837,O3:57,NO2:0.039,RO2/CFACTOR:0.047}
Solved in 0.410105 seconds
```

* Example input files can be downloaded from https://github.com/barronh/pykpp/tree/main/pykpp/models and are also in the models folder of pykpp source

Install Instructions
--------------------

Instructions
------------

1. Install Anaconda (https://www.continuum.io/downloads)
2. Download release version zip file from github. (or source)
3. Extract zip file
4. navigate to extracted folder pykpp-X.Y (the unzipped folder from pykpp-X.Y.zip)
5. Install pykpp from DOS terminal or terminal on Mac/Linux

    DOS: `python.exe setup.py install`
    Terminal: `python setup.py install`

6. Run `simple_organic` test case
  1. Copy `simple_organic.kpp` and `ekma.py` from test folder (either through github or from source
  2. navigate to where you downloaded `simple_organic.kpp` and `ekma.py`
7. Run pykpp

    DOS: `python.exe -m pykpp simple_organic.kpp -c ekma.py`
    Terminal: `python -m pykpp simple_organic.kpp -c ekma.py`

8. Compare ekma.png, which is an EKMA diagram, to `Figure 6.12 from Seinfeld & Pandis ACP 2nd ed. 2006` (search Google Books for that phrase to compare).
9. Do something more complicated (see Using pykpp with MCM)



Using pykpp with a MCM Extraction
---------------------------------

Adapt mcm kpp output for pykpp following these steps.

1. Fix line endings (mixed linux/windows)
2. Remove header comment (First 16 lines)
3. Remove everything above `#INCLUDE F90_RCONST`
4. replace `F90_RCONST` with `PY_UTIL`
5. replace `USE constants` with 2 lines

    `add_world_updater(func_updater(Update_M, incr = 300))`
    `add_world_updater(func_updater(Update_THETA, incr = 300))`
    
6. remove comment lines (preceded by `!`)
7. replace `RO2 = & C(ind_XXX1) + C(ind_XXX2) + ....` with `RO2 = y[[ind_XXX1, ind_XXX2, ...]].sum()`
8. remove `CALL mcm_constants(time, temp, M, N2, O2, RO2, H2O) ` line
9. replace all `D-` and `D+` wit `e-` and `e+` (also check for `#D#`, which is an implicit +, and replace as necessary)
10. wrap code in `add_world_updater(code_updater("""<CODE GOES HERE>""", incr = 300))`

For example, methane chemistry changes from 

```
{********************************************************************* ;
* A citation to the MCM website and the relevant mechanism          * ;
* construction protocols should be given in any publication using   * ;
* information obtained from this source, using the following or     * ;
* comparable wording:                                               * ;
* The chemical mechanistic information was taken from the Master    * ;
* Chemical Mechanism, MCM v3.3.1 (ref), via website:                  * ;
* http://mcm.leeds.ac.uk/MCM.                                       * ;
* The reference should be: (Jenkin et al., Atmos. Environ., 31, 81, * ;
* 1997; Saunders et al., Atmos. Chem. Phys., 3, 161, 2003), for     * ;
* non aromatic schemes; (Jenkin et al., Atmos. Chem. Phys., 3,  * ;
* 181, 2003; Bloss et al., Atmos. Chem. Phys., 5, 641, 2005), for   * ;
* aromatic schemes; (Jenkin et al., Atmos. Chem. Phys.,  12, * ;
* 5275, 2012), for the beta-caryophyllene scheme and (Jenkin et al., ;
* Atmos. Chem. Phys., 15, 11433, 2015), for the isoprene scheme.    * ;
********************************************************************* ;}
#INLINE F90_GLOBAL 
 REAL(dp)::M, N2, O2, RO2, H2O 
 #ENDINLINE {above lines go into MODULE KPP_ROOT_Global}
#INCLUDE atoms 
#DEFVAR
HCHO = IGNORE ;
...
#INLINE F90_RCONST
 USE constants
 !end of USE statements 
 !
 ! start of executable statements
 RO2 = & 
 C(ind_CH3O2) 
 KRO2NO = 2.7D-12*EXP(360/TEMP)
KRO2HO2 = 2.91D-13*EXP(1300/TEMP)
...
#ENDINLINE
```

and become

```
#INLINE PY_UTIL
add_world_updater(func_updater(Update_M, incr = 300))
add_world_updater(func_updater(Update_THETA, incr = 300))
add_world_updater(code_updater("""
EXP = exp
LOG10 = log10
try:
    RO2 = y[[ind_CH3O2]].sum()
except:
    warn('Using RO2 = 1e-32 @ %s' % t)
    RO2 = 1e-32

KRO2NO = 2.7e-12*EXP(360/TEMP)
KRO2HO2 = 2.91e-13*EXP(1300/TEMP)
KBPPN = (KPPN0*KPPNI)*FCPPN/(KPPN0+KPPNI)
""", incr = 300, message = 'MCM')
```

11. Replace `= :` with `= DUMMY :` in equations
12. Replace `J(##)` with `MCMJ(##, THETA)`
13. Move any species you want as constants from reactants to rate constant and comment them out as products.
14. Add intialization code. For example,

```
#INLINE PY_INIT
DT = 300
TSTART = 5 * 3600.
TEND = TSTART + 8 * 3600
P = 101325.
TEMP = 298.
H2O = h2o_from_rh_and_temp(80., TEMP)
StartDate = 'datetime(2013, 6, 1)'
Latitude_Degrees = 35.
Longitude_Degrees = 0.00E+00
#ENDINLINE

#INITVALUES
CFACTOR = P * Avogadro / R / TEMP * centi **3 * nano {ppb-to-molecules/cm3}
ALL_SPEC=1e-32;
CH4 = 1850
NO = 1.
NO2 = 7.
O3 = 60.
```

The mcm_ch4.kpp file under pykpp/models was extracted by selecting CH4 from the MCM website, and exporting KPP format with inorganic reactions and rate constants on 2016-01-28 at 10:00am and converted following the steps described above. It can be used as a reference.
