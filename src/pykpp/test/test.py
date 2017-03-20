from __future__ import print_function
__all__ = ['testit']
from matplotlib.mlab import csv2rec, rec2csv
import numpy as np
import os
import shutil
from StringIO import StringIO
from glob import glob
kpphome = os.environ.get('KPP_HOME', None)

_modelconfigs = dict(cbm4 = """
#DEFVAR
    NO          = IGNORE;       {nitric oxide}
    NO2         = IGNORE;      {nitrogen dioxide}
    NO3         = IGNORE;      {nitrogen trioxide}
    N2O5        = IGNORE;     {dinitrogen pentoxide}
    HONO        = IGNORE;  {nitrous acid}
    HNO3        = IGNORE;  { nitric acid }
    PNA         = IGNORE; {HO2NO2 peroxynitric acid}
    O1D         = IGNORE;           {oxygen atomic first singlet state}
    O           = IGNORE;           {oxygen atomic ground state (3P)}
    OH          = IGNORE;       {hydroxyl radical}
    O3          = IGNORE;          {ozone}
    HO2         = IGNORE;      {perhydroxyl radical}
    H2O2        = IGNORE;     {hydrogen peroxide}
    HCHO        = IGNORE;  {formalydehyde}
    ALD2        = IGNORE;      {high molecular weight aldehides}
    C2O3        = IGNORE; {CH3CO(O)OO peroxyacyl radical}
    PAN         = IGNORE; {CH3C(O)OONO2, peroxyacyl nitrate}
    PAR         = IGNORE;      {parafin carbon bond}
    ROR         = IGNORE;      {secondary organic oxy radical}
    OLE         = IGNORE;      {olefinic carbon bond}
    ETH         = IGNORE;     {CH2=CH2 ethene}
    TOL         = IGNORE;     {C6H5-CH3 toluene}
    CRES        = IGNORE;      {cresol and h.m.w. phenols}
    TO2         = IGNORE;      {toluene-hydroxyl radical adduct}
    CRO         = IGNORE;      {methylphenoxy radical}
    OPEN        = IGNORE;      {h.m.w. aromatic oxidation ring fragment}
    XYL         = IGNORE;    {C6H4-(CH3)2 xylene}
    MGLY        = IGNORE; {CH3C(O)C(O)H methylglyoxal}
    ISOP        = IGNORE;      {isoprene}
    XO2         = IGNORE;      {NO-to-NO2 operation}
    XO2N        = IGNORE;      {NO-to-nitrate operation}          
    CO          = IGNORE;       {carbon monoxide}

#DEFFIX
    H2O         = IGNORE;      {water}
    H2          = IGNORE;          {molecular hydrogen}
    O2          = IGNORE;          {molecular oxygen}
    N2          = IGNORE;          {molecular nitrogen}
    CH4         = IGNORE;      {methane}
    M           = IGNORE;      {third body}

#EQUATIONS {of the CBM-IV mechanism}

{ 1.} NO2 {+ hv} =  NO + O              :  8.89E-3*SUN ;
{ 2.} O {+ O2 + M} = O3               :  ARR2(1.4E+3, 1175.0) ;
{ 3.} O3 + NO = NO2                   :  ARR2(1.8E-12, -1370.0) ;
{ 4.} O + NO2 = NO                    :  9.3E-12 ;
{ 5.} O + NO2   = NO3                 :  ARR2(1.6E-13, 687.0) ;
{ 6.} O + NO   =  NO2                 :  ARR2(2.2E-13, 602.0) ;
{ 7.} O3 + NO2 =  NO3                 :  ARR2(1.2E-13, -2450.0) ;
{ 8.} O3 {+ hv} = O                     :  3.556E-04*SUN ; {4.0E-2*RCONST(1) ;}
{ 9.} O3 {+ hv} = O1D                   :  2.489E-05*SUN ; {2.8E-3*RCONST(1) ;}
{10.} O1D   = O                       :  ARR2(1.9E+8, 390.0)  ;

{11.} O1D + H2O = 2OH                 :  2.2E-10 ;
{12.} O3 + OH = HO2                   :  ARR2(1.6E-12, -940.0) ;
{13.} O3 + HO2 = OH                   :  ARR2(1.4E-14, -580.0) ;
{14.} NO3 {+ hv} = 0.89 NO2 + 0.89 O 
                 + 0.11 NO            :  1.378E-01*SUN ; {15.5*RCONST(1);}
{15.} NO3 + NO = 2 NO2                :  ARR2(1.3E-11, 250.0) ;
{16.} NO3 + NO2 = NO + NO2            :  ARR2(2.5E-14, -1230.0) ;
{17.} NO3 + NO2   =  N2O5             :  ARR2(5.3E-13, 256.0) ;
{18.} N2O5 + H2O = 2 HNO3             :  1.3E-21 ;
{19.} N2O5   =  NO3 + NO2             :  ARR2(3.5E+14, -10897.0) ;
{20.} 2 NO  =  2 NO2                  :  ARR2(1.8E-20, 530.0) ;


{21.} NO + NO2 + H2O = 2 HONO         :  4.4E-40 ;
{22.} OH + NO   =  HONO               :  ARR2(4.5E-13, 806.0) ;
{23.} HONO {+ hv} =  OH + NO            :  1.511e-03*SUN ; {0.17*RCONST(1);}
{24.} OH + HONO =  NO2                :  6.6E-12 ;
{25.} 2 HONO = NO + NO2               :  1.0E-20 ;
{26.} OH + NO2   =  HNO3              :  ARR2(1.0E-12, 713.0) ;
{27.} OH + HNO3   =  NO3              :  ARR2(5.1E-15, 1000.0) ;
{28.} HO2 + NO = OH + NO2             :  ARR2(3.7E-12, 240.0) ;
{29.} HO2 + NO2   =  PNA              :  ARR2(1.2E-13, 749.0) ;
{30.} PNA   = HO2 + NO2               :  ARR2(4.8E+13, -10121.0) ;

{31.} OH + PNA = NO2                  :  ARR2(1.3E-12, 380.0) ;
{32.} 2 HO2 = H2O2                    :  ARR2(5.9E-14, 1150.0)  ;
{33.} 2 HO2 + H2O = H2O2              :  ARR2(2.2E-38, 5800.0) ;
{34.} H2O2 {+ hv} = 2 OH                :  6.312E-06*SUN ; {7.1E-4*RCONST(1);}
{35.} OH + H2O2 = HO2                 :  ARR2(3.1E-12, -187.0) ;
{36.} OH + CO  = HO2                  :  2.2E-13 ;
{37.} HCHO + OH  =  HO2 + CO          :  1.0E-11 ;
{38.} HCHO {+ hv + 2 O2} = 2 HO2 + CO :  2.845E-05*SUN ; {3.2E-3*RCONST(1);}
{39.} HCHO +  hv = CO                 :  3.734E-05*SUN ; {4.2E-3*RCONST(1);}
{40.} HCHO + O = OH + HO2 + CO        :  ARR2(3.0E-11, -1550.0) ;

{41.} HCHO + NO3  = HNO3 
                        + HO2 + CO    :  6.3E-16 ;
{42.} ALD2 + O  =  C2O3 + OH          :  ARR2(1.2E-11, -986.0) ;
{43.} ALD2 + OH = C2O3                :  ARR2(7.0E-12, 250.0) ;
{44.} ALD2 + NO3  = C2O3 + HNO3       :  2.5E-15  ;
{45.} ALD2 {+ hv + 2 O2} = HCHO + XO2 
                         + CO + 2 HO2 :  4.00E-06*SUN ; {4.5E-4*RCONST(1);}
{46.} C2O3 + NO  = HCHO + XO2 
                       + HO2 + NO2    :  ARR2(5.4E-12, 250.0) ;
{47.} C2O3 + NO2 = PAN                :  ARR2(8.0E-20, 5500.0) ;
{48.} PAN = C2O3 + NO2                :  ARR2(9.4E+16, -14000.0) ;
{49.} 2 C2O3 = 2 HCHO + 2 XO2 + 2 HO2 :  2.0E-12  ;
{50.} C2O3 + HO2 = 0.79 HCHO 
   + 0.79 XO2 + 0.79 HO2 + 0.79 OH    :  6.5E-12 ;

{51.} OH = HCHO + XO2 + HO2           :  ARR2(1.1E+2, -1710.0) ;
{52.} PAR + OH = 0.87 XO2 + 0.13 XO2N 
               + 0.11 HO2 + 0.11 ALD2
               + 0.76 ROR - 0.11 PAR  :  8.1E-13 ;
{53.} ROR = 1.1 ALD2 + 0.96 XO2 
               + 0.94 HO2 + 0.04 XO2N
               + 0.02 ROR - 2.10 PAR  :  ARR2(1.0E+15, -8000.0) ;
{54.} ROR = HO2                       :  1.6E+03 ;
{55.} ROR + NO2 =  PROD               :  1.5E-11  ;
{56.} O + OLE = 0.63 ALD2 + 0.38 HO2 
                + 0.28 XO2 + 0.3 CO
                + 0.2 HCHO + 0.02 XO2N 
                + 0.22 PAR + 0.2 OH   :  ARR2(1.2E-11, -324.0) ;
{57.} OH + OLE = HCHO + ALD2 + XO2 
                 + HO2 - PAR          :  ARR2(5.2E-12, 504.0) ;
{58.} O3 + OLE = 0.5 ALD2 + 0.74 HCHO 
                 + 0.33 CO + 0.44 HO2 
                 + 0.22 XO2
                 + 0.1 OH - PAR       :  ARR2(1.4E-14, -2105.0)  ;
{59.} NO3 + OLE = 0.91 XO2 + HCHO 
                  + ALD2 + 0.09 XO2N                          
                  + NO2 - PAR         :  7.7E-15 ;
{60.} O + ETH = HCHO + 0.7 XO2 + CO
                + 1.7 HO2 + 0.3 OH    :  ARR2(1.0E-11, -792.0) ;

{61.} OH + ETH = XO2 + 1.56 HCHO + HO2 + 0.22 ALD2        :  ARR2(2.0E-12, 411.0)  ;
{62.} O3 + ETH = HCHO + 0.42 CO + 0.12 HO2                :  ARR2(1.3E-14, -2633.0) ;
{63.} OH + TOL = 0.08 XO2 + 0.36 CRES 
                 + 0.44 HO2 + 0.56 TO2                    :  ARR2(2.1E-12, 322.0) ;
{64.} TO2 + NO =  0.9 NO2 + 0.9 OPEN + 0.9 HO2            :  8.1E-12 ;
{65.} TO2 = HO2 + CRES                                    :  4.20 ;
{66.} OH + CRES = 0.4 CRO + 0.6 XO2 + 0.6 HO2 + 0.3 OPEN  :  4.1E-11  ;
{67.} NO3 + CRES = CRO + HNO3                             :  2.2E-11 ;
{68.} CRO + NO2 = PROD                                    :  1.4E-11 ;
{69.} OH + XYL = 0.7 HO2 + 0.5 XO2 + 0.2 CRES + 0.8 MGLY
                 + 1.10 PAR + 0.3 TO2                     :  ARR2(1.7E-11, 116.0) ;
{70.} OH + OPEN = XO2 + C2O3 + 2 HO2 + 2 CO + HCHO        :  3.0E-11 ;

{71.} OPEN {+ hv} = C2O3 + CO + HO2                         :  5.334E-05*SUN ; {6.0E-3*RCONST(1);}
{72.} O3 + OPEN = 0.03 ALD2 + 0.62 C2O3 
                  + 0.7 HCHO + 0.03 XO2 + 0.69 CO 
                  + 0.08 OH + 0.76 HO2 + 0.2 MGLY         :  ARR2(5.4E-17, -500.0)  ;
{73.} OH + MGLY =  XO2 + C2O3                             :  1.70E-11 ;
{74.} MGLY {+ hv} = C2O3 + CO + HO2                         :  1.654E-04*SUN ; {1.86E-2*RCONST(1);}
{75.} O + ISOP =  0.6 HO2 + 0.8 ALD2 + 0.55 OLE + 0.5 XO2
                  + 0.5 CO + 0.45 ETH + 0.9 PAR           :  1.80E-11  ;
{76.} OH + ISOP = HCHO + XO2 + 0.67 HO2 
                  + 0.4 MGLY + 0.2 C2O3
                  + ETH + 0.2 ALD2 + 0.13 XO2N            :  9.6E-11  ;
{77.} O3 + ISOP = HCHO + 0.4 ALD2 + 0.55 ETH + 0.2 MGLY 
                  + 0.06 CO + 0.1 PAR + 0.44 HO2 + 0.1 OH :  1.2E-17 ;
{78.} NO3 + ISOP =  XO2N                                  :  3.2E-13 ;
{79.} XO2 + NO = NO2                                      :  8.1E-12 ;
{80.} 2 XO2 =  PROD                                       :  ARR2(1.7E-14, 1300.0) ;
{81.} XO2N + NO =  PROD                                   :  6.8E-13 ;
 
#HESSIAN OFF
#INTEGRATOR kpp_lsode
#LANGUAGE   Fortran90 
#DRIVER     general 

#LOOKATALL

#MONITOR O3; NO; NO2; HO2;

#INITVALUES
  CFACTOR = 2.55E+10; {ppb-to-mcm}
  ALL_SPEC = 1.0E-8;
{Variable species}
  NO   = 50.0;
  NO2  = 20.0;
  HONO = 1.0;
  O3   = 100.0;
  HCHO = 10.0;
  ALD2 = 10;
  PAN  = 1.0;
  PAR  = 50.0;
  OLE  = 10.0;
  ETH  = 10.0;
  TOL  = 10.0;
  XYL  = 10.0;
  ISOP = 10.0;
  CO   = 300.0;
{Fixed species}
  H2O  = 1.25E+8; {30 %}


#INLINE F90_INIT
        TSTART = 12.D0*3600.D0
        TEND = TSTART + 24.D0*3600.D0 * 5
        DT = 3600.D0
        TEMP = 288.15
#ENDINLINE

#INLINE F90_GLOBAL
        REAL(dp) M
#ENDINLINE
""",
small_strato = """
#DEFVAR
O   = IGNORE;
O1D = IGNORE;
O3  = IGNORE;
NO  = IGNORE;
NO2 = IGNORE;    { Nitrogen dioxide }


#DEFFIX
M   = IGNORE;{ Atmospheric generic molecule }
O2  = IGNORE;        { Molecular oxygen }

#EQUATIONS { Small Stratospheric Mechanism }


<R1>  O2   + hv = 2O		: (2.643E-10) * SUN*SUN*SUN;
<R2>  O    + O2 = O3		: (8.018E-17);
<R3>  O3   + hv = O   + O2 	: (6.120E-04) * SUN;
<R4>  O    + O3 = 2O2		: (1.576E-15);
<R5>  O3   + hv = O1D + O2	: (1.070E-03) * SUN*SUN;
<R6>  O1D  + M  = O   + M	: (7.110E-11);
<R7>  O1D  + O3 = 2O2 		: (1.200E-10);
<R8>  NO   + O3 = NO2 + O2 	: (6.062E-15);
<R9>  NO2  + O  = NO  + O2	: (1.069E-11);
<R10> NO2  + hv = NO  + O	: (1.289E-02) * SUN;
 
#LOOKATALL
#MONITOR O3;O2;O;NO;O1D;NO2; {Screen Output}

#INITVALUES                    {Initial Values}
 
CFACTOR = 1.    ;              {Conversion Factor} 
O1D = 9.906E+01 ; 
O   = 6.624E+08 ; 
O3  = 5.326E+11 ; 
O2  = 1.697E+16 ;
NO  = 8.725E+08 ; 
NO2 = 2.240E+08 ; 
M   = 8.120E+16 ;
 
#INLINE F90_INIT
        TSTART = (12*3600)
        TEND = TSTART + (3*24*3600)
        DT = 0.25*3600  
        TEMP = 270
#ENDINLINE         
 
#INLINE F90_GLOBAL
        REAL(dp) M
#ENDINLINE         

#HESSIAN OFF
#INTEGRATOR kpp_lsode
#LANGUAGE   Fortran90
#DRIVER     general
""")
saprc99 = """
#MODEL      saprc99
#INTEGRATOR kpp_lsode
#LANGUAGE   Fortran90 
#DRIVER     general 
#HESSIAN OFF
#include saprc99.spc
#include saprc99.eqn

#LOOKATALL

#MONITOR O3; NO; NO2; ETHENE;

#INITVALUES

   CFACTOR = 2.4476e+13;

   ALL_SPEC = 0.0e0;
   NO = 1.0e-1;
   NO2 = 5.0e-2;
   HONO = 1.e-3;
   SO2 = 5.e-2;
   HCHO = 1.121e-2;
   CCHO = 2.316e-3;
   RCHO = 1.72e-3;
   ACET = 5.07e-3;
   MEK = 3.26e-3;
   MEOH = 5.89e-3;
   GLY = 1.21e-4;
   MGLY = 8.37e-5;
   PHEN = 6.06e-4;
   CRES = 5.60e-4;
   BALD = 7.51e-5; 
   METHACRO = 1.30e-3; 
   ISOPROD = 8.93e-5; 
   PROD2 = 1.93e-3; 
   ETHENE = 1.89e-2;                 
   ISOPRENE = 4.33e-4;
   ALK1 = 1.167e-2;
   ALK2 = 1.88e-2;
   ALK3 = 4.69e-2;
   ALK4 = 4.17e-2;
   ALK5 = 3.06e-2;   
   ARO1 = 1.18e-2;
   ARO2 = 8.74e-3;
   OLE1 = 1.04e-2;
   OLE2 = 7.97e-3;   
   TERP = 8.20e-4;
   XC   = 0.2E0;
   CCO_OH = 1.16e-3; 
   RCO_OH = 3.92e-4;     
   HCOOH  = 6.77e-4;
   O3P    = 7.843e-9;
   H2O = 2.0e+04; 
   O2  = 2.09e+5;
   AIR = 1.0e+6;
   CH4 = 1.0e0;
     
#INLINE F77_INIT
        TSTART = 12.0D0*3600.0D0
        TEND   = TSTART + 120.0D0*3600.0D0
        DT     = 3600.D0
        TEMP   = 300.0D0
#ENDINLINE   
   
#INLINE F90_INIT
        M = CFACTOR * 1e6
        TSTART = 12.0d0*3600.0d0
        TEND   = TSTART + 120.0d0*3600.0d0
        DT     = 3600.d0
        TEMP   = 300.0d0
#ENDINLINE   

#INLINE F90_GLOBAL
        REAL*8 M
#ENDINLINE F90_GLOBAL

#INLINE MATLAB_INIT
   global TSTART TEND DT TEMP
   TSTART = 12*3600;
   TEND   = TSTART + 120*3600;
   DT     = 3600;
   TEMP   = 300;
#ENDINLINE   
     
#INLINE C_INIT
   TSTART = 12.0*3600.0;
   TEND   = TSTART + 120.0*3600.0;
   DT     = 3600.0;
   TEMP   = 300.0;
#ENDINLINE  
"""
_allmodels = _modelconfigs.keys()

def runmodels(models = _allmodels, pykpp = True, kpp = True, verbose = False):
    for model in models:
        if os.path.isfile(model):
            modeldef = model
            model = os.path.splitext(os.path.basename(modeldef))[0]
            os.system('rm -rf %(model)s' % locals())
            os.mkdir(model)
            os.system('cp ' + ' '.join(glob(os.path.splitext(modeldef)[0] + '.*')) + ' ./%s' % model)
        else:
            modeldef = model + '.kpp'
            os.system('rm -rf %(model)s' % locals())
            os.mkdir(model)
            if model in _modelconfigs:
                print('Test uses its own defintion of %s' % model)
                file(os.path.join(model, modeldef), 'w').write(_modelconfigs[model])
            elif kpp:
                path = os.path.join(kpphome, 'examples', model + '_f90.kpp')
                if not os.path.exists(path):
                    path = os.path.join(kpphome, 'models', model + '.def')
                if not os.path.exists(path):
                    raise IOError('No kpp or def file could be found for %(model)s' % locals())
                shutil.copy(path, os.path.join(model, modeldef))
            else:
                raise IOError('Can only run cbm4, small_strato, and saprc99 tests without KPP installed and KPP_HOME defined')

        if kpp:
            exit_code = os.system('cd %(model)s && kpp %(modeldef)s && make -f Makefile_%(model)s && ./%(model)s.exe' % locals())
            if exit_code != 0:
                raise Exception('KPP failed; see above')
        if pykpp: os.system('cd %(model)s && python -m pykpp -j %(modeldef)s' % locals())


def makediffs(models = _allmodels, verbose = False, kpp = True):
    for model in models:
        model = os.path.splitext(os.path.basename(model))[0]
        if kpp:
            kppdat = csv2rec(os.path.join(model, model + '.dat'), delimiter = ' ')
        else:
            if model not in _modelconfigs:
                raise IOError('If KPP is not properly installed, you cannot run tests on mechanisms other than cbm4, saprc99, and small_strato.')
            kppdat = csv2rec(os.path.join(os.path.dirname(__file__), model + '.dat'), delimiter = ' ')
        pykppdat = csv2rec(os.path.join(model, model + '.pykpp.dat'), delimiter = ',')
        diff = pykppdat.copy()
        pct = pykppdat.copy()
        keys = set(kppdat.dtype.names).intersection(pykppdat.dtype.names)
        notkeys = set(pykppdat.dtype.names).difference(kppdat.dtype.names)
        notkeys.remove('t')
        for k in notkeys:
            diff[k] = np.nan
            pct[k] = np.nan
    
        for k in keys:
            diff[k] = pykppdat[k] - kppdat[k][:]
            pct[k] = diff[k] / kppdat[k][:] * 100
        diff['t'] = pykppdat['t'] - (kppdat['time'] * 3600. + pykppdat['t'][0])
        pct['t'] = diff['t'] / (kppdat['time'] * 3600. + pykppdat['t'][0]) * 100
        
        rec2csv(diff, os.path.join(model, model + '.diff.csv'), delimiter = ',')
        rec2csv(pct, os.path.join(model, model + '.pct.csv'), delimiter = ',')

def checkmodels(models = _allmodels, verbose = False, kpp = True):
    for model in models:
        model = os.path.splitext(os.path.basename(model))[0]
        if kpp:
            kppdat = csv2rec(os.path.join(model, model + '.dat'), delimiter = ' ')
        else:
            kppdat = csv2rec(os.path.join(os.path.dirname(__file__), model + '.dat'), delimiter = ' ')
        pykppdat = csv2rec(os.path.join(model, model + '.pykpp.dat'), delimiter = ',')
        diff = csv2rec(os.path.join(model, model + '.diff.csv'), delimiter = ',')
        pct = csv2rec(os.path.join(model, model + '.pct.csv'), delimiter = ',')
        keys = set(kppdat.dtype.names).intersection(pykppdat.dtype.names)
        keys.add('time')
        output = StringIO('')
        outputs = {}
        print('%20s,%s,%s,%s,%s,%10s,%10s,%10s,%10s,%10s' % ('Key', 'Type  ', 'Time pct', 'Time    ', 'Idx', 'KPP', 'PYKPP', 'Absolute', 'Percent', 'Check'), file =   output)
        trans = dict(time = 't')
        check = True
        checks = {False: [], True: []}
        for k in keys:
            thisdiff = diff[trans.get(k, k)]
            thispct = pct[trans.get(k, k)]
            thispykpp = pykppdat[trans.get(k, k)]
            thiskpp = kppdat[k]
            if k == 'time':
                thiskpp *= 3600.
                thiskpp += thispykpp[0]
            thismaxdiff = thispykpp.max() - thiskpp.max()
            thismaxpct = thismaxdiff / thiskpp.max() * 100
            absidx = np.abs(thisdiff).argmax()
            pctidx = np.abs(thispct).argmax()
            abspass = ('PASS' if np.abs(thispct[absidx]) < 2 else 'FAIL') if np.abs(thisdiff[absidx]) > .001 else 'PASS'
            pctpass = ('PASS' if np.abs(thispct[pctidx]) < 2 else 'FAIL') if np.abs(thisdiff[pctidx]) > .001 else 'PASS'
            endpass = ('PASS' if np.abs(thispct[-1]) < 2 else 'FAIL') if np.abs(thisdiff[-1]) > .001 else 'PASS'
            maxpass = ('PASS' if np.abs(thismaxpct) < 2 else 'FAIL') if np.abs(thismaxdiff) > .001 else 'PASS'
            thischeck = all([abspass == 'PASS', pctpass == 'PASS', endpass == 'PASS', maxpass == 'PASS'])
            check = check and thischeck
            checks[thischeck].append(k)
            outputs.setdefault(k, []).append('%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%,%s' % (k.upper(), 'MaxAbs', 100 * absidx / len(pykppdat['t']), pykppdat['t'][absidx], absidx, thiskpp[absidx], thispykpp[absidx], thisdiff[absidx], thispct[absidx], abspass))
            outputs[k].append('%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%,%s' % (k.upper(), 'MaxPct', 100 * pctidx / len(pykppdat['t']), pykppdat['t'][pctidx], absidx, thiskpp[pctidx], thispykpp[pctidx], thisdiff[pctidx], thispct[pctidx], pctpass))
            outputs[k].append('%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%,%s' % (k.upper(), 'Max   ', -1, -1, -1, thiskpp.max(), thispykpp.max(), thismaxdiff, thismaxpct, maxpass))
            outputs[k].append('%20s,%s,%7.1f%%,%8.0f,%3d,%10.3g,%10.3g,%10.3g,%8.2f%%,%s' % (k.upper(), 'Ending', 100, pykppdat['t'][-1], -1, thiskpp[-1], thispykpp[-1], thisdiff[-1], thispct[-1], endpass))
        print(keys)
        if not check or verbose:
            fails = checks[False]
            passes = checks[True]
            fails.sort()
            passes.sort()
            for k in passes + fails:
                for l in outputs[k]:
                    print(l, file = output)
            
            output.seek(0, 0)
            print(output.read())
            print("""
**Notes:
    1)  If using a non-rosenbrock solver, you must modify the 
        solver to update SUN and RCONST during integration
    2)  Mechanism using double precision numbers
        in real precision calculations will fails
        a) saprc99 - H2O2 is a known example of this
        b) to fix this problem with KPP, duplicate rate functions
           with double precision inputs (EP3 now has a DP3 analog).
           Change rates who raise underflow warnings during compilation
""")
        print(model, 'PASS' if check else 'FAIL')

def testit(*models, **kwds):
    verbose = kwds.pop('verbose', False)
    if len(models) == 0:
        models = _allmodels
    runmodels(models = models, verbose = verbose, kpp = not kpphome is None)
    makediffs(models = models, verbose = verbose, kpp = not kpphome is None)
    print('\n\n\n\n')
    checkmodels(models = models, verbose = verbose, kpp = not kpphome is None)

if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.set_usage("""Usage: python -m pykpp.test [-v]""")

    parser.add_option("-v", "--verbose", dest = "verbose", action = "store_true", default = False, help = "Show extended output")
    options, args = parser.parse_args()
    testit(verbose = options.verbose)