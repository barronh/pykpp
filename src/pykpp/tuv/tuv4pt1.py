import numpy as np
tuv4pt1_data = {100: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   100.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          0.000E+00
   2 O3 -> O2 + O(1D)                          0.000E+00
   3 O3 -> O2 + O(3P)                          0.000E+00
   4 NO2 -> NO + O(3P)                         0.000E+00
   5 NO3 -> NO + O2                            0.000E+00
   6 NO3 -> NO2 + O(3P)                        0.000E+00
   7 N2O5 -> NO3 + NO + O(3P)                  0.000E+00
   8 N2O5 -> NO3 + NO2                         0.000E+00
   9 N2O + hv -> N2 + O(1D)                    0.000E+00
  10 HO2 + hv -> OH + O                        0.000E+00
  11 H2O2 -> 2 OH                              0.000E+00
  12 HNO2 -> OH + NO                           0.000E+00
  13 HNO3 -> OH + NO2                          0.000E+00
  14 HNO4 -> HO2 + NO2                         0.000E+00
  15 CH2O -> H + HCO                           0.000E+00
  16 CH2O -> H2 + CO                           0.000E+00
  17 CH3CHO -> CH3 + HCO                       0.000E+00
  18 CH3CHO -> CH4 + CO                        0.000E+00
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     0.000E+00
  21 CHOCHO -> products                        0.000E+00
  22 CH3COCHO -> products                      0.000E+00
  23 CH3COCH3                                  0.000E+00
  24 CH3OOH -> CH3O + OH                       0.000E+00
  25 CH3ONO2 -> CH3O+NO2                       0.000E+00
  26 PAN + hv -> products                      0.000E+00
  27 ClOO + hv -> Products                     0.000E+00
  28 ClONO2 + hv -> Cl + NO3                   0.000E+00
  29 ClONO2 + hv -> ClO + NO2                  0.000E+00
  30 CH3Cl + hv -> Products                    0.000E+00
  31 CCl2O + hv -> Products                    0.000E+00
  32 CCl4 + hv -> Products                     0.000E+00
  33 CClFO + hv -> Products                    0.000E+00
  34 CF2O + hv -> Products                     0.000E+00
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     0.000E+00
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     0.000E+00
  37 CF3CF2Cl (CFC-115) + hv -> Products       0.000E+00
  38 CCl3F (CFC-11) + hv -> Products           0.000E+00
  39 CCl2F2 (CFC-12) + hv -> Products          0.000E+00
  40 CH3CCl3 + hv -> Products                  0.000E+00
  41 CF3CHCl2 (HCFC-123) + hv -> Products      0.000E+00
  42 CF3CHFCl (HCFC-124) + hv -> Products      0.000E+00
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     0.000E+00
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     0.000E+00
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  0.000E+00
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  0.000E+00
  47 CHClF2 (HCFC-22) + hv -> Products         0.000E+00
  48 BrONO2 + hv -> Br + NO3                   0.000E+00
  49 BrONO2 + hv -> BrO + NO2                  0.000E+00
  50 CH3Br + hv -> Products                    0.000E+00
  51 CHBr3                                     0.000E+00
  52 CF3Br (Halon-1301) + hv -> Products       0.000E+00
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  0.000E+00
  54 CF2Br2 (Halon-1202) + hv -> Products      0.000E+00
  55 CF2BrCl (Halon-1211) + hv -> Products     0.000E+00
  56 Cl2 + hv -> Cl + Cl                       0.000E+00""",
90: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   90.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          5.716E-40
   2 O3 -> O2 + O(1D)                          3.121E-08
   3 O3 -> O2 + O(3P)                          5.786E-06
   4 NO2 -> NO + O(3P)                         1.433E-04
   5 NO3 -> NO + O2                            2.817E-04
   6 NO3 -> NO2 + O(3P)                        2.447E-03
   7 N2O5 -> NO3 + NO + O(3P)                  5.171E-14
   8 N2O5 -> NO3 + NO2                         3.757E-07
   9 N2O + hv -> N2 + O(1D)                    3.569E-37
  10 HO2 + hv -> OH + O                        1.782E-34
  11 H2O2 -> 2 OH                              5.406E-08
  12 HNO2 -> OH + NO                           3.074E-05
  13 HNO3 -> OH + NO2                          2.124E-09
  14 HNO4 -> HO2 + NO2                         1.342E-08
  15 CH2O -> H + HCO                           1.462E-07
  16 CH2O -> H2 + CO                           4.612E-07
  17 CH3CHO -> CH3 + HCO                       1.050E-08
  18 CH3CHO -> CH4 + CO                        2.410E-16
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     5.911E-08
  21 CHOCHO -> products                        1.209E-06
  22 CH3COCHO -> products                      1.708E-06
  23 CH3COCH3                                  1.269E-09
  24 CH3OOH -> CH3O + OH                       4.583E-08
  25 CH3ONO2 -> CH3O+NO2                       2.819E-09
  26 PAN + hv -> products                      3.987E-09
  27 ClOO + hv -> Products                     4.431E-33
  28 ClONO2 + hv -> Cl + NO3                   5.223E-07
  29 ClONO2 + hv -> ClO + NO2                  5.209E-08
  30 CH3Cl + hv -> Products                    6.990E-38
  31 CCl2O + hv -> Products                    1.287E-35
  32 CCl4 + hv -> Products                     3.371E-35
  33 CClFO + hv -> Products                    7.878E-36
  34 CF2O + hv -> Products                     2.067E-37
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     8.581E-37
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     4.942E-38
  37 CF3CF2Cl (CFC-115) + hv -> Products       0.000E+00
  38 CCl3F (CFC-11) + hv -> Products           8.191E-36
  39 CCl2F2 (CFC-12) + hv -> Products          3.556E-37
  40 CH3CCl3 + hv -> Products                  1.384E-35
  41 CF3CHCl2 (HCFC-123) + hv -> Products      7.326E-37
  42 CF3CHFCl (HCFC-124) + hv -> Products      7.437E-39
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     1.047E-36
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     6.254E-39
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  1.683E-36
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  7.990E-39
  47 CHClF2 (HCFC-22) + hv -> Products         0.000E+00
  48 BrONO2 + hv -> Br + NO3                   5.664E-06
  49 BrONO2 + hv -> BrO + NO2                  1.387E-05
  50 CH3Br + hv -> Products                    5.555E-35
  51 CHBr3                                     4.399E-09
  52 CF3Br (Halon-1301) + hv -> Products       1.058E-35
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  9.342E-35
  54 CF2Br2 (Halon-1202) + hv -> Products      2.515E-13
  55 CF2BrCl (Halon-1211) + hv -> Products     1.490E-22
  56 Cl2 + hv -> Cl + Cl                       3.278E-05""",
75: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   75.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          4.802E-39
   2 O3 -> O2 + O(1D)                          1.214E-06
   3 O3 -> O2 + O(3P)                          1.897E-04
   4 NO2 -> NO + O(3P)                         2.372E-03
   5 NO3 -> NO + O2                            1.076E-02
   6 NO3 -> NO2 + O(3P)                        7.840E-02
   7 N2O5 -> NO3 + NO + O(3P)                  4.295E-13
   8 N2O5 -> NO3 + NO2                         7.553E-06
   9 N2O + hv -> N2 + O(1D)                    3.108E-36
  10 HO2 + hv -> OH + O                        1.478E-33
  11 H2O2 -> 2 OH                              1.140E-06
  12 HNO2 -> OH + NO                           4.971E-04
  13 HNO3 -> OH + NO2                          6.171E-08
  14 HNO4 -> HO2 + NO2                         4.503E-07
  15 CH2O -> H + HCO                           4.098E-06
  16 CH2O -> H2 + CO                           8.872E-06
  17 CH3CHO -> CH3 + HCO                       3.939E-07
  18 CH3CHO -> CH4 + CO                        1.988E-15
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     1.842E-06
  21 CHOCHO -> products                        2.303E-05
  22 CH3COCHO -> products                      2.964E-05
  23 CH3COCH3                                  4.399E-08
  24 CH3OOH -> CH3O + OH                       9.225E-07
  25 CH3ONO2 -> CH3O+NO2                       8.578E-08
  26 PAN + hv -> products                      9.552E-08
  27 ClOO + hv -> Products                     2.487E-32
  28 ClONO2 + hv -> Cl + NO3                   9.259E-06
  29 ClONO2 + hv -> ClO + NO2                  1.178E-06
  30 CH3Cl + hv -> Products                    6.228E-37
  31 CCl2O + hv -> Products                    1.077E-34
  32 CCl4 + hv -> Products                     2.869E-34
  33 CClFO + hv -> Products                    6.656E-35
  34 CF2O + hv -> Products                     1.788E-36
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     7.477E-36
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     4.311E-37
  37 CF3CF2Cl (CFC-115) + hv -> Products       0.000E+00
  38 CCl3F (CFC-11) + hv -> Products           7.097E-35
  39 CCl2F2 (CFC-12) + hv -> Products          3.121E-36
  40 CH3CCl3 + hv -> Products                  1.196E-34
  41 CF3CHCl2 (HCFC-123) + hv -> Products      6.366E-36
  42 CF3CHFCl (HCFC-124) + hv -> Products      6.427E-38
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     9.137E-36
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     5.453E-38
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  1.457E-35
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  7.973E-38
  47 CHClF2 (HCFC-22) + hv -> Products         0.000E+00
  48 BrONO2 + hv -> Br + NO3                   1.009E-04
  49 BrONO2 + hv -> BrO + NO2                  2.471E-04
  50 CH3Br + hv -> Products                    4.687E-34
  51 CHBr3                                     1.234E-07
  52 CF3Br (Halon-1301) + hv -> Products       8.898E-35
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  7.828E-34
  54 CF2Br2 (Halon-1202) + hv -> Products      3.071E-12
  55 CF2BrCl (Halon-1211) + hv -> Products     8.945E-22
  56 Cl2 + hv -> Cl + Cl                       5.674E-04""",
60: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   60.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          3.946E-38
   2 O3 -> O2 + O(1D)                          7.319E-06
   3 O3 -> O2 + O(3P)                          3.124E-04
   4 NO2 -> NO + O(3P)                         5.286E-03
   5 NO3 -> NO + O2                            1.697E-02
   6 NO3 -> NO2 + O(3P)                        1.296E-01
   7 N2O5 -> NO3 + NO + O(3P)                  1.470E-11
   8 N2O5 -> NO3 + NO2                         2.065E-05
   9 N2O + hv -> N2 + O(1D)                    2.567E-35
  10 HO2 + hv -> OH + O                        1.210E-32
  11 H2O2 -> 2 OH                              3.215E-06
  12 HNO2 -> OH + NO                           1.136E-03
  13 HNO3 -> OH + NO2                          2.249E-07
  14 HNO4 -> HO2 + NO2                         1.656E-06
  15 CH2O -> H + HCO                           1.284E-05
  16 CH2O -> H2 + CO                           2.288E-05
  17 CH3CHO -> CH3 + HCO                       1.658E-06
  18 CH3CHO -> CH4 + CO                        6.753E-14
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     6.643E-06
  21 CHOCHO -> products                        4.670E-05
  22 CH3COCHO -> products                      6.407E-05
  23 CH3COCH3                                  1.864E-07
  24 CH3OOH -> CH3O + OH                       2.524E-06
  25 CH3ONO2 -> CH3O+NO2                       3.182E-07
  26 PAN + hv -> products                      2.983E-07
  27 ClOO + hv -> Products                     1.282E-31
  28 ClONO2 + hv -> Cl + NO3                   2.221E-05
  29 ClONO2 + hv -> ClO + NO2                  3.518E-06
  30 CH3Cl + hv -> Products                    5.158E-36
  31 CCl2O + hv -> Products                    8.832E-34
  32 CCl4 + hv -> Products                     2.364E-33
  33 CClFO + hv -> Products                    5.476E-34
  34 CF2O + hv -> Products                     1.476E-35
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     6.172E-35
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     3.559E-36
  37 CF3CF2Cl (CFC-115) + hv -> Products       0.000E+00
  38 CCl3F (CFC-11) + hv -> Products           5.857E-34
  39 CCl2F2 (CFC-12) + hv -> Products          2.575E-35
  40 CH3CCl3 + hv -> Products                  9.872E-34
  41 CF3CHCl2 (HCFC-123) + hv -> Products      5.255E-35
  42 CF3CHFCl (HCFC-124) + hv -> Products      5.298E-37
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     7.545E-35
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     4.501E-37
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  1.203E-34
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  6.354E-37
  47 CHClF2 (HCFC-22) + hv -> Products         0.000E+00
  48 BrONO2 + hv -> Br + NO3                   2.237E-04
  49 BrONO2 + hv -> BrO + NO2                  5.477E-04
  50 CH3Br + hv -> Products                    3.855E-33
  51 CHBr3                                     4.471E-07
  52 CF3Br (Halon-1301) + hv -> Products       7.314E-34
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  6.414E-33
  54 CF2Br2 (Halon-1202) + hv -> Products      8.625E-11
  55 CF2BrCl (Halon-1211) + hv -> Products     1.221E-20
  56 Cl2 + hv -> Cl + Cl                       1.333E-03""",
45: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   45.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          5.886E-34
   2 O3 -> O2 + O(1D)                          1.885E-05
   3 O3 -> O2 + O(3P)                          3.869E-04
   4 NO2 -> NO + O(3P)                         7.423E-03
   5 NO3 -> NO + O2                            2.024E-02
   6 NO3 -> NO2 + O(3P)                        1.573E-01
   7 N2O5 -> NO3 + NO + O(3P)                  3.010E-10
   8 N2O5 -> NO3 + NO2                         3.351E-05
   9 N2O + hv -> N2 + O(1D)                    3.000E-31
  10 HO2 + hv -> OH + O                        1.981E-28
  11 H2O2 -> 2 OH                              5.319E-06
  12 HNO2 -> OH + NO                           1.621E-03
  13 HNO3 -> OH + NO2                          4.378E-07
  14 HNO4 -> HO2 + NO2                         3.191E-06
  15 CH2O -> H + HCO                           2.232E-05
  16 CH2O -> H2 + CO                           3.549E-05
  17 CH3CHO -> CH3 + HCO                       3.414E-06
  18 CH3CHO -> CH4 + CO                        1.547E-12
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     1.266E-05
  21 CHOCHO -> products                        6.292E-05
  22 CH3COCHO -> products                      8.868E-05
  23 CH3COCH3                                  4.001E-07
  24 CH3OOH -> CH3O + OH                       4.095E-06
  25 CH3ONO2 -> CH3O+NO2                       6.231E-07
  26 PAN + hv -> products                      5.296E-07
  27 ClOO + hv -> Products                     2.807E-25
  28 ClONO2 + hv -> Cl + NO3                   3.311E-05
  29 ClONO2 + hv -> ClO + NO2                  6.052E-06
  30 CH3Cl + hv -> Products                    5.060E-32
  31 CCl2O + hv -> Products                    1.372E-29
  32 CCl4 + hv -> Products                     3.204E-29
  33 CClFO + hv -> Products                    7.840E-30
  34 CF2O + hv -> Products                     1.800E-31
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     7.230E-31
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     4.111E-32
  37 CF3CF2Cl (CFC-115) + hv -> Products       0.000E+00
  38 CCl3F (CFC-11) + hv -> Products           7.086E-30
  39 CCl2F2 (CFC-12) + hv -> Products          2.908E-31
  40 CH3CCl3 + hv -> Products                  1.211E-29
  41 CF3CHCl2 (HCFC-123) + hv -> Products      6.209E-31
  42 CF3CHFCl (HCFC-124) + hv -> Products      6.607E-33
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     8.676E-31
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     5.253E-33
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  1.461E-30
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  4.875E-33
  47 CHClF2 (HCFC-22) + hv -> Products         0.000E+00
  48 BrONO2 + hv -> Br + NO3                   3.171E-04
  49 BrONO2 + hv -> BrO + NO2                  7.763E-04
  50 CH3Br + hv -> Products                    5.571E-29
  51 CHBr3                                     8.756E-07
  52 CF3Br (Halon-1301) + hv -> Products       1.078E-29
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  1.659E-28
  54 CF2Br2 (Halon-1202) + hv -> Products      6.557E-10
  55 CF2BrCl (Halon-1211) + hv -> Products     2.476E-17
  56 Cl2 + hv -> Cl + Cl                       1.948E-03""",
30: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   30.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          9.140E-30
   2 O3 -> O2 + O(1D)                          3.169E-05
   3 O3 -> O2 + O(3P)                          4.350E-04
   4 NO2 -> NO + O(3P)                         8.827E-03
   5 NO3 -> NO + O2                            2.220E-02
   6 NO3 -> NO2 + O(3P)                        1.740E-01
   7 N2O5 -> NO3 + NO + O(3P)                  1.353E-09
   8 N2O5 -> NO3 + NO2                         4.347E-05
   9 N2O + hv -> N2 + O(1D)                    4.703E-27
  10 HO2 + hv -> OH + O                        3.085E-24
  11 H2O2 -> 2 OH                              6.972E-06
  12 HNO2 -> OH + NO                           1.943E-03
  13 HNO3 -> OH + NO2                          6.325E-07
  14 HNO4 -> HO2 + NO2                         4.565E-06
  15 CH2O -> H + HCO                           2.994E-05
  16 CH2O -> H2 + CO                           4.461E-05
  17 CH3CHO -> CH3 + HCO                       5.043E-06
  18 CH3CHO -> CH4 + CO                        8.316E-12
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     1.796E-05
  21 CHOCHO -> products                        7.371E-05
  22 CH3COCHO -> products                      1.050E-04
  23 CH3COCH3                                  6.151E-07
  24 CH3OOH -> CH3O + OH                       5.306E-06
  25 CH3ONO2 -> CH3O+NO2                       9.018E-07
  26 PAN + hv -> products                      7.261E-07
  27 ClOO + hv -> Products                     1.445E-21
  28 ClONO2 + hv -> Cl + NO3                   4.094E-05
  29 ClONO2 + hv -> ClO + NO2                  8.125E-06
  30 CH3Cl + hv -> Products                    8.050E-28
  31 CCl2O + hv -> Products                    2.139E-25
  32 CCl4 + hv -> Products                     4.966E-25
  33 CClFO + hv -> Products                    1.216E-25
  34 CF2O + hv -> Products                     2.805E-27
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     1.136E-26
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     6.435E-28
  37 CF3CF2Cl (CFC-115) + hv -> Products       0.000E+00
  38 CCl3F (CFC-11) + hv -> Products           1.107E-25
  39 CCl2F2 (CFC-12) + hv -> Products          4.610E-27
  40 CH3CCl3 + hv -> Products                  1.888E-25
  41 CF3CHCl2 (HCFC-123) + hv -> Products      9.680E-27
  42 CF3CHFCl (HCFC-124) + hv -> Products      1.036E-28
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     1.358E-26
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     8.255E-29
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  2.283E-26
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  1.110E-28
  47 CHClF2 (HCFC-22) + hv -> Products         0.000E+00
  48 BrONO2 + hv -> Br + NO3                   3.806E-04
  49 BrONO2 + hv -> BrO + NO2                  9.318E-04
  50 CH3Br + hv -> Products                    8.640E-25
  51 CHBr3                                     1.275E-06
  52 CF3Br (Halon-1301) + hv -> Products       1.672E-25
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  1.831E-24
  54 CF2Br2 (Halon-1202) + hv -> Products      1.893E-09
  55 CF2BrCl (Halon-1211) + hv -> Products     1.239E-15
  56 Cl2 + hv -> Cl + Cl                       2.371E-03""",
15: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   15.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          9.120E-28
   2 O3 -> O2 + O(1D)                          4.139E-05
   3 O3 -> O2 + O(3P)                          4.626E-04
   4 NO2 -> NO + O(3P)                         9.627E-03
   5 NO3 -> NO + O2                            2.329E-02
   6 NO3 -> NO2 + O(3P)                        1.832E-01
   7 N2O5 -> NO3 + NO + O(3P)                  2.791E-09
   8 N2O5 -> NO3 + NO2                         4.966E-05
   9 N2O + hv -> N2 + O(1D)                    4.714E-25
  10 HO2 + hv -> OH + O                        3.086E-22
  11 H2O2 -> 2 OH                              8.006E-06
  12 HNO2 -> OH + NO                           2.127E-03
  13 HNO3 -> OH + NO2                          7.653E-07
  14 HNO4 -> HO2 + NO2                         5.489E-06
  15 CH2O -> H + HCO                           3.474E-05
  16 CH2O -> H2 + CO                           5.004E-05
  17 CH3CHO -> CH3 + HCO                       6.154E-06
  18 CH3CHO -> CH4 + CO                        1.965E-11
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     2.148E-05
  21 CHOCHO -> products                        7.998E-05
  22 CH3COCHO -> products                      1.144E-04
  23 CH3COCH3                                  7.711E-07
  24 CH3OOH -> CH3O + OH                       6.056E-06
  25 CH3ONO2 -> CH3O+NO2                       1.092E-06
  26 PAN + hv -> products                      8.550E-07
  27 ClOO + hv -> Products                     7.855E-20
  28 ClONO2 + hv -> Cl + NO3                   4.564E-05
  29 ClONO2 + hv -> ClO + NO2                  9.454E-06
  30 CH3Cl + hv -> Products                    8.136E-26
  31 CCl2O + hv -> Products                    2.140E-23
  32 CCl4 + hv -> Products                     4.946E-23
  33 CClFO + hv -> Products                    1.212E-23
  34 CF2O + hv -> Products                     2.802E-25
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     1.140E-24
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     6.444E-26
  37 CF3CF2Cl (CFC-115) + hv -> Products       1.144E-36
  38 CCl3F (CFC-11) + hv -> Products           1.107E-23
  39 CCl2F2 (CFC-12) + hv -> Products          4.652E-25
  40 CH3CCl3 + hv -> Products                  1.886E-23
  41 CF3CHCl2 (HCFC-123) + hv -> Products      9.672E-25
  42 CF3CHFCl (HCFC-124) + hv -> Products      1.039E-26
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     1.359E-24
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     8.288E-27
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  2.284E-24
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  1.320E-26
  47 CHClF2 (HCFC-22) + hv -> Products         7.774E-37
  48 BrONO2 + hv -> Br + NO3                   4.175E-04
  49 BrONO2 + hv -> BrO + NO2                  1.022E-03
  50 CH3Br + hv -> Products                    8.614E-23
  51 CHBr3                                     1.551E-06
  52 CF3Br (Halon-1301) + hv -> Products       1.666E-23
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  1.658E-22
  54 CF2Br2 (Halon-1202) + hv -> Products      3.214E-09
  55 CF2BrCl (Halon-1211) + hv -> Products     7.801E-15
  56 Cl2 + hv -> Cl + Cl                       2.617E-03""",
0: """ INPUT PARAMETERS:

 RADIATION SCHEME: 2 streams

 w-grid: 141  120.  735.
 z-grid: 81  0.  80.
  measurement point: index 1 altitude=  0.
 DATAE1/SUN/susim_hi.flx
 DATAE1/SUN/atlas3_1994_317_a.dat        
 DATAE1/SUN/neckel.flx                   
 air temperature: USSA, 1976
 air density: USSA, 1976
     old sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     new sea level air column =  2.1518E+25 # cm-2  =  1014.94 mbar
     surface air column =  2.1518E+25 # cm-2  =  1014.94 mbar
 ozone profile: USSA, 1976
     old O3 Column =  9.3811E+18 # cm-2  =   349.13  Dobson Units 
     new O3 Column =  8.0610E+18 # cm-2  =   300.00  Dobson Units 
 SO2:  1 ppb in lowest 1 km, 0 above
     old SO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new SO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 NO2:  1 ppb in lowest 1 km, 0 above
     old NO2 Column =  2.6900E+15 # cm-2  =     0.10  Dobson Units 
     new NO2 Column =  0.0000E+00 # cm-2  =     0.00  Dobson Units 
 Cloud:  4levels, tot opt. dep. =   0.
 aerosols:  Elterman (1968)
   Total aerosol od at 340 nm =   0.380021363
 wavelength-independent albedo =   0.100000001
 idate =  630 esfact =   0.966632068
 solar zenith angle =   0.
 
PHOTOLYSIS RATES (1/sec):

   1 O2 + hv -> O + O                          3.602E-27
   2 O3 -> O2 + O(1D)                          4.496E-05
   3 O3 -> O2 + O(3P)                          4.717E-04
   4 NO2 -> NO + O(3P)                         9.887E-03
   5 NO3 -> NO + O2                            2.364E-02
   6 NO3 -> NO2 + O(3P)                        1.862E-01
   7 N2O5 -> NO3 + NO + O(3P)                  3.473E-09
   8 N2O5 -> NO3 + NO2                         5.175E-05
   9 N2O + hv -> N2 + O(1D)                    1.929E-24
  10 HO2 + hv -> OH + O                        1.216E-21
  11 H2O2 -> 2 OH                              8.355E-06
  12 HNO2 -> OH + NO                           2.187E-03
  13 HNO3 -> OH + NO2                          8.119E-07
  14 HNO4 -> HO2 + NO2                         5.811E-06
  15 CH2O -> H + HCO                           3.636E-05
  16 CH2O -> H2 + CO                           5.184E-05
  17 CH3CHO -> CH3 + HCO                       6.543E-06
  18 CH3CHO -> CH4 + CO                        2.564E-11
  19 CH3CHO -> CH3CO + H                       0.000E+00
  20 C2H5CHO -> C2H5 + HCO                     2.270E-05
  21 CHOCHO -> products                        8.203E-05
  22 CH3COCHO -> products                      1.175E-04
  23 CH3COCH3                                  8.275E-07
  24 CH3OOH -> CH3O + OH                       6.309E-06
  25 CH3ONO2 -> CH3O+NO2                       1.158E-06
  26 PAN + hv -> products                      8.996E-07
  27 ClOO + hv -> Products                     2.573E-19
  28 ClONO2 + hv -> Cl + NO3                   4.720E-05
  29 ClONO2 + hv -> ClO + NO2                  9.908E-06
  30 CH3Cl + hv -> Products                    3.504E-25
  31 CCl2O + hv -> Products                    8.477E-23
  32 CCl4 + hv -> Products                     1.958E-22
  33 CClFO + hv -> Products                    4.792E-23
  34 CF2O + hv -> Products                     1.130E-24
  35 CF2ClCFCl2 (CFC-113) + hv -> Products     4.707E-24
  36 CF2ClCF2Cl (CFC-114) + hv -> Products     2.674E-25
  37 CF3CF2Cl (CFC-115) + hv -> Products       7.951E-28
  38 CCl3F (CFC-11) + hv -> Products           4.482E-23
  39 CCl2F2 (CFC-12) + hv -> Products          1.981E-24
  40 CH3CCl3 + hv -> Products                  7.583E-23
  41 CF3CHCl2 (HCFC-123) + hv -> Products      3.979E-24
  42 CF3CHFCl (HCFC-124) + hv -> Products      4.261E-26
  43 CH3CFCl2 (HCFC-141b) + hv -> Products     5.638E-24
  44 CH3CF2Cl (HCFC-142b) + hv -> Products     3.446E-26
  45 CF3CF2CHCl2 (HCFC-225ca) + hv -> Product  9.298E-24
  46 CF2ClCF2CHFCl (HCFC-225cb) + hv -> Produ  6.782E-26
  47 CHClF2 (HCFC-22) + hv -> Products         5.478E-28
  48 BrONO2 + hv -> Br + NO3                   4.296E-04
  49 BrONO2 + hv -> BrO + NO2                  1.052E-03
  50 CH3Br + hv -> Products                    3.404E-22
  51 CHBr3                                     1.649E-06
  52 CF3Br (Halon-1301) + hv -> Products       6.573E-23
  53 CF2BrCF2Br (Halon-2402) + hv -> Products  6.415E-22
  54 CF2Br2 (Halon-1202) + hv -> Products      3.780E-09
  55 CF2BrCl (Halon-1211) + hv -> Products     1.349E-14
  56 Cl2 + hv -> Cl + Cl                       2.698E-03"""
}

angles = [k for k in tuv4pt1_data.keys()]
angles.sort()
angles = np.array(angles, dtype = 'f')

jvalues_byidx = {}
jvalues_bykey = {}
jidxs = []
jlabels = []
grab = True
for angle in angles:
    for line in tuv4pt1_data[angle].split('\n')[33:]:
        idx = int(line[:4])
        key = line[5:46].strip()
        if grab:
            jidxs.append(idx)
            jlabels.append(key)
        val = eval(line[47:])
        jvalues_byidx.setdefault(idx, []).append(val)
        jvalues_bykey.setdefault(key, []).append(val)
    grab = False

for jvalues in [jvalues_byidx, jvalues_bykey]:
    for k, v in jvalues.items():
        v = jvalues[k] = np.array(v)
        assert(v.size == angles.size)

def TUV_J(idx, zenithangle, scale = 1.):
    """
        idx = TUV 4.1 reaction string (e.g., jlabel) or TUV 4.1 numeric index
        zenithangle = angle of the sun from zenith in degrees
    """
    assert(np.all(zenithangle >= 0))
    if np.all(zenithangle > 100):
        return zenithangle * 0.
    
    try:
        jvals = jvalues_byidx[idx]
    except KeyError:
        try:
            jvals = jvalues_bykey[idx]
        except KeyError:
            raise KeyError('Not in tuv data (idx and jlabels follow) -- idx: %s -- jlabel: %s' % (', '.join([str(i_) for i_ in jvalues_byidx.keys()]), ', '.join(jvalues_bykey.keys())))
    
    return np.interp(zenithangle, angles, jvals)*scale

TUV_J.__doc__ +=  '\n' + '\n'.join(['%s %s' % ik_ for ik_ in zip(jidxs, jlabels)])
TUV_J4pt1 = TUV_J
