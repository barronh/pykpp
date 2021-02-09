__all__ = ['MCMJ']
from warnings import warn
from pykpp.tuv.tuv5pt0 import TUV_J


def MCMJ(idx, THETA):
    """
    Parameters:
        idx - index from MCM export (v3.3.1)
        THETA - zenith angle in degrees

    Returns:
        photolysis frequency (s**-1) from TUV_J
        for closest surrogate.
    """
    if idx == 1:
        """
        MCM
           1 = O3 -> O(1D) +O2
        TUV v5.0
           2 = O3 -> O2 + O(1D)
        """
        return TUV_J(2, THETA)
    elif idx == 2:
        """
        MCM
           2 = O3 - O(3P) +O2
        TUV v5.0
           3 = O3 -> O2 + O(3P)
        """
        return TUV_J(3, THETA)
    elif idx == 3:
        """
        MCM
           3 = H2O2 -> OH + OH
        TUV v5.0
           5 = H2O2 -> 2 OH
        """
        return TUV_J(5, THETA)
    elif idx == 4:
        """
        MCM
           4 = NO2 -> NO + O(3P)
        TUV v5.0
           6 = NO2 -> NO + O(3P)
        """
        return TUV_J(6, THETA)
    elif idx == 5:
        """
        MCM
           5 = NO3 -> NO + O2
        TUV v5.0
           7 = NO3 -> NO + O2
        """
        return TUV_J(7, THETA)
    elif idx == 6:
        """
        MCM
           6 = NO3 -> NO2 + O(3P)
        TUV v5.0
           8 = NO3 -> NO2 + O(3P)
        """
        return TUV_J(8, THETA)
    elif idx == 7:
        """
        MCM
           7 = HONO -> NO + OH
        TUV v5.0
          12 = HNO2 -> OH + NO
        """
        return TUV_J(12, THETA)
    elif idx == 8:
        """
        MCM
           8 = HNO3 -> NO2 + OH
        TUV v5.0
          13 = HNO3 -> OH + NO2
        """
        return TUV_J(13, THETA)
    elif idx == 11:
        """
        MCM
          11 = HCHO -> H + HCO
        TUV v5.0
          17 = CH2O -> H + HCO
        """
        return TUV_J(17, THETA)
    elif idx == 12:
        """
        MCM
          12 = HCHO -> H2 + CO
        TUV v5.0
          18 = CH2O -> H2 + CO
        """
        return TUV_J(18, THETA)
    elif idx == 13:
        """
        MCM
          13 = CH3CHO -> CH3 + HCO
        TUV v5.0
          19 = CH3CHO -> CH3 + HCO
          20 = CH3CHO -> CH4 + CO
          21 = CH3CHO -> CH3CO + H

        Notes: Assuming the MCM uses same total QY with only one product set.
        """
        return TUV_J(19, THETA)+TUV_J(20, THETA)+TUV_J(21, THETA)
    elif idx == 14:
        """
        MCM
          14 = C2H5CHO -> C2H5 + HCO
        TUV v5.0
          22 = C2H5CHO -> C2H5 + HCO
        """
        return TUV_J(22, THETA)
    elif idx == 15:
        """
        MCM
          15 = C3H7CHO -> n-C3H7 + HCO (QY: 0.21)
          16 = C3H7CHO -> C2H4 + CH3CHO (QY: 0.1)
        TUV v5.0
          22 = C2H5CHO -> C2H5 + HCO

        Notes: Using C3 aldehyde as a surrogate.
        """
        return .667*TUV_J(22, THETA)
    elif idx == 16:
        """
        MCM
          15 = C3H7CHO -> n-C3H7 + HCO (QY: 0.21)
          16 = C3H7CHO -> C2H4 + CH3CHO (QY: 0.1)
        TUV v5.0
          22 = C2H5CHO -> C2H5 + HCO

        Notes: Using C3 aldehyde as a surrogate.
        """
        return .333*TUV_J(22, THETA)
    elif idx == 17:
        """
        MCM
          17 = IPRCHO -> n-C4H9 + HCO
        TUV v5.0
          22 = C2H5CHO -> C2H5 + HCO

        Notes: Using C3 aldehyde as a surrogate.
        """
        return TUV_J(22, THETA)
    elif idx == 18:
        """
        MCM
          18 = MACR -> CH2=CCH3 + HCO (QY = 1.9e-3)
          19 = MACR -> CH2=C(CH3)CO + H (QY = 1.9e-3)
        TUV v5.0
          25 = CH2=C(CH3)CHO -> Products

        Notes: Assuming equal distribution among QY channels.
        """
        return .5*TUV_J(25, THETA)
    elif idx == 19:
        """
        MCM
          18 = MACR -> CH2=CCH3 + HCO (QY = 1.9e-3)
          19 = MACR -> CH2=C(CH3)CO + H (QY = 1.9e-3)
        TUV v5.0
          25 = CH2=C(CH3)CHO -> Products

        Notes: Assuming equal distribution among QY channels.
        """
        return .5*TUV_J(25, THETA)
    elif idx == 20:
        """
        MCM
          20 = C5HPALD1 -> CH3C(CHO)=CHCH2O + OH
        TUV v5.0
          25 = CH2=C(CH3)CHO -> Products

        Notes: This compound cleaves at the peroxide, so using CH3OOH as a
               surrogate.
        """
        return TUV_J(31, THETA)
    elif idx == 21:
        """
        MCM
          21 = CH3COCH3 -> CH3CO + CH3
        TUV v5.0
          26 = CH3COCH3 -> CH3CO + CH3
        """
        return TUV_J(26, THETA)
    elif idx == 22:
        """
        MCM
          22 = MEK -> CH3CO + C2H5
        TUV v5.0
          28 = CH3COCH2CH3 -> CH3CO + CH2CH3
        """
        return TUV_J(28, THETA)
    elif idx == 23:
        """
        MCM
          23 = MVK -> CH3CH=CH2 + CO
          24 = MVK -> CH3CO + CH2=CH
        TUV v5.0
          27 = CH3COCHCH2 -> Products

        Note: Assuming equal split.
        """
        return .5*TUV_J(27, THETA)
    elif idx == 24:
        """
        MCM
          23 = MVK -> CH3CH=CH2 + CO
          24 = MVK -> CH3CO + CH2=CH
        TUV v5.0
          27 = CH3COCHCH2 -> Products

        Note: Assuming equal split.
        """
        return .5*TUV_J(27, THETA)
    elif idx == 31:
        """
        MCM
          31 = GLYOX -> CO + CO + H2
        TUV v5.0
          unavailable for this channel
        """
        return 0.
    elif idx == 32:
        """
        MCM
          32 = GLYOX -> HCHO + CO
        TUV v5.0
          45 = CHOCHO -> CH2O + CO
        """
        return TUV_J(45, THETA)
    elif idx == 33:
        """
        MCM
          33 = GLYOX -> HCO + HCO
        TUV v5.0
          44 = CHOCHO -> HCO + HCO
        """
        return TUV_J(44, THETA)
    elif idx == 34:
        """
        MCM
          34 = MGLYOX -> CH3CO + HCO
        TUV v5.0
          46 = CH3COCHO -> CH3CO + HCO
        """
        return TUV_J(46, THETA)
    elif idx == 35:
        """
        MCM
          35 = BIACET -> CH3CO + CH3CO
        TUV v5.0
          47 = CH3COCOCH3 -> Products
        """
        return TUV_J(47, THETA)
    elif idx == 41:
        """
        MCM
          41 = CH3OOH -> CH3O + OH
        TUV v5.0
          31 = CH3OOH -> CH3O + OH
        """
        return TUV_J(31, THETA)
    elif idx == 51:
        """
        MCM
          51 = CH3NO3 -> CH3O + NO2
        TUV v5.0
          34 = CH3ONO2 -> CH3O + NO2
        """
        return TUV_J(34, THETA)
    elif idx == 52:
        """
        MCM
          52 = C2H5NO3 -> C2H5O + NO2
        TUV v5.0
          35 = CH3CH2ONO2 -> CH3CH2O + NO2
        """
        return TUV_J(35, THETA)
    elif idx == 53:
        """
        MCM
          53 = NC3H7NO3 -> n-C3H7O + NO2
        TUV v5.0
          35 = CH3CH2ONO2 -> CH3CH2O + NO2
        """
        return TUV_J(35, THETA)
    elif idx == 54:
        """
        MCM
          54 = IC3H7NO3 -> CH3C(O.)CH3 + NO2
        TUV v5.0
          36 = CH3CHONO2CH3 -> CH3CHOCH3 + NO2
        """
        return TUV_J(36, THETA)
    elif idx == 55:
        """
        MCM
          55 = TC4H9NO3 -> t-C4H9O + NO2
        TUV v5.0
          39 = C(CH3)3(ONO2) -> C(CH3)3(O.) + NO2
        """
        return TUV_J(39, THETA)
    elif idx == 56:
        """
        MCM
          56 = NOA -> CH3C(O)CH2(O.) + NO2
          56 = NOA -> CH3CO + HCHO + NO2
        TUV v5.0
          38 = CH3COCH2(ONO2) -> CH3COCH2(O.) + NO2
        """
        return TUV_J(38, THETA)
    elif idx == 57:
        """
        57 = existed in v3.1, but has been removed
        MCM v3.1
          C51NO32CO = PEN2ONE1O + NO2 :  J(56)  ;
          C51NO32CO = C3H7CO3 + HCHO + NO2 :  J(57)  ;
          C52NO31CO = C4CHO2O + NO2 :  J(56)  ;
          C52NO31CO = C3H7CHO + CO + HO2 + NO2 :  J(57)  ;
          NO3CH2CHO = NO2 + HCOCH2O :  J(56)  ;
          NO3CH2CHO = HO2 + CO + HCHO + NO2 :  J(57)  ;
          NOA = CH3CO3 + HCHO + NO2 :  J(57)  ;
          CHOPRNO3 = HO2 + CO + CH3CHO + NO2 :  J(57)  ;
          C51NO324CO = HCHO + CO2C3CO3 + NO2 :  J(57)  ;
          C5NO3CO4OH = HO2CO4C5O + NO2 :  J(56)  ;
          C5NO3CO4OH = IPROPOLO2 + HCHO + CO + NO2 :  J(57)  ;
          C5NO3OAOOH = C5NO3COAO + OH :  J(41)+J(56)+J(57)  ;
          C3NO3COOOH = C3NO3COO + OH :  J(41)+J(56)+J(57)  ;

        """
        warn('mcm TUV has been updated to v3.3.1 and uses J(56) only')
        return 0
    elif idx == 61:
        """
        MCM
          61 = existed in v3.1, but has been removed
        """
        return 0
