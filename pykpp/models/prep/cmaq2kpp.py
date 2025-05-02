from __future__ import print_function
import sys, re
from PseudoNetCDF.cmaqfiles import jtable

def print_usage():
    print("""

  Usage: %s inputpath

    Converts CMAQ formatted reactions to KPP formatted
    reactions. Output is sent to stdout. Photolysis reactions
    must be updated manually after the conversion.
    
    Example:
      %s mech.def > mech.kpp
    
    """ % (sys.argv[0], sys.argv[0]))

jkey2tuvkey = {
    'NO2_IUPAC10': "TUV_J5pt0('NO2 -> NO + O(3P)', THETA);",
    'O3_O3P_IUPAC10': "TUV_J5pt0('O3 -> O2 + O(3P)', THETA);",
    'O3_O1D_IUPAC10': "TUV_J5pt0('O3 -> O2 + O(1D)', THETA);",
    'H2O2_IUPAC10': "TUV_J5pt0('H2O2 -> 2 OH', THETA);",
    'NO3NO2_06': "TUV_J5pt0('NO3 -> NO2 + O(3P)', THETA);",
    'NO3NO_06': "TUV_J5pt0('NO3 -> NO + O2', THETA);",
    # what about N2O5 -> NO3 + NO + O3P; BHH TUV_J5pt0('N2O5 -> NO3 + NO + O(3P)', THETA)+
    'N2O5_IUPAC10': "TUV_J5pt0('N2O5 -> NO3 + NO2', THETA);",
    'HONO_IUPAC10': "TUV_J5pt0('HNO2 -> OH + NO', THETA);",
    'HNO3_IUPAC10': "TUV_J5pt0('HNO3 -> OH + NO2', THETA);",
    'PNA_IUPAC10': "JHNO4_NEAR_IR(TUV_J5pt0('HNO4 -> HO2 + NO2', THETA));",
    'MEPX_IUPAC10': "TUV_J5pt0('CH3OOH -> CH3O + OH', THETA);",
    'PAN_IUPAC10': "TUV_J5pt0('CH3CO(OONO2) -> CH3CO(OO) + NO2', THETA) + TUV_J5pt0('CH3CO(OONO2) -> CH3CO(O) + NO3', THETA);",
    'NTR_IUPAC10': "TUV_J5pt0('CH3CHONO2CH3 -> CH3CHOCH3 + NO2', THETA);",
    'FORM_R_IUPAC13': "TUV_J5pt0('CH2O -> H + HCO', THETA);",
    'FORM_M_IUPAC13': "TUV_J5pt0('CH2O -> H2 + CO', THETA);",
    'ALD2_R_IUPAC13': "TUV_J5pt0('CH3CHO -> CH3 + HCO', THETA)+TUV_J5pt0('CH3CHO -> CH3CO + H', THETA);",
    'ALDX_R_IUPAC13': "TUV_J5pt0('C2H5CHO -> C2H5 + HCO', THETA);",
    'MGLY_IUPAC10': "TUV_J5pt0('CH3COCHO -> CH3CO + HCO', THETA);",
    'ACET_IUPAC10': "TUV_J5pt0('CH3COCH3 -> CH3CO + CH3', THETA);",
    'ISPD': "0.0036 * TUV_J5pt0('CH2=CHCHO -> Products', THETA);",
    'CL2_IUPAC04': "TUV_J5pt0('Cl2 -> Cl + Cl', THETA);",
    'FMCL_IUPAC04': "0; {HOCL_IUPAC04 not sure what to do}",
    'HOCL_IUPAC04': "0; {HOCL_IUPAC04 not sure what to do}",
    'GLY_R_IUPAC13': "TUV_J5pt0('CHOCHO -> HCO + HCO', THETA);",
    'GLYD_IUPAC13': "TUV_J5pt0('CH2(OH)CHO -> Products', THETA);",
    'KET_IUPAC10': "TUV_J5pt0('CH3COCH2CH3 -> CH3CO + CH2CH3', THETA);",
    'CLONO2_1': "TUV_J5pt0('ClONO2 -> ClO + NO2', THETA);",
    'CLONO2_2': "TUV_J5pt0('ClONO2 -> Cl + NO3', THETA);",
    'CLNO2_IUPAC13': "TUV_J5pt0('ClNO2 -> Cl + NO2', THETA);",
    'ACRO_09': "TUV_J5pt0('CH2=CHCHO -> Products', THETA);",
    'BR2_IUPAC10': "TUV_J5pt0('Br2 -> Br + Br', THETA);",
    'HOBR_IUPAC10': "TUV_J5pt0('HOBr -> OH + Br', THETA);",
    'BRO_IUPAC10': "TUV_J5pt0('BrO -> Br + O', THETA);",
    'BRNO2_IUPAC10': " 0; {BRNO2_IUPAC10 not sure what to do}",
    'BRONO2_M_IUPAC10': "TUV_J5pt0('BrONO2 -> BrO + NO2', THETA);",
    'BRONO2_R_IUPAC10': "TUV_J5pt0('rONO2 -> Br + NO3', THETA);",
    'BRCL_IUPAC10': " 0; {BRCL_IUPAC10 not sure what to do}",
    'COHBR_JPL2010': " 0; {COHBR_JPL2010 not sure what to do}",
    'MB3_IUPAC10': "TUV_J5pt0('CHBr3 -> Products', THETA);",
    'MB2C_BLIDE98': " 0; {MB2C_BLIDE98 not sure what to do}",
    'MBC2_BLIDE98': " 0; {MBC2_BLIDE98 not sure what to do}",
    'I2_IUPAC10': " 0; {I2_IUPAC10 not sure what to do}",
    'HOI_IUPAC10': " 0; {HOI_IUPAC10 not sure what to do}",
    'IO_IUPAC10': " 0; {IO_IUPAC10 not sure what to do}",
    'OIO_06': " 0; {OIO_06 not sure what to do}",
    'INO_06': " 0; {INO_06 not sure what to do}",
    'INO2_06': "0; {INO2_06 not sure what to do}",
    'IONO2_06': "0; {IONO2_06 not sure what to do}",
    'ICL_IUPAC10': "0; {ICL_IUPAC10 not sure what to do}",
    'IBR_IUPAC10': "0; {IBR_IUPAC10 not sure what to do}",
    'CH3I_IUPAC10': "0; {CH3I_IUPAC10 not sure what to do}",
    'CH3I_IUPAC10': "0; {CH3I_IUPAC10 not sure what to do}",
    'MI2_IUPAC10': "0; {MI2_IUPAC10 not sure what to do}",
    'MIB_IUPAC10': "0; {MIB_IUPAC10 not sure what to do}",
    'MIC_IUPAC10': "0; {MIC_IUPAC10 not sure what to do}",
    'HPALD': "0; {HPALD not sure what to do}",
    'IC3ONO2': "0; {IC3ONO2 not sure what to do}",
    # '': "TUV_J5pt0('', THETA)",
}

if __name__ == '__main__':
    if len(sys.argv) != 2:
        print_usage()
        exit()
    try:
        constant_spcs = "O2 M N2 H2 H2O".split()

        constant_spc_re = re.compile(r'((?:(?:\+\s*)?\b(?:' + '|'.join(constant_spcs) + r')\b(?:\s*\+\s*)?)+)', re.M)
        constant_spc_rct = re.compile(r'\{((?:(?:\+\s*)?\b(?:' + '|'.join(constant_spcs) + r')\b(?:\s*\+\s*)?)+)\}.*?=.*?\n', re.M)

        newtxt = mechtxt = open(sys.argv[1], 'r').read()

        comment = re.compile(r'![^\n]*\n', re.M)
        newtxt = comment.sub('', newtxt)

        newline = re.compile(r'\n[ ]+', re.M)
        newtxt = newline.sub(' ', newtxt)

        ratedelim = re.compile('([&^@;])', re.M)
        newtxt = ratedelim.sub(r' \1 ', newtxt)

        spaces = re.compile('[ ]+', re.M)
        newtxt = spaces.sub(' ', newtxt)

        newline = re.compile(r'\n+', re.M)
        newtxt = newline.sub(r'\n', newtxt)



        adddecimal = re.compile('(^-?\d+$)')

        def fillwithzero(template, matcho):
            groups = [g or '0.0' for g in matcho.groups()]
            groups = [adddecimal.sub(r'\1.0', g) for g in groups]
            return template % tuple(groups)

        def hetk(matcho):
            groups = matcho.groups()
            if len(groups) > 2:
                raise ValueError("Oops too many HET groups")
            retval = ': 0 {{{}~{}}};\n'.format(*groups)
            return retval
            
        def tuvj(matcho):
            groups = matcho.groups()
            if len(groups) > 2:
                raise ValueError("Oops too many TUV groups")
            jkey = groups[1]
            jtxt = jkey2tuvkey.get(groups[1], None)
            if jtxt is None:
                jtxt = '?TUV_J5pt0("{}", THETA)'.format(groups[1])
                print(jtxt)
            jcmt = '{}/{}'.format(*groups)
            retval = f'{jtxt} {{{jcmt}}}'
    
            if groups[0] is not None and groups[0] != '1.0':
                retval = ('%s * ' % groups[0]) + retval
            retval = ': ' + retval + '\n'
            return retval        
    
        cmaq_1to4 = re.compile('(?<!%\d\s)#\s?(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s([\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?)){1}\n', re.M)
        newtxt = cmaq_1to4.sub(lambda m: fillwithzero(': CMAQ_1to4(%s, %s, %s);\n', m), newtxt)

        cmaq_8 = re.compile(r'%\d\s#\s?(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))\n')
        newtxt = cmaq_8.sub(lambda m: fillwithzero(': CMAQ_8(%s, %s, %s, %s, %s, %s);\n', m), newtxt)

        cmaq_9 = re.compile(r'%\d\s#\s?(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))\n')
        newtxt = cmaq_9.sub(lambda m: fillwithzero(': CMAQ_9(%s, %s, %s, %s);\n', m).replace('e', 'd').replace('E', 'd'), newtxt)

        cmaq_10 = re.compile(r'(?<!%\d\s)#\s?(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s([\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s([\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s(?:[\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s(?:[\d.Ee\-+]*))?(?:\s[&;]\s?))\n')
        newtxt = cmaq_10.sub(lambda m: fillwithzero(': CMAQ_10(%s, %s, %s, %s, %s, %s, %s, %s);\n', m), newtxt)

        tuv_j = re.compile(r'(?<!%\d\s)#\s?([\d.Ee\-+]*)\s*/\s*<\s*(\S+)\s*>\s*;\s*\n', re.M)
        newtxt = tuv_j.sub(tuvj, newtxt)

        het_k = re.compile(r'(?<!%\d\s)#\s?([\d.Ee\-+]*)\s*~\s*<\s*(\S+)\s*>\s*;\s*\n', re.M)
        newtxt = het_k.sub(hetk, newtxt)

        label = re.compile(r'^<(\s*\S*\s*)>', re.M)
        newtxt = label.sub(r'<\1>', newtxt)

        header = re.compile(r'.*REACTIONS\s*\[CM\] =\s*\n', re.M | re.S)
        newtxt = header.sub('', newtxt)

        footer = re.compile(r'\n(endmech|END).*', re.M | re.I | re.S)
        newtxt = footer.sub('', newtxt)

        newlabels = [re.compile('<(\s*\S+\s*)>').search(line).groups()[0].strip() for line in newtxt.split('\n') if line != '']

        def deref_rate(matcho):
            groups = matcho.groups()
            rxnkey = groups[1].strip()
            try:
                ratetxt = rconst[rxnkey].replace(';', '')
                ratetxt = ': %s * (%s) {%s};\n' % (groups[0], ratetxt, groups[1])
            except Exception:
                ratetxt = '# %s*K<%s>;\n' % (groups[0], groups[1])
                print(ratetxt, groups, file=sys.stderr)
            return ratetxt

        rateref = re.compile('(?<!%\d\s)#\s([\d.Ee\-+]*)\*K<(\S+)>\s*;\s*\n', re.M)
        rconst = {k.strip(): v.strip().replace(';', '') for k, v in re.compile('<(.+?)>.+?:(.+)\n').findall(newtxt)}
        newtxt = rateref.sub(deref_rate, newtxt)
        rconst = {k.strip(): v.strip().replace(';', '') for k, v in re.compile('<(.+?)>.+?:(.+)\n').findall(newtxt)}
        newtxt = rateref.sub(deref_rate, newtxt)
        rconst = {k.strip(): v.strip().replace(';', '') for k, v in re.compile('<(.+?)>.+?:(.+)\n').findall(newtxt)}
        newtxt = rateref.sub(deref_rate, newtxt)

        newtxt = '#EQUATIONS\n' + newtxt

        newtxt = re.compile('(?<![eE])([+])(?!\d)', re.M).sub(r' \1 ', newtxt)
        newtxt = re.compile('[ ]+', re.M).sub(' ', newtxt)

        def removeasterisks(matcho):
            groups = matcho.groups()
            if len(groups) != 1:
                raise ValueError("Should just be 1")
    
            return groups[0].replace('*', ' ')

        rxn_no_rate = re.compile('(?<=})([^:]+)(?=:)', re.M)
        newtxt = rxn_no_rate.sub(removeasterisks, newtxt)

        max_length = max([len(lab) for lab in newlabels])
        for i in range(1, max_length + 1):
            newtxt = re.compile('^{\s*(\S{%d})}' % i, re.M).sub('{' + ((max_length - i) * ' ') + r'\1} ', newtxt)

        def scaleconstants(matcho):
            groups = [[vs.strip() for vs in v.split('+') if vs.strip() != ''] for v in matcho.groups() if v is not None]
            outtext = matcho.group()

            for group in groups:
                outtext = outtext.replace(':', ': ' + ' * '.join(group) + ' *')
            return outtext
    
        newtxt = constant_spc_re.sub(r'{\1}', newtxt)
        newtxt = constant_spc_rct.sub(scaleconstants, newtxt)
        newtxt = newtxt.replace('= :', '= DUMMY :')
        def rctprd(matcho):
            groups = matcho.groups()
            if len(groups) > 3:
                raise ValueError("Oops too many rctprd groups")
            retval = groups[0] + groups[1].replace('*', ' ') + groups[2]
            return retval
            

        newtxt = re.compile('(>)(.+?)(:)').sub(rctprd, newtxt)
        print('Has CMAQ_1to4', 'CMAQ_1to4' in newtxt, file=sys.stderr)
        print('Has CMAQ_8', 'CMAQ_8' in newtxt, file=sys.stderr)
        print('Has CMAQ_9', 'CMAQ_9' in newtxt, file=sys.stderr)
        print('Has CMAQ_10', 'CMAQ_10' in newtxt, file=sys.stderr)
        print(newtxt, file=sys.stdout)
    except KeyError as e:
        print(str(e))
        print_usage()
        import pdb; pdb.set_trace()
        raise e
