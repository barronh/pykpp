from __future__ import print_function
import re
import sys

definitions = dict(
    A3O2 = 'CH3CH2CH2OO',
    ACET = 'CH3C(O)CH3',
    ACTA = 'CH3C(O)OH',
    ALD2 = 'CH3CHO',
    ALK4 = 'CH3CH2CH2CH3',
    ATO2 = 'CH3C(O)CH2O2',
    B3O2 = 'CH3CH(OO)CH3',
    C2H6 = 'C2H6',
    C3H8 = 'C3H8',
    CH2O = 'CH2O',
    CH4 = 'CH4',
    CO = 'CO',
    CO2 = 'CO2',
    EOH = 'C2H5OH',
    ETO2 = 'CH3CH2OO',
    ETP = 'CH3CH2OOH',
    GLYC = 'HOCH2CHO',
    GLYX = 'CHOCHO',
    H2 = 'H2',
    H2O = 'H2O',
    H2O2 = 'H2O2',
    HAC = 'HOCH2C(O)CH3',
    HCOOH = 'HCOOH',
    HNO2 = 'HONO',
    HNO3 = 'HNO3',
    HNO4 = 'HNO4',
    HO2 = 'HO2',
    IALD = 'HOCH2C(CH3)=CHCHO',
    IAP = 'HOCH2C(CH3)(OOH)CH(OH)CHO',
    INO2 = 'O2NOCH2C(OO)(CH3)CH=CH2',
    INPN = 'O2NOCH2C(OOH)(CH3)CH=CH2',
    ISN2 = 'CH2=C(CH3)CH(ONO2)CH2OH',
    ISNO3 = 'RONO2',
    ISNP = 'HOCH2C(OOH)(CH3)CH(ONO2)CH2OH',
    ISOP = 'CH2=C(CH3)CH=CH2',
    KO2 = 'CH3CH2C(O)CH2O2', 
    MACR = 'CH2=C(CH3)CHO',
    MAN2 = 'HOCH2C(ONO2)(CH3)CHO',
    MAO3 = 'CH2=C(CH3)C(O)OO',
    MAOP = 'CH2=C(CH3)C(O)OOH',
    MAP = 'CH3C(O)OOH',
    MCO3 = 'CH3C(O)OO',
    MEK = 'CH3C(O)CH2CH3',
    MGLY = 'CH3COCHO',
    MNO3 = 'CH3ONO2',
    MO2 = 'CH3O2',
    MOH = 'CH3OH',
    MP = 'CH3OOH',
    MRO2 = 'HOCH2C(OO)(CH3)CHO',
    MRP = 'HOCH2C(OOH)(CH3)CHO',
    MVK = 'CH2=CHC(O)CH3',
    N2 = 'N2',
    N2O = 'N2O',
    N2O5 = 'N2O5',
    NH2 = 'NH2',
    NH3 = 'NH3',
    NO = 'NO',
    NO2 = 'NO2',
    NO3 = 'NO3',
    O = 'O',
    O1D = 'O1D',
    O2 = 'O2',
    O2CH2OH = 'O2CH2OH',
    O3 = 'O3',
    OH = 'OH',
    PAN = 'CH3C(O)OONO2',
    PMN = 'CH2=C(CH3)C(O)OONO2',
    PO2 = 'HOCH2CH(OO)CH3',
    PP = 'HOCH2CH(OOH)CH3',
    PPN = 'CH3CH2C(O)OONO2',
    PRN1 = 'O2NOCH2CH(OO)CH3',
    PRPE = 'C3H6',
    PRPN = 'O2NOCH2CH(OOH)CH3',
    R4N1 = 'CH3CH2CH2CHO2NOO',
    R4N2 = 'CH3CH2CH2CH2O2NO',
    R4O2 = 'CH3CH2CH2CH2O2',
    R4P = 'CH3CH2CH2CH2OOH',
    RA3P = 'CH3CH2CH2OOH',
    RB3P = 'CH3CH(OOH)CH3',
    RCHO = 'CH3CH2CHO',
    RCO3 = 'CH3CH2C(O)OO',
    RCOOH = 'C2H5C(O)OH',
    RIO1 = 'HOCH2C(OO)(CH3)CH=CHOH',
    RIO2 = 'HOCH2C(OO)(CH3)CH=CH2',
    RIP = 'HOCH2C(OOH)(CH3)CH=CH2',
    ROH = 'C3H7OH',
    RP = 'CH3CH2C(O)OOH',
    VRO2 = 'HOCH2CH(OO)C(O)CH3',
    VRP = 'HOCH2CH(OOH)C(O)CH3',
    DMS = 'CH3CH3S',
    SO2 = 'SO2',
    SO4 = 'SO4',
    MSA = 'CH4SO3',
    HBrMCO3 = 'HBrCH3C(O)O2',
    HBrA3O2 = 'HBrCH3C(O)CH2O2',
    HC5 = 'C5H12',
    HC5OO = 'C5H11O2',
    MPN = 'CH3OONO3',
    MVKN = 'HOCH2CH(ONO2)C(=O)CH3',
    ISOPNBO2 = 'C5H8NO4O2',
    ISOPNDO2 = 'C5H8NO4O2',
    ATOOH = 'CH3C(O)CH2OOH',
    MACRNO2 = 'C4H5NO5',
    MACRN = 'C4H7NO5',
    DHMOB = 'HOCH2C(CH3)(OH)C(=O)CHO',
    ETHLN = 'CHOCH2ONO2',
    HBrATO2 = 'HBrCH3C(O)CH2O2',
    ISOPND = 'C5H9NO4',
    ISOPNB = 'C5H9NO4',
    MOBAOO = ' HOC(=O)C(CH3)=CHCOOO',
    IEPOX = 'HOCH2C(O)(CH3)CHCH2OH',
    IEPOXOO = 'C5H9O4',# A - C5H9O3, B - C5H9O4 C - C5H9O5
    MAOPO2 = 'CH2OH-CHOO*CH3C(O)OOH',
    MOBA = 'HOC(=O)C(CH3)=CHCHO',
    HBrETO2 = 'HBrCH3CH2OO',
    PROPNN = 'CH3C(=O)CH2ONO2',
    PYAC = 'CH3COCOOH',
    PMNN = 'C4H5NO5NO2',
    LISOPOH = 'IGNORE',
    DIBOO = 'C5H8OHO2',
    ISN1 = 'C5H8NO6',#http://mcm.leeds.ac.uk/MCM/browse.htt?species=NISOPNO3
    ISNOHOO = 'C5H8NO7',#http://mcm.leeds.ac.uk/MCM/browse.htt?species=C510O2
    ISNOOB = 'C5H7NO6NO3',
    ISNOOA = 'C5H6NO6',#http://mcm.leeds.ac.uk/MCM/browse.htt?species=NC4CO3
    )


check = re.compile('^((?:[A-Z][a-z]?\\d?)+) = (\S+);', flags = re.MULTILINE)
repl = re.compile('(?:([A-Z][a-z]?)(\d?))')
def atoms(matcho):
     count = matcho.groups()[1]
     if count == '':
        count = '1'
     atom = matcho.groups()[0]
     #return count * (atom + ' + ')
     return '%s%s + ' % (count, atom)

def formula(matcho):
    first = matcho.groups()[0]
    defn = matcho.groups()[1]
    if first.strip() in definitions:
        defn = definitions[first]
    elif defn == 'IGNORE':
        defn = first
        print('Guessing structure from name: ' + first, file = sys.stderr)
    else:
        defn = defn
    
    if defn != 'IGNORE':
        defn = defn.replace('(', '').replace(')', '').replace('=', '').replace('*','').replace('-', '')
        defn = re.sub(repl, atoms, defn)
        if defn[-3:] == ' + ':
            defn = defn[:-3]
    
    out = first + ' = ' + defn + ';'
    return out



def counteri():
    i = 0
    while True:
        i += 1
        if i == 133:
            i+=1
        yield i
counter = counteri()

foundord = []
def counter(ord):
    global foundord
    ord = str(ord).strip()
    foundord.append(ord)
    nord = foundord.count(ord)
    ord = ord + 'ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz'[nord-1]
    return ord

scinot = r'[-+]?[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?'
change_ord = re.compile(r'^(?P<type>[AD])\s+(?P<ORD>\d+)', re.M)
geos_spec_pattern = r'^(?P<active>[AD])\s+(?P<ORD>\d+[A-Za-z]?)\s+(?P<AR0>' + scinot + r')\s+(?P<BR0>' + scinot + r')\s+(?P<CR0>' + scinot + r') [012] (?P<type>(?:[ A-Z]|HR|QB)?)\s+(?P<FC>' + scinot + r')\s+(?P<FCT1>' + scinot + r')\.\s+(?P<FCT2>' + scinot + r')\.(?P<COMMENT>[ A-Za-z]{0,20})?[ ]*h?(?P<BRKEY>\d)?[ ]*(\n\s+(?P<AR1>' + scinot + r')\s+(?P<BR1>' + scinot + r')\s+(?P<CR1>' + scinot + r')\s+0( [A-Z])?\s+0\.00\s+0\.\s+0\.\s+(\n\s+(?P<AR2>' + scinot + r')\s+(?P<BR2>' + scinot + r')\s+(?P<CR2>' + scinot + r')\s+0( [A-Z])?\s+0\.00\s+0\.\s+0\.\s+)?)?(?=\n )'
geos_spec = re.compile(geos_spec_pattern, re.M)
M_RCT = re.compile(' (\+ )?M =')
geos_e_pattern = r'(?<=: )(?P<static>(?P<prate>GEOS_P\(%(scinot)s, %(scinot)s, %(scinot)s, %(scinot)s, %(scinot)s, %(scinot)s, %(scinot)s, %(scinot)s, %(scinot)s\));\n)(?P<key>{\s*(-)?\d+[A-Za-z]?}) (?P<rate>GEOS_E\(%(scinot)s, %(scinot)s, %(scinot)s, )lastrate\)\n\s+(?P<rxn>[^\n]+)\n' % globals()
geos_e_reorder = re.compile(geos_e_pattern, re.M)
rxn_continued = re.compile(r'\n(?=[+=])', re.M)
superfluous = re.compile(r'\+[ ]*([\n+=])', re.M)
geos_std_reorder_pattern = r'(?P<key>{\s*(-)?\d+[A-Za-z]?})[ ]+(?P<rate>(0\.00E\+00|TUV_J\(' + scinot + ', THETA\)|GEOS_Q\(TUV_J\(' + scinot + ', THETA\)\)|GEOS_[A-Z]{1,3}\((' + scinot + ', ){2,11}' + scinot + '\)))[ ]*\n\s+(?P<rxn>[^\n]+)\n'
geos_std_reorder = re.compile(geos_std_reorder_pattern, re.M)
geos_std_rate = re.compile(r'(?P<first>GEOS_STD\(' + scinot + ', ' + scinot + ', )(?P<CR>' + scinot + ')\)')
clean_up_trail = re.compile(r'([^{} \d\n]\S{1,8})[ ]{2,100}(?=\S)')
clean_up_head = re.compile(r'\+\s+')
stoic_sep = re.compile(r'(?<=[=+])(?P<stoic>' + scinot + r')(?=[^,\d )])')
plus_sep = re.compile(r'(\D[A-Z1-9]{1,8})[ ]+(?P<sign>[+-=])(?=\d)')
dry_dep = re.compile(r'[^\n]+\n[^\n]+\n=[^\n]+DRYDEP[^\n]+\n[^\n]+\n[^\n]+\n[^\n]+\n', re.M)
emission = re.compile(r'^A[^\n]+\n\s+EMISSION\s+\+\s+\n[^\n]+\n[^\n]+\n[^\n]+\n[^\n]+\n', re.M)
kinetic_prelude = re.compile(r'^[\n\s\S]+#=============================================================================\n# Kinetic reactions\n#=============================================================================\n\nBEGIN\n', re.M)
phot_prelude = re.compile(r'  9999 0\.00E-00   0\.0         0 0     0\.00     0\.     0\.     \n      END KINETIC', re.M)
section_comment = re.compile(r'(?:^#[=]+|^# \S+ reactions|^BEGIN)\s*\n', re.M)
postlude = re.compile(r'  9999 0\.00E-00 0\.0 0 0     0\.00 0\.     0\.[ ]+\n      END PHOTOLYSIS\n', re.M)
newscino2 = re.compile(r'(?<=\.\d{2})E')
newscino1 = re.compile(r'(?<=\.\d{1})E')
comment_spc_and_tail = re.compile(r'(\d.\d\d\d \b(O2|H2O|M|H2|CO2|H2O)\b [+{])')
comment_spc_and_head = re.compile(r'([+}] \d.\d\d\d \b(O2|H2O|M|H2|CO2|H2O)\b)')
double_comment_in = re.compile(r'{}')
double_comment_out = re.compile(r'} {')
template_dict = {' ': '{%(ORD)4s} GEOS_STD(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'A': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e)',
                 'B': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e)',
                 'C': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'D': '{%(ORD)4s} 0.00E+00',
                 'E': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, lastrate)',
                 'F': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'G': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e)',
                 'HR': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e)',
                 'K': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'L': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'N': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'O': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'P': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e, %(FC).6e, %(FCT1).6e, %(FCT2).6e)',
                 'Q': '{%(ORD)4s} GEOS_%(type)s(TUV_J(%(AR0).6e, THETA))',
                 'R': '{%(ORD)4s} TUV_J(%(AR0).6e, THETA)',
                 'T': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'V': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e)',
                 'W': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'X': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e, %(AR2).6e, %(BR2).6e, %(CR2).6e)',
                 'Y': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e)',
                 'Z': '{%(ORD)4s} GEOS_%(type)s(%(AR0).6e, %(BR0).6e, %(CR0).6e, %(AR1).6e, %(BR1).6e, %(CR1).6e, 1.400000e-21, 0.000000e+00, 2.200000e+02)',
                 'QB': '{%(ORD)4s} TUV_J(%(AR0).6e, THETA)',
                 }
def subrxn(matcho):
    global template_dict
    found = {}
    found.update(matcho.groupdict())
        
    typ = found['type']
    for k, v in found.iteritems():
        if k == 'COMMENT':
            found[k] = v
        elif k[:1] in 'ABCF' and v is not None:
            found[k] = float(v)
        if k == 'ORD':
            found[k] = str(v)

    if found['active'] == 'D':
        typ = found['type'] = 'D'
        
            
    template = template_dict[typ]
     
    return template % found

def print_usage():
    print("""

  Usage: %s inputpath

    Converts GEOS-Chem formatted reactions to KPP formatted
    reactions. Output is sent to stdout. Photolysis reactions
    must be updated manually after the conversion.
    
    Example:
      %s globchem.dat > geoschem.kpp
    
    """ % (sys.argv[0], sys.argv[0]))
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print_usage()
        exit()
    try:
        kpp_text_all = file(sys.argv[1], 'r').read()
        out_texts = []
        kinetic_text, photo_text = [si.strip() + '\n' for si in phot_prelude.split(kpp_text_all)]
        #kendpos = kpp_text_all.find('END KINETIC') + 12
        #kinetic_text = kpp_text_all[:kendpos]
        #photo_text = kpp_text_all[kendpos:]
        for kpp_text in [kinetic_text, photo_text]:
            kpp_text = emission.sub(r'', kpp_text)
            kpp_text = dry_dep.sub(r'', kpp_text)
            kpp_text = change_ord.sub(lambda match: '%s %s' % (match.groupdict()['type'], counter(match.groupdict()['ORD'])), kpp_text)
            kpp_text = rxn_continued.sub('', kpp_text)
            kpp_text = kinetic_prelude.sub(r'', kpp_text)
            kpp_text = section_comment.sub(r'', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = superfluous.sub(r'\1', kpp_text)
            kpp_text = geos_spec.sub(subrxn, kpp_text)
            kpp_text = geos_std_reorder.sub(r'\g<key>  \g<rxn> : \g<rate>;\n', kpp_text)
            kpp_text = geos_e_reorder.sub(r'\g<static>\g<key>  \g<rxn>: \g<rate>\g<prate>);\n', kpp_text)
            kpp_text = clean_up_trail.sub(r'\1 ', kpp_text)
            kpp_text = clean_up_head.sub(r'+ ', kpp_text)
            kpp_text = stoic_sep.sub(r'\g<stoic> ', kpp_text)
            kpp_text = plus_sep.sub(r'\1 \g<sign> ', kpp_text)
            kpp_text = phot_prelude.sub(r'', kpp_text)
            kpp_text = postlude.sub(r'', kpp_text)
            kpp_text = comment_spc_and_head.sub(r'{\1}', kpp_text)
            kpp_text = comment_spc_and_tail.sub(r'{\1}', kpp_text)
            kpp_text = double_comment_in.sub(r'', kpp_text)
            kpp_text = double_comment_out.sub(r' ', kpp_text)
            kpp_text = M_RCT.sub(r' =', kpp_text)
            template_dict[' '] = '{%(ORD)4s} TUV_J(%(AR0).6e, THETA)'
            out_texts.append(kpp_text)
        out_text = ''.join(out_texts).replace('\n\n', '\n').replace('\n\n', '\n').replace('\n\n', '\n').replace('\n\n', '\n').replace('\n\n', '\n').replace('\n\n', '\n')
        spc_text = re.sub(r'{[^}]+}\s*(.+?)\s*:.*', r'\1', out_text, flags = re.M)
        spc_text = re.sub(r'\s+' + scinot + '\s+', '', spc_text, flags = re.M)
        spc_text = re.sub(r'\d\.\d+', '', spc_text, flags = re.M)
        spc_text = re.sub(r'{[^}]+}', '', spc_text, flags = re.M)
        spc_text = re.sub(r'[ ]*', '', spc_text, flags = re.M).strip()
        spc_text = re.sub(r'(\s*[+=]\s*|\n+)', ',', spc_text, flags = re.M)
        spc_text = re.sub(r',+', ',', spc_text, flags = re.M)
        spc_text = ' = IGNORE;\n'.join(list(set(spc_text.split(',')))) + ' = IGNORE;\n'
        print('#DEFVAR')
        print(re.sub(check, formula, spc_text))
        print('#EQUATIONS')
        print(out_text)
    except IOError, e:
        print(str(e))
        print()
        print_usage()