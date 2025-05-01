#!/usr/bin/env python
from __future__ import print_function
import re
multispace = re.compile('[ ]{2,100}', re.M)
continuation = re.compile('\n&', re.M)
comment = re.compile('![^\n]*\n', re.M)
thermal = re.compile('\nTHERMAL A-FACT ([\d.E\-+]+) E/R ([\d.E\-+]+)', re.M)
special = re.compile('\nSPECIAL RWK\(\#\) = (.*?)(?=\n)', re.M)
thermalt2 = re.compile('\nTHERMAL-T2 A-FACT ([\d.E\-+]+) E/R ([\d.E\-+]+)', re.M)
photolysis = re.compile(r'\nPHOTOLYSIS', re.M)
stoic = re.compile(r'(?<=[=+}])\s*(?:([\d.]+)\s+)?([A-Z]\S+)[ ]*(?=[^:+*]+:)')
troe = re.compile(r'\nTROE KO ([\d.Ee\-]+) N ([\d.Ee\-]+) KINF ([\d.Ee\-]+) M ([\d.Ee\-]+)', re.M)
troeequil = re.compile(r'\nTROE-EQUIL KO ([\d.Ee\-]+) N ([\d.Ee\-]+) KINF ([\d.Ee\-]+) M ([\d.Ee\-]+)\n A-FACT ([\d.Ee\-]+) B ([\d.Ee\-]+)', re.M)
rxnn = re.compile('(^\d{3})\s', re.M)
constant_spcs = "O2 M N2 H2 H2O".split()

constant_spc_re = re.compile(r'((?:(?:\+\s*)?\b(?:' + '|'.join(constant_spcs) + r')\b(?:\s*\+\s*)?)+)', re.M)
constant_spc_rct = re.compile(r'\{((?:(?:\+\s*)?\b(?:' + '|'.join(constant_spcs) + r')\b(?:\s*\+\s*)?)+)\}.*?=.*?\n', re.M)


def stoic_repl(matchobj):
    stoic, spc = matchobj.groups()
    if stoic is None:
        return ' %s + ' % spc
    else:
        return ' %s %s + ' % (stoic, spc)

def sbox2kpp(text):
    text = rxnn.sub(r'<\1> ', text)
    text = comment.sub('', text)
    text = continuation.sub(r' ', text)
    text = re.sub(r'->', '= ', text)
    text = multispace.sub(' ', text)
    text = thermal.sub(r' : \1 * exp(-(\2) / TEMP );', text)
    text = thermalt2.sub(r' : \1 * TEMP**2 * exp(-(\2) / TEMP);', text)
    text = photolysis.sub(r' : TUV_J(, THETA);', text)
    text = troe.sub(r' : RACM_TROE(\1, \2, \3, \4);', text)
    text = troeequil.sub(r' : RACM_TROE_EQUIL(\1, \2, \3, \4, \5, \6);', text)
    text = special.sub(r' : \1;', text)
    new_text = text
    old_text = ''
    while old_text != new_text:
        old_text = new_text
        new_text = stoic.sub(stoic_repl, old_text)

    text = re.sub(r'\+\s+(?==)', '', new_text)
    text = re.sub(r'\+\s+(?=:)', '', text)
    text = re.sub(r'=\s+:', '= DUMMY :', text)
    text = re.sub(r'\n\n', '\n', text)
    text = re.sub(r'} 1.000', '} ', text)
    text = re.sub(r' \* exp\(-\(0\.0\) / TEMP \)', '', text)
    text = re.sub(r'\n\s*\n', r'\n', text)
    def scaleconstants(matcho):
        groups = [[vs.strip() for vs in v.split('+') if vs.strip() != ''] for v in matcho.groups() if v is not None]
        outtext = matcho.group()
    
        for group in groups:
            outtext = outtext.replace(':', ': ' + ' * '.join(group) + ' *')
        return outtext
    
    text = constant_spc_re.sub(r'{\1}', text)
    text = constant_spc_rct.sub(scaleconstants, text)
    text = text.replace('= :', '= DUMMY :')
    return text

def print_usage():
    print( """

  Usage: %s inputpath

    Converts SBOX (RACM2) formatted reactions to KPP formatted
    reactions. Output is sent to stdout. Photolysis reactions
    must be updated manually after the conversion.
    
    Example:
      %s RACM2_sbox > RACM2_sbox.kpp
    
    """ % (sys.argv[0], sys.argv[0]))
    
if __name__ == '__main__':
    import sys
    if len(sys.argv) != 2:
        print_usage()
        exit()
    try:
        text = file(sys.argv[1], 'r').read()
        text = sbox2kpp(text)
        print(text)
    except Exception, e:
        print_usage()
        raise e