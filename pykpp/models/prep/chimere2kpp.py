#!/usr/bin/env python
from __future__ import print_function

import re
import sys

def print_usage():
    print("""

  Usage: %s inputpath

    Converts CHIMERE formatted reactions to KPP formatted
    reactions. Output is sent to stdout. Photolysis reactions
    must be updated manually after the conversion.
    
    Example:
      %s REACTIONS.univ.melchior1 > REACTIONS.univ.melchior1.kpp
    
    """ % (sys.argv[0], sys.argv[0]))
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print_usage()

    try:
        input = file(sys.argv[1], 'r').read()
        input = input.replace('->', ' = ')
        stoic_delim = re.compile(r'([\d\.])\*([A-Za-z])')
        spc_delim = re.compile(r'([A-Z0-9])\+([0-9A-Za-z])')
        output = stoic_delim.sub(r'\1 \2', input)
        output = spc_delim.sub(r'\1 + \2', output)

        comments = re.compile(r'^#.*$', re.M)
        output = comments.sub('', output)
        output = "#EQUATIONS\n" + output
        blank = re.compile(r'\n\s*\n', re.M)
        output = blank.sub('\n', output)

        def addlinenumber(matchobj):
            global output
            lineno = output.count("\n", 0, matchobj.start())+1
            return "\n<%d> " % lineno
    
        newline = re.compile(r'\n(?=[a-zA-Z0-9])', re.M)
        output = newline.sub(addlinenumber, output)
        std = re.compile(r'\s+(k\(T\)=Aexp\(-B/T\),A=([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),B=([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?))\s*$', re.M)
        output = std.sub(r' : (\2) * exp(-(\3) / TEMP) {\1};', output)

        k = re.compile(r'\s+(k=([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?))\s*$', re.M)
        output = k.sub(r' : (\2) {\1};', output)

        ks = re.compile(r'\s+(ks=([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?) depo\(.*\))\s*$', re.M)
        output = ks.sub(r' : 0. {\1};', output)


        kovert = re.compile(r'\s+(k\(1/T\)=A/T,A=([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?))$', re.M)
        output = kovert.sub(r' : \2 / TEMP {\1};', output)
        stdj = re.compile(r'\s+(J\(Z\)=photorate\((.+)\))\s*$', re.M)
        output = stdj.sub(r' : TUV_J(\2, THETA) {\1};', output)

        o3j = re.compile(r'\s+(J\(T,Z,H2O\)=photorate\((.+)\))\s*$', re.M)
        output = o3j.sub(r' : GEOS_JO3(TUV_J(\2, THETA)) {\1};', output)

        mtroe = re.compile(r'\s+(k\(T,M\)=mtroe\(([-+]{0,1}\d+(?:\.\d+(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:\.\d+(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:\.\d+(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:\.\d+(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:\.\d+(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:\.\d+(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:\.\d+(?:[eE][-+]{0,1}\d+)?)?)\))$', re.M)
        output = mtroe.sub(r' : CHIMERE_MTROE(\2, \3, \4, \5, \6, \7, \8) {\1};', output)

        troe = re.compile(r'\s+(k\(T,M\)=troe\(([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?)\))$', re.M)
        output = troe.sub(r' : CHIMERE_TROE(\2, \3, \4, \5, \6, \7, \8) {\1};', output)
        over = re.compile(r'\s+(k=([-+]{0,1}\d+(?:(?:\.\d+)?(?:[eE][-+]{0,1}\d+)?)?))$', re.M)
        output = over.sub(r' : \2 {\1};', output)

        param3 = re.compile(r'\s+(k\(T\)=Aexp\(-B/T\)\(300/T\)\*\*N,A=([-+]{0,1}\d+(?:(?:\.\d+)?(?:[eE][-+]{0,1}\d+)?)?),B=([-+]{0,1}\d+(?:(?:\.\d+)?(?:[eE][-+]{0,1}\d+)?)?),N=([-+]{0,1}\d+(?:(?:\.\d+)?(?:[eE][-+]{0,1}\d+)?)?))$', re.M)
        output = param3.sub(r': (\2) * exp(-(\3) / TEMP) * (300. / TEMP)**(\4) {\1};', output)

        special1 = re.compile(r'(k\(T\)=SPECIAL_1\(([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?)\))\s*$', re.M)
        output = special1.sub(r': CHIMERE_SPECIAL_4(\2, \3, \4, \5) {\1};', output)

        special2 = re.compile(r'(k\(T\)=SPECIAL_2\(([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?)\))\s*$', re.M)
        output = special2.sub(r': CHIMERE_SPECIAL_4(\2, \3, \4, \5) {\1};', output)

        special3 = re.compile(r'(k\(T\)=SPECIAL_3\(([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?)\))\s*$', re.M)
        output = special3.sub(r': CHIMERE_SPECIAL_3(\2, \3, \4, \5, \6, \7, \8, \9) {\1};', output)

        special4 = re.compile(r'(k\(T\)=SPECIAL_4\(([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?),([-+]{0,1}\d+(?:(?:\.)?(?:\d+)?(?:[eE][-+]{0,1}\d+)?)?)\))\s*$', re.M)
        output = special4.sub(r': CHIMERE_SPECIAL_4(\2, \3, \4, \5, \6, \7, \8, \9) {\1};', output)

        print(output)
    except Exception, e:
        print(e)
        print()
        print_usage()