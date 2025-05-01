from __future__ import print_function
import sys, re
from PseudoNetCDF.cmaqfiles.jtable import jtable

def print_usage():
    print("""

  Usage: %s inputpath

    Converts CMAQ formatted reactions to KPP formatted
    reactions. Output is sent to stdout. Photolysis reactions
    must be updated manually after the conversion.
    
    Example:
      %s mech.def > mech.kpp
    
    """ % (sys.argv[0], sys.argv[0]))
if __name__ == '__main__':
    if len(sys.argv) != 2:
        print_usage()
        exit()
    try:
        #j = jtable(path = '/Volumes/LaCie/JTABLE_1985172')
        constant_spcs = "O2 M N2 H2 H2O".split()

        constant_spc_re = re.compile(r'((?:(?:\+\s*)?\b(?:' + '|'.join(constant_spcs) + r')\b(?:\s*\+\s*)?)+)', re.M)
        constant_spc_rct = re.compile(r'\{((?:(?:\+\s*)?\b(?:' + '|'.join(constant_spcs) + r')\b(?:\s*\+\s*)?)+)\}.*?=.*?\n', re.M)

        newtxt = mechtxt = file(sys.argv[1], 'r').read()

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
        def fillwithzerocmaq10(template, matcho):
            groups = [g for g in matcho.groups()]
            if groups[-1] is None:
                groups[-1] = '1.0'
            if groups[-2] is None:
                groups[-2] = '0.6'
            groups = [g or '0.0' for g in groups]
            groups = map(lambda x: ('%9.2E' % float(x)).strip().replace('-0.00E+00', '0.00E+00'), [groups[0], groups[2], str(-float(groups[1])), groups[3], groups[5], str(-float(groups[4]))] + groups[6:])
            groups = [adddecimal.sub(r'\1.0', g) for g in groups]
            return template % tuple(groups)
            
        def fillwithzerocmaq1to4(template, matcho):
            groups = [i for i in matcho.groups()]
            groups = [groups[0], groups[2], groups[1]] + groups[3:]
            prefix = ''
            expr = ''
            suffix = ''
            
            if groups[0] is not None:
                prefix += 'k'
                if groups[1] is None and groups[2] is None:
                    expr += '%s' % groups[0]
                else:
                    expr += 'A'
                    suffix += ',A=' + ('%9.2E' % float(groups[0])).strip().replace('-0.00E+00', '0.00E+00')
            if groups[1] is not None:
                prefix += '(T)'
                expr += 'exp(-B/T)'
                suffix += ',B=' + ('%9.2E' % float(groups[1])).strip().replace('-0.00E+00', '0.00E+00')
            if groups[2] is not None:
                if '(T)' not in prefix:
                    prefix += '(T)'
                expr += '(300/T)**N'
                suffix += ',N=' + ('%9.2E' % (-float(groups[2]))).strip().replace('-0.00E+00', '0.00E+00')
            return ' : ' + prefix + '=' + expr +  suffix + '\n'
        def tuvj(matcho):
            groups = matcho.groups()
            if len(groups) > 2:
                raise ValueError, "Oops too many TUV groups"
            retval = 'J(Z)=photorate(%s)\n' % (groups[1])
    
            if groups[0] is not None and groups[0] != '1.0':
                retval = ('%s * ' % groups[0]) + retval
            retval = ' : ' + retval
            return retval        
    
        cmaq_1to4 = re.compile('(?<!%\d\s)#\s(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s([\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?)){1}\n', re.M)
        newtxt = cmaq_1to4.sub(lambda m: fillwithzerocmaq1to4('      k(T)=Aexp(-B/T)(300/T)**N,A=%s,B=%s,N=%s\n', m), newtxt)

        cmaq_8 = re.compile(r'%\d\s#\s(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))\n')
        newtxt = cmaq_8.sub(lambda m: fillwithzero(' : k(T,M)=%%2tro(%s, %s, %s, %s, %s, %s)\n', m), newtxt)

        cmaq_9 = re.compile(r'%\d\s#\s(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))\n')
        newtxt = cmaq_9.sub(lambda m: fillwithzero(' : k(T,M)=%%3tro(%s, %s, %s, %s)\n', m).replace('e', 'd').replace('E', 'd'), newtxt)

        cmaq_10 = re.compile(r'(?<!%\d\s)#\s(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s([\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s([\d.Ee\-+]*))?(?:\s@\s([\d.Ee\-+]*))?(?:\s[&;]\s?))(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s(?:[\d.Ee\-+]*))?(?:\s[&;]\s?))?(?:([\d.Ee\-+]*)(?!.*/)(?:\s\^\s(?:[\d.Ee\-+]*))?(?:\s@\s(?:[\d.Ee\-+]*))?(?:\s[&;]\s?))?\n')
        newtxt = cmaq_10.sub(lambda m: fillwithzerocmaq10(' : tro_cb05(%s, %s, %s, %s, %s, %s, %s, %s)\n', m), newtxt)

        tuv_j = re.compile(r'(?<!%\d\s)#\s([\d.Ee\-+]*)\s*/\s*<\s*(\S+)\s*>\s*;\s*\n', re.M)
        newtxt = tuv_j.sub(tuvj, newtxt)

        label = re.compile(r'^<(\s*\S*\s*)>', re.M)
        newtxt = label.sub(r'<\1>', newtxt)

        header = re.compile(r'.*REACTIONS\s*\[CM\] =\s*\n', re.M | re.S)
        newtxt = header.sub('', newtxt)

        footer = re.compile(r'\n(endmech|END).*', re.M | re.I | re.S)
        newtxt = footer.sub('', newtxt)

        newlabels = [re.compile('<(\s*\S+\s*)>').search(line).groups()[0].strip() for line in newtxt.split('\n') if line != '']

        def deref_rate(matcho):
            groups = matcho.groups()
    
            rate_num = newlabels.index(groups[1].strip()) + 1
            return ': %s * RCONST(%d) {%s};\n' % (groups[0], rate_num, groups[1])

        rateref = re.compile('(?<!%\d\s)#\s([\d.Ee\-+]*)\*K<(\S+)>\s*;\s*\n', re.M)
        newtxt = rateref.sub(deref_rate, newtxt)

        newtxt = '#EQUATIONS\n' + newtxt

        newtxt = re.compile('(?<![eE])([+])(?!\d)', re.M).sub(r' \1 ', newtxt)
        newtxt = re.compile('[ ]+', re.M).sub(' ', newtxt)

        def removeasterisks(matcho):
            groups = matcho.groups()
            if len(groups) != 1:
                raise ValueError, "Should just be 1"
    
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
        newtxt = newtxt.replace(' = ', '->')
        newtxt = newtxt.replace(' ', '')
        newtxt = newtxt.replace(':', '      ')
        newtxt = newtxt.replace('}', '')
        newtxt = newtxt.replace('{', '')
        newtxt = re.compile('^<.+?>', re.M).sub(' ', newtxt)
        print('Has CMAQ_1to4', 'CMAQ_1to4' in newtxt, file = sys.stderr)
        print('Has CMAQ_8', 'CMAQ_8' in newtxt, file = sys.stderr)
        print('Has CMAQ_9', 'CMAQ_9' in newtxt, file = sys.stderr)
        print('Has CMAQ_10', 'CMAQ_10' in newtxt, file = sys.stderr)
        print(newtxt, file = sys.stdout)
    except KeyError as e:
        print_usage()
        import pdb; pdb.set_trace()
        raise e