import sys
import re
intxt = file(sys.argv[1], 'r').read()

counter = 0
commentonly = re.compile('(=\s*{[^}]+}\s*):')
intertspcs = r'(?:\b(?:M|H2|CO2|H2O|hv)\b)'
stoic = r'(?:\d*\.?\d+\*)'
first = '(?:(?<=)\s*' + stoic + '?' + intertspcs + '\s*\+?)'
second = r'(?:\+[ ]*' + stoic + '?' + intertspcs + ')'
inert = re.compile('((?:(?:' + first + '|' + second + ')\s*)+)', flags = re.M)
def makedble(ra_in):
    ra_in = ra_in.strip()
    if 'd' in ra_in or 'D' in ra_in or ra_in == '':
        ra_out = ra_in
    elif 'e' in ra_in:
        ra_out = ra_in.replace('e', 'd')
    elif 'E' in ra_in:
        ra_out = ra_in.replace('E', 'D')
    elif '.' in ra_in:
        ra_out = ra_in + 'd0'
    else:
        ra_out = ra_in + '.d0'
    return ra_out
        
def rateform(matcho):
    global counter
    gd = matcho.groupdict()
    kwds = {}
    kwds.update(gd)
    if kwds['more'] is None:
        kwds['more'] = ''
    if kwds['rateargs'] is None:
        kwds['rateargs'] = ''
    rateargs = kwds['rateargs'].strip()
    rateargs = [makedble(ra) for ra in rateargs.split(',')]
    kwds['rateargs'] = rateargs
    label = gd['label']
    ratefunc = kwds['ratefunc'] = 'am3_' + {None: 'standard'}.get(label, label)
    nargs = len(rateargs)
    minarg = dict(am3_standard = 3).get(ratefunc, nargs)
    rateargs = rateargs + ['0.0d0'] * (minarg - nargs)
    kwds['rateargs'] = ', '.join(rateargs)
    counter += 1
    kwds['counter'] = counter
    reaction = (kwds['reaction'] + kwds['more']).replace('\t', ' ')
    reactants, products = reaction.split('->')
    reactants = inert.sub(r'{\1}', reactants)
    products = inert.sub(r'{\1}', products)
    reaction = reactants + '->' + products
    kwds['factor'] = ''
    kwds['reaction'] = reaction
    
    out = '%(reaction)s : %(factor)s%(ratefunc)s(%(rateargs)s);\n' % kwds
    out = out.replace('*', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ')
    if kwds['comment'] == '*':
        out = '{active?}' + out
    return '\n' + out


species_g = re.search(r'^\s*Explicit\s*?\n(?P<explicit>.+)^\s*End Explicit\s*Implicit\s*\n(?P<implicit>.+)^\s*End Implicit', intxt, re.M + re.DOTALL).groupdict()
species_txt = '#DEFVAR\n' + ' = IGNORE;\n'.join((species_g['explicit'] + '\n'+ species_g['implicit']).strip().replace('\n', ',').replace(' ', '').replace(',,', ',').replace(',,', ',').replace(',,', ',').split(',')) + ' = IGNORE;\n'
photolysis_txt = re.search(r'^\s*Photolysis\s*\n(?P<photolysis>.+)^\s*End Photolysis', intxt, re.M + re.DOTALL).groupdict()['photolysis']
kinetic_txt = re.search(r'^\s*Reactions\s*\n(?P<kinetic>.+)^\s*End Reactions', intxt, re.M + re.DOTALL).groupdict()['kinetic']

n = -1
while n != 0:
    photolysis_txt, n = re.subn(r'(?:^|\n)\s*(?:\[(\S+)\])?\s*([0-9a-zA-Z* >+-.\t]+)', r'\n{\1}: \2: TUV_J(\1, THETA, 1);', photolysis_txt, flags = re.M)
photolysis_txt = inert.sub(r'{\1}', photolysis_txt)

n = -1
while n != 0:
    kinetic_txt, n = re.subn(r'(?:^|\n)(?P<comment>\*)?\s*(?:\[(?P<label>\S+)\])?\s*(?P<reaction>[0-9a-zA-Z][0-9a-zA-Z* >+-.\t]+)(?:;\s*(?P<rateargs>[^\n]+))?\n(?:(?P<more>\s+[+][0-9a-zA-Z* >+-.\t]+)\n)*', rateform, kinetic_txt, re.M)

kinetic_txt = '\n'.join(['{%d} ' % (linei + 1) + line.strip() for linei, line in enumerate(kinetic_txt.strip().split('\n'))])

kinetic_txt = commentonly.sub(r'\1 DUMMY :', kinetic_txt)
photolysis_txt = commentonly.sub(r'\1 DUMMY :', photolysis_txt)
print species_txt
print '#EQUATIONS'
print photolysis_txt.strip().replace('->', '=').replace('*', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ')
print kinetic_txt.replace('->', '=')
