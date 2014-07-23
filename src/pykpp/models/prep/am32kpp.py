import sys
import re
intxt = file(sys.argv[1], 'r').read()

counter = 0
def rateform(matcho):
    global counter
    gd = matcho.groupdict()
    kwds = {}
    kwds.update(gd)
    if kwds['more'] is None:
        kwds['more'] = ''
    if kwds['rateargs'] is None:
        kwds['rateargs'] = ''
    kwds['rateargs'] = kwds['rateargs'].strip()
    label = gd['label']
    kwds['ratefunc'] = 'am3_' + {None: 'standard'}.get(label, label)
    counter += 1
    kwds['counter'] = counter
    out = '%(reaction)s %(more)s : %(ratefunc)s(%(rateargs)s);\n' % kwds
    out = out.replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ').replace('  ', ' ')
    if kwds['comment'] == '*':
        out = '{active?}' + out
    return '\n' + out.replace('*', ' ')

photolysis_txt = re.search(r'^\s*Photolysis\s*\n(?P<photolysis>.+)^\s*End Photolysis', intxt, re.M + re.DOTALL).groupdict()['photolysis']
kinetic_txt = re.search(r'^\s*Reactions\s*\n(?P<kinetic>.+)^\s*End Reactions', intxt, re.M + re.DOTALL).groupdict()['kinetic']

n = -1
while n != 0:
    photolysis_txt, n = re.subn(r'(?:^|\n)\s*(?:\[(\S+)\])?\s*([0-9a-zA-Z* >+-.\t]+)', r'\n{\1}: \2: TUV_J(\1, THETA, 1);', photolysis_txt, re.M)

n = -1
while n != 0:
    kinetic_txt, n = re.subn(r'(?:^|\n)(?P<comment>\*)?\s*(?:\[(?P<label>\S+)\])?\s*(?P<reaction>[0-9a-zA-Z][0-9a-zA-Z* >+-.\t]+)(?:;\s*(?P<rateargs>[^\n]+))?\n(?:(?P<more>\s+[+][0-9a-zA-Z* >+-.\t]+)\n)*', rateform, kinetic_txt, re.M)

kinetic_txt = '\n'.join(['{%d} ' % (linei + 1) + line.strip() for linei, line in enumerate(kinetic_txt.strip().split('\n'))])
print photolysis_txt.strip()
print kinetic_txt
