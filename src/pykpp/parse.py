__all__ = '_parsefile _reactionstoic _stoicreaction _allspc'.split()
from pyparsing import *
import os
import re
from warnings import warn

RegexM = lambda expr: Regex(expr, flags = re.MULTILINE + re.I)
RegexMA = lambda expr: Regex(expr, flags = re.MULTILINE + re.DOTALL + re.I)

spcname = Word(alphas, bodyChars = alphanums + '_')

real = Combine(Optional(oneOf("+ -")) + Word(nums) +
               Optional("." + Optional(Word(nums)) +
               Optional(oneOf("e E d D") + Optional(oneOf("+ -")) +
               Word(nums)))
              ).setName("real").setParseAction( lambda toks: '1.' if toks[0] == '' else toks[0].replace('D', 'E').replace('d', 'e'))
              
#stoic = Group(Optional(real, default = '1.') + spcname).setResultsName('stoic')


stoic = Group(Optional(Combine(Optional(oneOf('+ -') + Optional(White().setParseAction(lambda toks: toks[0].replace('\n', ''))), default = '') + Optional(real, default = '1.'))) + spcname).setResultsName('stoic')

inlinecomment = Suppress('{' + RegexM('[^}]*').setResultsName('inline') + '}')
lbl = Optional(Suppress('<') + Suppress(ZeroOrMore(' ')) + Regex('[^>]+').setResultsName("label") + Suppress('>'))
rcts = Group(delimitedList(Optional(inlinecomment) +stoic + Optional(inlinecomment), '+')).setResultsName("reactants")
rcts = Group(OneOrMore(Optional(inlinecomment) +stoic + Optional(inlinecomment))).setResultsName("reactants")
prods = Group(OneOrMore(Optional(inlinecomment) + stoic + Optional(inlinecomment))).setResultsName("products")
rate = RegexM('[^;]+').setResultsName('rate')

linecomment = Suppress('//' + RegexM('.*?$'))
reaction = Group(Suppress(Optional(White())) + lbl + rcts + Suppress('=') + prods + Suppress(':') + rate + Suppress(';' + Optional(White()) + Optional(inlinecomment) + Optional(White())))

def reactions_parse(s, loc, toks):
    try:
        out = (ZeroOrMore(inlinecomment) + OneOrMore(reaction) + ZeroOrMore(inlinecomment)).parseString(toks[0][0], parseAll = True)
    except Exception, e:
        import pdb; pdb.set_trace()
        raise ParseFatalException('Error in parsing EQUATIONS (reactions start on character %d; lines numbered within): ' % loc + str(e))
    return out

reactions = Group(Suppress(RegexM('^#EQUATIONS')) + RegexMA('.+?(?=^#)')).setResultsName('EQUATIONS').addParseAction(reactions_parse)

#reactions = Group(Suppress(lineStart + '#EQUATIONS') + ZeroOrMore(inlinecomment) + OneOrMore(reaction) + ZeroOrMore(inlinecomment) + Optional(FollowedBy(lineStart + '#'))).setResultsName('reactions')


assignment = Group(Optional(inlinecomment) + spcname + Optional(' ') + Suppress('=') + RegexM('[^;]+') + Suppress(';') + Optional(inlinecomment))

def initvalues_parse(loc, toks):
    try:
        return [linecomment.addParseAction(lambda : '').transformString(inlinecomment.addParseAction(lambda : '').transformString(real.transformString('\n'.join([l_.strip() for l_ in toks[0][0].split('\n')]))))]
    except Exception, e:
        raise ParseFatalException('Error in parsing INITVALUES (reactions start on character %d; lines numbered within): ' % loc + str(e))

initvalues = Group(Suppress(RegexM('^#INITVALUES')) + RegexMA('.+?(?=^#)')).setResultsName('INITVALUES' ).addParseAction(initvalues_parse)

worldupdater = Optional(Group(Suppress(lineStart + '#WORLDUPDATER') + Regex('[^#]+') + Optional(FollowedBy(lineStart + '#')))).setResultsName('WORLDUPDATER')

def code_func(loc, toks):
    if toks[0][0] == 'F90_INIT':
        return ParseResults([linecomment.addParseAction(lambda : '').transformString(inlinecomment.addParseAction(lambda : '').transformString(real.transformString('\n'.join([l_.strip() for l_ in toks[0][1].split('\n')]))))], name = 'INIT')

codeseg = Group(Suppress(RegexM("^#INLINE ")) + Regex('(F90|F77|C|MATLAB)_(INIT|GLOBAL|RCONST|RATES|UTIL)') + RegexMA('[^#]+') + Suppress(RegexM('^#ENDINLINE.*'))).setResultsName('CODESEG').addParseAction(code_func)

monitor = Optional(Group(Suppress(RegexM('^#MONITOR')) + RegexM('.+;') + Optional(inlinecomment)).setResultsName('MONITOR'))

def ignoring(toks):
    print 'Ignoring', toks[0][0]
lookat = Optional(Group(Or([Suppress(RegexM('^#LOOKAT')) + 'ALL' + Suppress(RegexM('.*')), Suppress(RegexM('^#LOOKAT')) + RegexM('.+;')])).setResultsName('LOOKAT'))

check = Optional(Group(RegexM('^#CHECK') + RegexM('.+')).setResultsName('CHECK').addParseAction(ignoring))

atoms = Optional(Group(RegexM('^#ATOMS') + RegexMA('.+?(?=^#)')).setResultsName('ATOMS').addParseAction(ignoring))

defvar = Optional(Group(RegexM('^#DEFVAR') + RegexMA('.+?(?=^#)')).setResultsName('DEFVAR').addParseAction(ignoring))

deffix = Optional(Group(RegexM('^#DEFFIX') + RegexMA('.+?(?=^#)')).setResultsName('DEFFIX').addParseAction(ignoring))

integrator = Optional(Group(RegexM('^#INTEGRATOR') + RegexM('.+')), default = 'lsoda').setResultsName('INTEGRATOR')

language = Optional(Group(RegexM('^#LANGUAGE') + RegexM('.+')).setResultsName('LANGUAGE').addParseAction(ignoring))

driver = Optional(Group(RegexM('^#DRIVER') + RegexM('.+')).addParseAction(ignoring))

elements = [language, Optional(initvalues), atoms, deffix, defvar, reactions, lookat, monitor, check, integrator, driver, ZeroOrMore(codeseg), ZeroOrMore(linecomment), ZeroOrMore(inlinecomment)]

for i in elements:
    i.verbose_stacktrace = True
parser = Each(elements)

includepaths = ['.']
def includeit(matchobj):
    fname = matchobj.groups()[0]
    for dirtxt in includepaths:
        ipath = os.path.join(dirtxt, fname)
        if os.path.exists(ipath):
            print ipath
            return file(ipath, 'r').read()
    else:
        raise IOError('Unable to find %s; looked in (%s)' % (fname, ', '.join(includepaths)))
    
def _parsefile(path):
    reinclude = re.compile('^#include (.+)', re.M)
    lcomment = re.compile('^//.*', re.M)
    ilcomment = re.compile('^{[^}]*}', re.M)
    includepaths.insert(0, os.path.dirname(path))
    old = ''
    deftext = file(path).read()
    while old != deftext:
        old = deftext
        deftext = reinclude.sub(includeit, deftext)
    
    deftext = lcomment.sub('', deftext)
    deftext = ilcomment.sub('', deftext)
    
    #file('test.txt', 'w').write(deftext)
    return parser.parseString(deftext)

def _allspc(parsed):
    spc = set()
    for rct in parsed['EQUATIONS']:
        for v, k in rct['reactants']:
            spc.add(k)
        for v, k in rct['products']:
            spc.add(k)
    spc = list(spc)
    spc.sort()
    return spc

def _stoicreaction(spc, reaction):
    stoic = 0.
    for role, fac in [('reactants', -1), ('products', 1)]:
        for v, k in reaction[role]:
            if spc == k:
                stoic += fac * eval(v)
    return stoic

def _reactionstoic(spc, reactions):
    stoics = {}
    for ri, reaction in enumerate(reactions):
        stoic = _stoicreaction(spc, reaction)
        if stoic != 0.:
            stoics[ri] = stoic
    return stoics
