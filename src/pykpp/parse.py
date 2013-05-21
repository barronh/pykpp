__all__ = '_parsefile _reactionstoic _stoicreaction _allspc'.split()
from pyparsing import *
import os
import re

spcname = Word(alphas, bodyChars = alphanums + '_')
real = Combine(Optional(Optional(oneOf("+ -")) + Word(nums)) +
               Optional("." + Word(nums)) +
               Optional(oneOf("e E") + Optional(oneOf("+ -")) +
               Word(nums))
              ).setName("real").setParseAction( lambda toks: '1.' if toks[0] == '' else toks[0] )

stoic = Group(Optional(real, default = 1.) + spcname).setResultsName('stoic')

inlinecomment = Suppress('{' + Word(alphanums + '*+-/%.;)([] ').setResultsName('rate') + '}')
lbl = Optional(Suppress('<') + Suppress(ZeroOrMore(' ')) + Regex('[^>]+').setResultsName("label") + Suppress('>'))
rcts = Group(delimitedList(Optional(inlinecomment) + stoic + Optional(inlinecomment), '+')).setResultsName("reactants")
prods = Group(delimitedList(Optional(inlinecomment) + stoic + Optional(inlinecomment), Regex('[+]'))).setResultsName("products")
rate = Regex('[^;]+').setResultsName('rate')

linecomment = Suppress('//' + Regex('.*?$'))
reaction = Group(Suppress(Optional(White())) + lbl + rcts + Suppress('=') + prods + Suppress(':') + rate + Suppress(';' + Optional(White()) + Optional(inlinecomment) + Optional(White())))

def reactions_parse(s, loc, toks):
    try:
        out = (ZeroOrMore(inlinecomment) + OneOrMore(reaction) + ZeroOrMore(inlinecomment)).parseString(toks[0][0], parseAll = True)
    except Exception, e:
        raise ParseFatalException('Error in parsing EQUATIONS (reactions start on character %d; lines numbered within): ' % loc + str(e))
    return out

reactions = Group(Suppress(lineStart + '#EQUATIONS') + Regex('[^#]+', flags = 16 + 8) + FollowedBy(Or([LineStart() + '#', stringEnd]))).setResultsName('EQUATIONS').addParseAction(reactions_parse)

#reactions = Group(Suppress(lineStart + '#EQUATIONS') + ZeroOrMore(inlinecomment) + OneOrMore(reaction) + ZeroOrMore(inlinecomment) + Optional(FollowedBy(lineStart + '#'))).setResultsName('reactions')


assignment = Group(Optional(inlinecomment) + spcname + Optional(' ') + Suppress('=') + Regex('[^;]+', flags = 8) + Suppress(';') + Optional(inlinecomment))

def initvalues_parse(loc, toks):
    try:
        return OneOrMore(Or([inlinecomment, assignment])).parseString(toks[0][0], parseAll = True)
    except Exception, e:
        raise ParseFatalException('Error in parsing INITVALUES (reactions start on character %d; lines numbered within): ' % loc + str(e))

initvalues = Group(Suppress(lineStart + '#INITVALUES') + Regex('[^#]+', flags = 16 + 8) + FollowedBy(Or([lineStart + '#', stringEnd]))).setResultsName('INITVALUES').addParseAction(initvalues_parse)

worldupdater = Optional(Group(Suppress(lineStart + '#WORLDUPDATER') + Regex('[^#]*') + Optional(FollowedBy(lineStart + '#')))).setResultsName('WORLDUPDATER')

codeseg = Group(Suppress(lineStart + "#INLINE ") + Regex('(F90|F77|C|MATLAB)') + '_' + Regex('(INIT|GLOBAL|RCONST|RATES|UTIL)') + Regex('[^#]+', flags = 16 + 8) + Suppress('#ENDINLINE')).setResultsName('CODESEG')

lookat = Optional(Group(Suppress(lineStart + '#LOOKAT') + Regex('.+') + lineEnd).setResultsName('LOOKAT'))
monitor = Optional(Group(Suppress(lineStart + '#MONITOR') + Regex('.+') + lineEnd).setResultsName('MONITOR'))
check = Optional(Group(Suppress(lineStart + '#CHECK') + Regex('.+') + lineEnd).setResultsName('CHECK'))

atoms = Optional(Group(Suppress(lineStart + '#ATOMS') + Regex('[^#]+', flags = re.MULTILINE + re.DOTALL) + FollowedBy(Or([lineStart + '#', stringEnd]))).setResultsName('ATOMS'))

defvar = Optional(Group(Suppress(lineStart + '#DEFVAR') + Regex('[^#]+', flags = re.MULTILINE) + FollowedBy(Or([lineStart + '#', stringEnd]))).setResultsName('DEFVAR'))

deffix = Optional(Group(Suppress(lineStart + '#DEFFIX') + Regex('[^#]+', flags = re.MULTILINE) + FollowedBy(Or([lineStart + '#', stringEnd]))).setResultsName('DEFFIX'))

language = Optional(Suppress(Regex('^#LANGUAGE.+$')))
integrator = Optional(Suppress(Regex('^#INTEGRATOR.+$')))
driver = Optional(Suppress(Regex('^#DRIVER.+$')))

parser = Each([initvalues, atoms, deffix, defvar, reactions, lookat, monitor, check, ZeroOrMore(codeseg), ZeroOrMore(linecomment)])

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
    includepaths.append(os.path.dirname(path))
    old = ''
    deftext = file(path).read()
    while old != deftext:
        old = deftext
        deftext = reinclude.sub(includeit, deftext)
    
    deftext = lcomment.sub('', deftext)
    return parser.parseString(deftext)

def _allspc(parsed):
    spc = set()
    #for k, v in parsed['initvalues']:
    #    spc.add(k)
    
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
