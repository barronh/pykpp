__all__ = '_parsefile _reactionstoic _stoicreaction _allspc'.split()
from pyparsing import *

spcname = Word(alphas, bodyChars = alphanums)
real = Combine(Optional(Optional(oneOf("+ -")) + Word(nums)) +
               Optional("." + Word(nums)) +
               Optional(oneOf("e E") + Optional(oneOf("+ -")) +
               Word(nums))
              ).setName("real").setParseAction( lambda toks: '1.' if toks[0] == '' else toks[0] )

stoic = Group(Optional(real, default = 1.) + spcname).setResultsName('stoic')

lbl = Optional(Suppress('<') + Suppress(ZeroOrMore(' ')) + Regex('[^>]+').setResultsName("label") + Suppress('>'))
rcts = Group(delimitedList(stoic, '+')).setResultsName("reactants")
prods = Group(delimitedList(stoic, Regex('[+-]'))).setResultsName("products")
rate = Regex('[^;]+').setResultsName('rate')
comment = Suppress('{' + Word(alphanums + '*+-/% ').setResultsName('rate') + '}')

reaction = Group(Suppress(Optional(White())) + lbl + rcts + Suppress('=') + prods + Suppress(':') + rate + Suppress(';' + Optional(White())))

strrct = '<R1> 2A + B = 3e4 C : SUN * T**3;\n<R2> 2C + E = 3e4 D : A * exp(300/T);'

reactions = Group(Suppress('#EQUATIONS') + ZeroOrMore(comment) + OneOrMore(reaction) + Optional(FollowedBy(lineStart + '#'))).setResultsName('reactions')


assignment = Group(spcname + Optional(' ') + Suppress('=') + Regex('[^;]+') + Suppress(';'))

initvalues = Group(Suppress('#INITVALUES') + OneOrMore(Or([comment, assignment])) + Optional(FollowedBy(lineStart + '#'))).setResultsName('initvalues')

worldupdater = Group(Suppress('#WORLDUPDATER') + Regex('[^#]*') + Optional(FollowedBy(lineStart + '#'))).setResultsName('worldupdater')

parser = Each([initvalues, reactions, Optional(worldupdater)])

def _parsefile(path):
    return parser.parseFile(path)

def _allspc(parsed):
    spc = set()
    #for k, v in parsed['initvalues']:
    #    spc.add(k)
    
    for rct in parsed['reactions']:
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
