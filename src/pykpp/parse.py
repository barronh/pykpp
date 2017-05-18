from __future__ import print_function, unicode_literals
__all__ = '_parsefile _reactionstoic _stoicreaction _allspc'.split()
from pyparsing import *
import os
import sys
import re
from glob import glob
from warnings import warn    
trailing_zeros = re.compile('0{2,100}$')
RegexM = lambda expr: Regex(expr, flags = re.MULTILINE + re.I)
RegexMA = lambda expr: Regex(expr, flags = re.MULTILINE + re.DOTALL + re.I)
spcname = Word(alphas, bodyChars = alphanums + '_')
nextorend = FollowedBy(Or(lineStart + '#', stringEnd))
lcomment = re.compile('^//.*', re.M)
ilcomment = re.compile('{[^}]*}', re.M)

real = Combine(Optional(oneOf("+ -")) +
               Word(nums) + 
               Optional("." + Optional(Word(nums))) + 
               Optional(oneOf("e E d D") + Optional(oneOf("+ -")) +
               Word(nums))
              ).setName("real").setParseAction( lambda toks: '1.' if toks[0] == '' else toks[0].replace('D', 'E').replace('d', 'e'))

inlinecomment = Suppress('{' + RegexM('[^}]*').setResultsName('inline') + '}')
linecomment = Suppress('//' + RegexM('.*?$')).setResultsName('linecomment')
def cleanreal(toks):
    return inlinecomment.transformString(real.transformString(toks[0]))

kpphome = os.environ.get('KPP_HOME', '.')
includepaths = ['.', os.path.join(os.path.dirname(__file__), 'models')] + [kpphome] + glob(os.path.join(kpphome, '*'))

def includeit(matchobj):
    fname =     matchobj.groups()[0].strip()
    for dirtxt in includepaths:
        ipath = os.path.join(dirtxt, fname)
        if os.path.exists(ipath):
            print('Included', ipath)
            return '\n'.join(['', open(ipath, 'r').read(), ''])
    else:
        raise IOError('Unable to find %s; looked in (%s)' % (fname, ', '.join(includepaths)))

def includemodel(matchobj):
    mname = matchobj.groups()[0].strip()
    return '#include %s.def' % (mname,)

def reactions_parse(s, loc, toks):
    try:
        out = (ZeroOrMore(inlinecomment) + OneOrMore(reaction) + ZeroOrMore(inlinecomment)).parseString(toks[0][0], parseAll = True)
    except Exception as e:
        print(e.markInputline())
        raise ParseFatalException('Error in parsing EQUATIONS (reactions start on character %d; lines numbered within): ' % loc + str(e))
    return out

def initvalues_parse(loc, toks):
    try:
        return [linecomment.addParseAction(lambda : '').transformString(inlinecomment.addParseAction(lambda : '').transformString(real.transformString('\n'.join([l_.strip() for l_ in toks[0][0].split('\n')]))))]
    except Exception as e:
        raise ParseFatalException('Error in parsing INITVALUES (reactions start on character %d; lines numbered within): ' % loc + str(e))

def code_func(loc, toks):
    if toks[0][0] == 'PY_INIT':
        return ParseResults(toks[0][1], name = 'INIT')
    if toks[0][0] == 'PY_RCONST':
        return ParseResults(toks[0][1], name = 'RCONST')
    if toks[0][0] == 'PY_UTIL':
        return ParseResults(toks[0][1], name = 'UTIL')
    if toks[0][0] == 'PY_RATES':
        return ParseResults(toks[0][1], name = 'UTIL')
    
def ignoring(toks):
    warn('Ignoring' + toks[0][0])

def MakeOpt(pattern, name):
    return Optional(Group(pattern).setResultsName(name).addParseAction(ignoring))

stoic = Group(Optional(Combine(Optional(oneOf('+ -') + Optional(White().setParseAction(lambda toks: toks[0].replace('\n', ''))), default = '') + Optional(real, default = '1.')).setParseAction(lambda toks: trailing_zeros.sub('0', '1.' if toks[0] == '+ 1.' else toks[0]))) + spcname).setResultsName('stoic')

rate = RegexM('[^;]+').setResultsName('rate').addParseAction(cleanreal)
lbl = Optional(Suppress('<') + Suppress(ZeroOrMore(' ')) + Regex('[^>]+').setResultsName("label") + Suppress('>'))
#rcts = Group(delimitedList(Optional(inlinecomment) +stoic + Optional(inlinecomment), '+')).setResultsName("reactants")
rcts = Group(OneOrMore(stoic)).ignore(inlinecomment).setResultsName("reactants")
prods = Group(OneOrMore(stoic)).ignore(inlinecomment).setResultsName("products")
reaction = Group(Suppress(Optional(White())) + lbl + rcts + Suppress('=') + prods + Suppress(':') + rate + Suppress(';' + Optional(White()) + Optional(inlinecomment) + Optional(White())))
reactions = Group(Suppress(Literal('#EQUATIONS')) +
                  RegexMA('.+?(?=#|\Z)')
                 ).setResultsName('EQUATIONS').addParseAction(reactions_parse)

assignment = Group(Optional(inlinecomment) + spcname + Optional(' ') + Suppress('=') + RegexM('[^;]+') + Suppress(';') + Optional(inlinecomment))


def getdef(path):
    remodel = re.compile('^#MODEL (.+)', re.M + re.I)
    reinclude = re.compile('^#include (.+)', re.M + re.I)
    if isinstance(path, (str,)):
        includepaths.insert(0, os.path.dirname(path))
        if os.path.exists(path):
            deftext = open(path).read()
        else:
            for includepath in includepaths:
                possiblepath = os.path.join(includepath, path)
                if os.path.exists(possiblepath):
                    deftext = open(possiblepath).read()
                    break
            else:
                deftext = path
    elif hasattr(path, 'read'):
        deftext = path.read()
    
    old = ''
    while old != deftext:
        old = deftext
        deftext = remodel.sub(includemodel, deftext)
        deftext = reinclude.sub(includeit, deftext)
        deftext = lcomment.sub('', deftext)
        deftext = ilcomment.sub('', deftext)
    return deftext

def _parsefile(path):
    onelinekeys = 'MONITOR|LOOKAT|INTEGRATOR|LANGUAGE|STOICHMAT|HESSIAN|REORDER|DOUBLE|DRIVER|TRANSPORT|DUMMYINDEX|FUNCTION|CHECK'.split('|')
    sectionkeys = 'EQUATIONS|INITVALUES|ATOMS|DEFVAR|DEFFIX'.split('|')
    codesegs = re.compile('#INLINE\s+(\S+)(.+?)#ENDINLINE', re.M|re.DOTALL|re.UNICODE)
    sections = re.compile('#(' + '|'.join(sectionkeys) + ')(.+?)(?=#|\Z)', re.M|re.DOTALL|re.UNICODE)
    oneliner = re.compile('#(' + '|'.join(onelinekeys) + ')\s+(.+?)\s*;?\s*$', re.UNICODE | re.M)
    c = re.compile('{[^}]+}', re.M | re.UNICODE | re.DOTALL)
    
    deftext = getdef(path)
    
    retvals = {}

    for codename, code in codesegs.findall(deftext):
        if 'PY_' == codename[:3]:
            k = codename[3:]
            if k not in retvals:
                retvals[k] = code
            else:
                retvals[k] += '\n' + code

    for codename, code in sections.findall(deftext):
        # moved to getdef
        #oldcode = ''
        #while oldcode == code:
        #    oldcode = code
        #    code = c.sub('', code)
        
        if codename not in retvals:
            retvals[codename] = code
        else:
            retvals[codename] += '\n' + code

    for key, content in oneliner.findall(deftext):
        retvals[key] = content

    retvals = dict([(key, [value]) for key, value in retvals.items()])
    for key in onelinekeys:
        if key in retvals and key not in ('INTEGRATOR', 'MONITOR'):
            retvals[key] = retvals[key][0]
    retvals['EQUATIONS'] = reactions.parseString('#EQUATIONS\n' + retvals['EQUATIONS'][0])
    retvals.setdefault('INTEGRATOR', ['odeint'])
    return retvals

    
def _parsefile_old(path):
    """
    Moved pyparsing expressions inside function because they have memories...
    """
    # pyparsing and _parsefile_old_approach
    initvalues = Group(Suppress(RegexM('#INITVALUES')) + RegexMA('.+?(?=^#|\Z)')).setResultsName('INITVALUES' ).addParseAction(initvalues_parse)
    codeseg = Group(Suppress(RegexM("^#INLINE ")) + Regex('(PY|F90|F77|C|MATLAB)_(INIT|GLOBAL|RCONST|RATES|UTIL)') + RegexMA('[^#]+') + Suppress(RegexM('^#ENDINLINE.*'))).setResultsName('CODESEG').addParseAction(code_func)
    monitor = Optional(Group(Suppress(RegexM('^#MONITOR')) + RegexM('.+;') + Optional(inlinecomment)).setResultsName('MONITOR'))
    lookat = Optional(RegexM(r'^#LOOKAT.+').setResultsName('LOOKAT').addParseAction(lambda toks: toks[0].replace('#LOOKAT', '')))
    integrator = Optional(Suppress(RegexM('^#INTEGRATOR')) + RegexM('.+'), default = 'odeint').setResultsName('INTEGRATOR')
    defvar = MakeOpt(lineStart + RegexM('#DEFVAR') + RegexMA('.+') + FollowedBy(Or(lineStart + '#', stringEnd)), 'DEFVAR')
    deffix = MakeOpt(lineStart + RegexM('#DEFFIX') + RegexMA('.+') + FollowedBy(Or(lineStart + '#', stringEnd)), 'DEFFIX')
    check = MakeOpt(lineStart + '#CHECK' + restOfLine, 'CHECK')
    atoms = MakeOpt(lineStart + '#ATOMS' + RegexMA('.+?(?=^#)'), 'ATOMS')
    reorder = MakeOpt(lineStart + '#REORDER' + restOfLine, 'REORDER')
    double = MakeOpt(lineStart + '#DOUBLE' + restOfLine, 'DOUBLE')
    language = MakeOpt(lineStart + '#LANGUAGE' + restOfLine, 'LANGUAGE')
    driver = MakeOpt(lineStart + '#DRIVER' + restOfLine, 'DRIVER')
    hessian = MakeOpt(lineStart + '#HESSIAN' + restOfLine, 'HESSIAN')
    stoicmat = MakeOpt(lineStart + '#STOICMAT' + restOfLine, 'STOICHMAT')
    mex = MakeOpt(lineStart + '#MEX' + restOfLine, 'MEX')
    stochastic = MakeOpt(lineStart + '#STOCHASTIC' + restOfLine, 'STOCHASTIC')
    transportall = MakeOpt(lineStart + '#TRANSPORT' + restOfLine, 'TRANSPORT')
    dummyidx = MakeOpt(lineStart + '#DUMMYINDEX' + restOfLine, 'DUMMYINDEX')
    function = MakeOpt(lineStart + '#FUNCTION' + restOfLine, 'FUNCTION')

    elements = [language, Optional(initvalues),
                atoms, deffix, defvar, integrator,
                ZeroOrMore(reactions),
                lookat, monitor, check, driver, ZeroOrMore(codeseg),
                ZeroOrMore(linecomment), ZeroOrMore(inlinecomment),
                reorder, double, hessian, stoicmat, dummyidx, transportall, stochastic, mex, function]

    for i in elements:
        i.verbose_stacktrace = True

    parser = Each(elements)
    if isinstance(path, (str,)) and not os.path.exists(path):
        base, ext = os.path.splitext(path)
        if ext == 'kpp':
            path = os.path.join(kpphome, 'examples', path)
        elif ext == 'def':
            path = os.path.join(kpphome, 'models', path)
        elif ext == '':
            tpath = os.path.join(kpphome, 'examples', path + '.kpp')
            if not os.path.exists(path):
                tpath = os.path.join(kpphome, 'models', path + '.def')
            path = tpath
            del tpath
    
    deftext = getdef(path)
    try:
        del code_func.got_code
    except Exception as e:
        pass
    #open('test.txt', 'w').write(deftext)
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

def _prune_meta_species(rxns, *meta_species):
    remove = []
    for spci, (stc, spc) in enumerate(rxns['reactants']):
        if spc.strip() in meta_species:
            warn('Ignoring %s' % spc)
            remove.append(spci)
    [rxns['reactants'].pop(spci) for spci in remove[::-1]]
    remove = []
    for spci, (stc, spc) in enumerate(rxns['products']):
        if spc.strip() in meta_species:
            warn('Ignoring %s' % spc)
            remove.append(spci)
    [rxns['products'].pop(spci) for spci in remove[::-1]]

