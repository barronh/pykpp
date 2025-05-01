from __future__ import print_function
__all__ = ['mech', 'parse', 'plot', 'stdfuncs', 'main', 'funcs']
import warnings
import sys
from . import mech
from . import parse
from . import plot
from . import funcs
from . import stdfuncs
from . import main

warn = warnings.warn


def clean_showwarning(
    message, category, filename, lineno, file=None, line=None
):
    print(
        '**PYKPP:%s:%s:%s:\n  %s'
        % ((filename), lineno, category.__name__, message), file=sys.stderr
    )
    return


warnings.showwarning = clean_showwarning
