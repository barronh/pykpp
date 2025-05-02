from __future__ import print_function
try:
    from setuptools import setup
except ImportError:
    from distutils.core import setup

version = '0.0.0'
with open('pykpp/__init__.py', 'r') as initf:
   while True:
       _l = initf.readline()
       if _l.startswith('__version__'):
           version = eval(_l.split('=')[1].strip())
           break

setup(
    name='pykpp',
    version=version,
    author='Barron Henderson',
    author_email='barronh@gmail.com',
    maintainer='Barron Henderson',
    maintainer_email='barronh@gmail.com',
    url='https://github.com/barronh/pykpp/',
    download_url='https://github.com/barronh/pykpp/archive/main.zip',
    long_description=(
        "pykpp is a KPP-like chemical mechanism parser that produces a box"
        + " model solvable by SciPy's odeint solver"
    ),
    packages=[
        'pykpp', 'pykpp.tuv', 'pykpp.morpho', 'pykpp.models', 'pykpp.tests',
        'pykpp.funcs'
    ],
    package_data={'pykpp.models': ['*.eqn', '*.txt', '*.kpp', '*.def']},
    scripts=['scripts/pykpprun.py'],
    entry_points={
      'console_scripts': [
          'pykpp = pykpp.__main__:main'
      ]
    },
    install_requires=['numpy', 'scipy', 'pyparsing', 'matplotlib', 'pandas']
)
