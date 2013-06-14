from distutils.core import setup

setup(name = 'pykpp',
      version = '0.9',
      author = 'Barron Henderson',
      author_email = 'barronh@gmail.com',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@gmail.com',
      url='https://github.com/barronh/pykpp/',
      download_url='https://github.com/barronh/pykpp/archive/master.zip',
      long_description="pykpp is a KPP-like chemical mechanism parser that produces a box model solvable by SciPy's odeint solver",
      packages = ['pykpp', 'pykpp.tuv', 'pykpp.morpho', 'pykpp.models', 'pykpp.test'],
      package_dir = {'pykpp': 'src/pykpp'},
      package_data = {'pykpp.models': ['*.eqn', '*.txt']},
      requires = ['numpy', 'scipy', 'pyparsing']
      )
