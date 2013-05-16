from distutils.core import setup

help(setup)
exit()
setup(name = 'pykpp',
      version = '0.9',
      author = 'Barron Henderson',
      author_email = 'barronh@gmail.com',
      maintainer = 'Barron Henderson',
      maintainer_email = 'barronh@gmail.com',
      url='https://github.com/barronh/pykpp/',
      packages = ['pykpp'],
      package_dir = {'pykpp': 'src/pykpp'},
      package_data = {'pykpp': ['*.eqn', '*.txt']},
      requires = ['numpy', 'scipy', 'pyparsing']
      )
