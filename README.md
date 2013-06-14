pykpp
=====

pykpp is a KPP-like chemical mechanism parser that produces a box model solvable by SciPy's odeint solver

Running Instructions
--------------------

1. Open DOS Prompt on Windows or Terminal on Linux/Mac OSX
2. Type "python -m pykpp model.kpp" where model.kpp is a valid KPP input file.
3. Example below assumes Installation and that cbm4.eqn has been downloaded to the downloads folder

```
C:\Users\MyUserName\Downloads>python -m pykpp cbm4.eqn
C:\Python27\lib\site-packages\pykpp\mech.py:138: UserWarning: Ignoring PROD
  warn('Ignoring %s' % spc)
Species: 33
Reactions: 81
C:\Python27\lib\site-packages\pykpp\mech.py:286: UserWarning: Using SunRise and SunSet of 4.5 and 19.5 (approximately JulianDay 145 and Latitude 45 degrees N)
warn('Using SunRise and SunSet of 4.5 and 19.5 (approximately JulianDay 145 and Latitude 45 degrees N)')
0.0%; {T=4.320E+04 O3=1.000E+02 NO=5.000E+01 }
0.8%; {T=4.680E+04 O3=9.770E+01 NO=9.916E+00 }
1.7%; {T=5.045E+04 O3=1.341E+02 NO=4.216E+00 }
...
98.7%; {T=4.695E+05 O3=1.030E+02 NO=1.393E-01 }
99.5%; {T=4.732E+05 O3=1.032E+02 NO=1.411E-01 }
100.0%; {T=4.752E+05 O3=1.034E+02 NO=1.411E-01 }
Solved in 11.307000 seconds
```

* Example input files can be downloaded from https://github.com/barronh/pykpp/tree/master/src/pykpp/models and are also in the models folder of pykpp source

Install Instructions
--------------------

1. Windows
  1. Install python-xy using Windows executable
  2. Open DOS command prompt (start run -> type "cmd" -> hit enter
  3. type "pip install https://github.com/barronh/pykpp/archive/master.zip
2. Mac OSX
  1. Install the Enthough Python Distribution
  2. Install following "Get Pip" instructions at http://www.pip-installer.org/en/latest/installing.html
  3. Open a Terminal
  4. type "sudo pip install https://github.com/barronh/pykpp/archive/master.zip
5. Ubuntu Linux
  1. Open Terminal
  2. sudo apt-get install python-scipy
  3. sudo apt-get install python-pip
  3. sudo pip install http://cheeseshop.python.org/packages/source/p/pyparsing/pyparsing-1.5.5.tar.gz
  4. type "sudo pip install https://github.com/barronh/pykpp/archive/master.zip


