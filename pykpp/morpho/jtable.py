import os
from numpy import *
from datetime import datetime

months = dict(MY = 5,
              JN = 6,
              JL = 7,
              AU = 8,
              ST = 9,
              OC = 10,)

def read_jtable(path, seconds_since_start = True):
    fname = os.path.basename(path)
    month = months[fname[:2]]
    day = int(fname[2:4])
    year = int(fname[4:6])
    if year < 80:
        year += 2000
    
    lines = file(path).readlines()
    location = 'header'
    keys = []
    units = []
    temp = ''
    values = []
    times = []
    for line in lines:
        if location == 'header':
            if line.strip() == '<':
                location = 'meta'
        elif location == 'meta':
            if line.strip() == '>':
                location = 'data'
            else:
                key, unit = map(str.strip, line.split(';')[0].split(','))
                keys.append(key)
                units.append(unit)
        elif location == 'data':
            if line.strip() in ('{', '}'):
                continue
            else:
                temp += line
            if '}' in line:
                temp = temp.replace('}', '').replace('{', '').strip()
                vals = temp.split()
                timestr = vals.pop(0)
                time = datetime(year, month, day, int(timestr[:2]), int(timestr[2:4]))
                vals = map(float, vals)
                vals = [time] + vals
                times.append(time)
                values.append(vals)
                temp = ''
    names = keys
    formats = ['f'] * len(vals)
    if not seconds_since_start:
        formats[0] = 'object'
    output = zeros(len(values), dtype = dtype(dict(names = names, formats = formats)))
    start_time = values[0][keys.index('TIME')]
    for ri, row in enumerate(values):
        for i, (key, unit, val) in enumerate(zip(keys, units, row)):        
            factor = {'PER MIN': 60., '1000*PER MIN': 60. / 1000, 'HHMM': 1}[unit]
            if factor != 1:
                val *= factor
            if seconds_since_start:
                if key == 'TIME':
                    val = (val - start_time).total_seconds()
            output[key][ri] = val
    return output

def get_jtable_funcs(path):
    """
    Creates two functions
    Update_JTABLE(mech, world) - prepares JTABLE values for a 
                                 time (t) in world
    JTABLE(key) - Returns jvalue in 1/s for time at last Update_JTABLE
    """
    
    def Update_JTABLE(mech, world):
        global joutput
        global jtimes
        global jvalues
        try:
            len(joutput)
        except:
            joutput = read_jtable(path, seconds_since_start = True)
            jtimes = joutput['TIME']
        t = world['t']
        jvalues = dict([(k, interp(t, jtimes, joutput[k])) for k in joutput.dtype.names])
    
    def JTABLE(key):
        global jvalues
        out = jvalues[key]
        return out
    return Update_JTABLE, JTABLE