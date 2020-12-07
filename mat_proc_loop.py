# /usr/bin/env python3

import numpy as np
import os
import fileinput
import glob

PI_or_his = 1
obase = '/scratch/e14/rmh561/access-cm2/HCvar/'

if PI_or_his==1:
    base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';
    name = 'PIcontrol'
else:
    base = '/g/data/p66/cm2704/archive/bj594/history/ocn/';
    name = 'historical'

files = glob.glob(base + 'ocean_month.nc-*')
dates = [x[-8:] for x in files]

dates = sorted(dates,key=lambda x: int(x))

dates = dates[:10]
print(dates)
for i in range(len(dates)):
    fname = base + 'ocean_month.nc-' + dates[i]
    oname = obase + 'CM2_' + name + '_' + dates[i] + '.mat'
    fscr = 'fscripts/Process_ACCESS_' + dates[i]
    os.system('cp Process_ACCESS_CM2 ' + fscr)
    with fileinput.FileInput(fscr, inplace=True) as file:
        for line in file:
            line_out = line.replace('XXNAMEXX', 'P' + dates[i]).replace('XXFNAMEXX', fname).replace('XXONAMEXX', oname)
            print(line_out, end='')
    os.system('qsub ' + fscr)
            
                    