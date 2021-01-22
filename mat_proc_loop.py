# /usr/bin/env python3

import numpy as np
import os
import fileinput
import glob

# Files:
PI_or_his = 4
obase = '/scratch/e14/rmh561/access-cm2/HCvar/'

if PI_or_his==1:
    base = '/g/data/p66/cm2704/archive/bi889/history/ocn/';
    name = 'PIcontrol'
elif PI_or_his==0:
    base = '/g/data/p66/cm2704/archive/bj594/history/ocn/';
    name = 'historical'
elif PI_or_his==2:
    base = '/g/data/p66/cm2704/archive/by350/history/ocn/';
    name = 'hisNATe1'
elif PI_or_his==3:
    base = '/g/data/p66/cm2704/archive/by438/history/ocn/';
    name = 'hisNATe2'
elif PI_or_his==4:
    base = '/g/data/p66/cm2704/archive/by563/history/ocn/';
    name = 'hisNATe3'

files = glob.glob(base + 'ocean_month.nc-*')
dates = [x[-8:] for x in files]

dates = sorted(dates,key=lambda x: int(x))

# Mask:
msk = ''

# Split into sets of 600 runs for PI control (ss is 0,1):
# ss = 0
# dates = dates[ss*600:ss*600+600]

# Testing (do only 3):
# dates = dates[:3]

print(dates)
print(len(dates))#len(dates))
for i in range(len(dates)):
    fname = base + 'ocean_month.nc-' + dates[i]
    oname = obase + 'CM2_' + name + '_' + msk + '_' + dates[i] + '.mat'
    fscr = 'fscripts/Process_ACCESS_' + msk + '_' + dates[i]
    os.system('cp Process_ACCESS_CM2 ' + fscr)
    with fileinput.FileInput(fscr, inplace=True) as file:
        for line in file:
            line_out = line.replace('XXNAMEXX', 'P' + dates[i]).replace('XXFNAMEXX', fname).replace('XXONAMEXX', oname).replace('XXMSKXX', msk)
            print(line_out, end='')
    os.system('qsub ' + fscr)
            
                    
