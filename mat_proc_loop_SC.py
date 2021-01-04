# /usr/bin/env python3

import numpy as np
import os
import fileinput
import glob

# Files:
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

nyr = 100

# Variables:
Tvars3D = ['temp_tendency','temp_submeso',
           'temp_vdiffuse_diff_cbt', 'temp_nonlocal_KPP',
           'temp_vdiffuse_sbc','frazil_3d','sw_heat','temp_rivermix',
           'neutral_diffusion_temp','neutral_gm_temp','temp_vdiffuse_k33',
           'mixdownslope_temp','temp_sigma_diff']
Tvars2D = ['sfc_hflux_pme','temp_eta_smooth','pme_river']
Tvars = Tvars3D + Tvars2D

Svars3D = ['salt_tendency','salt_advection','salt_submeso',
           'salt_vdiffuse_diff_cbt', 'salt_nonlocal_KPP',
           'salt_vdiffuse_sbc','salt_rivermix',
           'neutral_diffusion_salt','neutral_gm_salt','salt_vdiffuse_k33',
           'mixdownslope_salt','salt_sigma_diff']
Svars2D = ['salt_eta_smooth']
Svars = Svars3D + Svars2D

Nvars = ['temp','salt','dht','area_t','geolon_t','geolat_t','st_ocean','st_edges_ocean']
vars = ','.join(Nvars) + ',' + ','.join(Tvars) + ',' + ','.join(Svars)
# print(vars)

set1 = 'ncea -v ' + vars + ' ' + ' '.join([base + 'ocean_month.nc-' + str(x) for x in dates[0:nyr:2]]) + ' ' + obase + 'ocean_month_SC_set1.nc'
set2 = 'ncea -v ' + vars + ' ' + ' '.join([base + 'ocean_month.nc-' + str(x) for x in dates[1:nyr:2]]) + ' ' + obase + 'ocean_month_SC_set2.nc'

# Set 1:
print(len(dates[0:nyr:2]))
print(set1)
print(len(dates[1:nyr:2]))
print(set2)

