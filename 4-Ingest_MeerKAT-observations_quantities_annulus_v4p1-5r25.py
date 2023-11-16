from mega_table_functions import *
import sys
sys.path.append('/Users/cosimaeibensteiner/Documents/GitHub/MegaTable/')
from mega_table.core import StatsTable
from mega_table import RadialMegaTable, TessellMegaTable
from mega_table.utils import nanaverage, nanrms, calc_pixel_per_beam

path_to_existing_base = '/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/scripts/0-data_base/5-Mega-Table_v4p1_5r25/'
path_to_new_base = '/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/scripts/0-data_base/5-Mega-Table_v4p1_5r25_Sigma/'

mhongoose_id = ['IC1954','NGC1566','NGC1672','NGC3511','NGC5068']
phangs_id = ['NGC1512','NGC4535','NGC7496']
adams_id = ['ngc1512','ngc4535','ngc7496']

cosi_mhongoose = [np.cos(57.1/ 180. * np.pi),
np.cos(29.5/ 180. * np.pi),
np.cos(42.6/ 180. * np.pi),
np.cos(75.1/ 180. * np.pi),
np.cos(35.7/ 180. * np.pi)
]

cosi_phangs = [np.cos(42.5/180.*np.pi),
np.cos(44.7/180.*np.pi),
np.cos(35.9/180.*np.pi)]


for m in range(len(mhongoose_id)):
    print(mhongoose_id[m])

    t = RadialMegaTable.read(path_to_existing_base+mhongoose_id[m]+'_base_annulus_1p5kpc.ecsv')

    t = calc_surf_dens_atom(
    t,
    colname='Sigma_atom_MHONGOOSE', unit='Msun pc-2',
    colname_e='e_Sigma_atom_MHONGOOSE', unit_e='Msun pc-2',
    I_HI=t['I_HI'],
    e_I_HI=t['e_I_HI'],
    cosi =cosi_mhongoose[m],
    snr_thresh=3)

    t.write(path_to_new_base+mhongoose_id[m]+'_base_annulus_sn3_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)


for p in range(len(phangs_id)):
    print(phangs_id[p])

    t = RadialMegaTable.read(path_to_existing_base+phangs_id[p]+'_base_annulus_1p5kpc.ecsv')

    t = calc_surf_dens_atom(
        t,
        colname='Sigma_atom_MeerKAT', unit='Msun pc-2',
        colname_e='e_Sigma_atom_MeerKAT', unit_e='Msun pc-2',
        I_HI=t['I_HI'],
        e_I_HI=t['e_I_HI'],
        cosi =cosi_phangs[p],
        snr_thresh=3
        )

    t.write(path_to_new_base+phangs_id[p]+'_base_annulus_sn3_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)
