from mega_table_functions import *
import sys
sys.path.append('/Users/cosimaeibensteiner/Documents/GitHub/MegaTable/')
from mega_table.core import StatsTable
from mega_table import RadialMegaTable, TessellMegaTable
from mega_table.utils import nanaverage, nanrms, calc_pixel_per_beam

add_area_average_for_image

path_to_moms_meerkat = '/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/Data/for_MHONGOOSE_targets/'
path_to_moms = '/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/Data/PHANGS_sample_cycle0/'
path_to_existing_base = '/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/scripts/0-data_base/4-Mega-Table_v4p0/'
path_to_new_base = '/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/scripts/0-data_base/4-Mega-Table_v4p0/'

mhongoose_id = ['IC1954','NGC1566','NGC1672','NGC3511','NGC5068']
phangs_id = ['NGC1512','NGC4535','NGC7496']
adams_id = ['ngc1512','ngc4535','ngc7496']

cosi_mhongoose = [np.cos(57.1/ 180. * np.pi),
np.cos(29.5/ 180. * np.pi),
np.cos(42.6/ 180. * np.pi),
np.cos(75.1/ 180. * np.pi),
np.cos(35.7/ 180. * np.pi)]

cosi_phangs = [np.cos(42.5/180.*np.pi),
np.cos(44.7/180.*np.pi),
np.cos(35.9/180.*np.pi)]


for m in range(len(mhongoose_id)):

    if m == 4:
        print(mhongoose_id[m],'through if')

        t = TessellMegaTable.read(path_to_existing_base+mhongoose_id[m]+'_base_hexagon_1p5kpc.ecsv')

        t.calc_image_stats(
            #path_to_moms_meerkat+"/old/NGC5068_r05-K_kms_fixed.fits",  # path to data file
            path_to_moms_meerkat+mhongoose_id[m]+"_r05-mom0_v3.fits",  # path to data file
            stat_func=np.nanmean,  # function that calculates stats in each hex
            colname="I_21cm_MHONGOOSE",  # the new column to store the result
            unit='K km/s',  # get physical unit of the result from the file header
            )

        t = add_area_average_for_image(
            t,
            colname='I_21cm_MHONGOOSE', unit='K km/s',
            colname_e='e_I_21cm_MHONGOOSE', unit_e='K km/s',
            #img_file=path_to_moms_meerkat+"/old/NGC5068_r05-K_kms_fixed.fits",
            img_file=path_to_moms_meerkat+mhongoose_id[m]+"_r05-mom0_v3.fits",
            err_file=path_to_moms_meerkat+"NGC5068_r05-emom0.fits")

        t = calc_surf_dens_atom(
            t,
            colname='Sigma_atom_MHONGOOSE', unit='Msun pc-2',
            colname_e='e_Sigma_atom_MHONGOOSE', unit_e='Msun pc-2',
            I_HI=t['I_21cm_MHONGOOSE'],
            e_I_HI=t['e_I_21cm_MHONGOOSE'],
            cosi =cosi_mhongoose[m],
            snr_thresh=0)
°
        #t.write(path_to_new_base+mhongoose_id[m]+'_base+MeerKAT_hexagon_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)


    else:
        print(mhongoose_id[m])

        t = TessellMegaTable.read(path_to_existing_base+mhongoose_id[m]+'_base_hexagon_1p5kpc.ecsv')

        t.calc_image_stats(
            path_to_moms_meerkat+mhongoose_id[m]+"_r05-mom0.fits",  # path to data file
            stat_func=np.nanmean,  # function that calculates stats in each hex
            colname="I_21cm_MHONGOOSE",  # the new column to store the result
            unit='K km/s',  # get physical unit of the result from the file header
            )

        t = add_area_average_for_image(
            t,
            colname='I_21cm_MHONGOOSE', unit='K km/s',
            colname_e='e_I_21cm_MHONGOOSE', unit_e='K km/s',
            img_file=path_to_moms_meerkat+mhongoose_id[m]+"_r05-mom0.fits",
            err_file=path_to_moms_meerkat+mhongoose_id[m]+"_r05-emom0.fits")

        t = calc_surf_dens_atom(
            t,
            colname='Sigma_atom_MHONGOOSE', unit='Msun pc-2',
            colname_e='e_Sigma_atom_MHONGOOSE', unit_e='Msun pc-2',
            I_HI=t['I_21cm_MHONGOOSE'],
            e_I_HI=t['e_I_21cm_MHONGOOSE'],
            cosi =cosi_mhongoose[m],
            snr_thresh=0)

        #t.write(path_to_new_base+mhongoose_id[m]+'_base+MeerKAT_hexagon_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)

    t.write(path_to_new_base+mhongoose_id[m]+'_base+MeerKAT_hexagon_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)
    print('done with ',mhongoose_id[m])

for p in range(len(phangs_id)):
    print(phangs_id[p])

    t = TessellMegaTable.read(path_to_existing_base+phangs_id[p]+'_base_hexagon_1p5kpc.ecsv')

    t.calc_image_stats(
        path_to_moms+adams_id[p]+'_meerkat_hi21cm_pbcorr_zoom_native_k_broad_mom0.fits',
        stat_func=np.nanmean,
        colname='I_21cm_MeerKAT',
        unit='K km/s'
        )

    t = add_area_average_for_image(
        t,
        colname='I_21cm_MeerKAT', unit='K km/s',
        colname_e='e_I_21cm_MeerKAT', unit_e='K km/s',
        img_file=path_to_moms+adams_id[p]+'_meerkat_hi21cm_pbcorr_zoom_native_k_broad_mom0.fits',
        err_file=path_to_moms+adams_id[p]+'_meerkat_hi21cm_pbcorr_zoom_native_k_broad_emom0.fits'
        )

    t = calc_surf_dens_atom(
        t,
        colname='Sigma_atom_MeerKAT', unit='Msun pc-2',
        colname_e='e_Sigma_atom_MeerKAT', unit_e='Msun pc-2',
        I_HI=t['I_21cm_MeerKAT'],
        e_I_HI=t['e_I_21cm_MeerKAT'],
        cosi =cosi_phangs[p],
        snr_thresh=0
        )

    t.write(path_to_new_base+phangs_id[p]+'_base+MeerKAT_hexagon_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)

n5068 = Table.read(path_to_new_base+'NGC5068_base+MeerKAT_hexagon_1p5kpc.ecsv')
print(n5068['Sigma_atom_MHONGOOSE'],np.nanmean(n5068['Sigma_atom_MHONGOOSE']))
