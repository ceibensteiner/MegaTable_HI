import sys
import json
import warnings
from pathlib import Path
import pandas as pd

import numpy as np
from astropy import units as u, constants as const
from astropy.io import fits
from astropy.table import Table, QTable


sys.path.append('/Users/ceibenst/Documents/GitHub/MegaTable_HI')
from mega_table.core import StatsTable
from mega_table.table import TessellMegaTable, RadialMegaTable
from mega_table.utils import nanaverage, nanrms, calc_pixel_per_beam

###############################################################################

# location of all relevant config files

config_dir = Path('/Users/ceibenst/Documents/GitHub/MegaTable_HI')
#config_dir = Path('/data/kant/0/sun.1608/PHANGS/mega-tables/code')

# location to save the output data tables
work_dir = Path('/Users/ceibenst/Desktop/home/1-Projects/MEERKAT/scripts/0-data_base/5-Mega-Table_v4p1_5r25')

# logging setting
logging = False

###############################################################################


class PhangsBaseMegaTable(StatsTable):

    """
    MegaTable for PHANGS base data.
    """

    def add_area_average_for_image(
            self, colname=None, unit=None, colname_e=None, unit_e=None,
            img_file=None, err_file=None):

        if not Path(img_file).is_file():
            self[colname] = np.nan * u.Unit(unit)
            if colname_e is not None:
                self[colname_e] = np.nan * u.Unit(unit_e)
            return

        # sample image
        with fits.open(img_file) as hdul:
            data = hdul[0].data.copy()
            hdr = hdul[0].header.copy()
        self.calc_image_stats(
            data, header=hdr, stat_func=nanaverage, weight=None,
            colname=colname, unit='header')
        self[colname] = self[colname].to(unit)

        # sample error image
        if colname_e is None:
            return
        if not Path(err_file).is_file():
            self[colname_e] = np.nan * u.Unit(unit_e)
            return
        with fits.open(err_file) as hdul:
            data_e = hdul[0].data.copy()
        nanflag = np.isnan(data)
        # calculate the direct RMS of errors in each pixel
        self.calc_image_stats(
            data_e, header=hdr, stat_func=nanrms, weight=~nanflag,
            colname=colname_e, unit='header')
        self[colname_e] = self[colname_e].to(unit_e)
        # account for correlated errors among pixels within each beam
        pix_per_beam = calc_pixel_per_beam(hdr)
        hdr.remove('BUNIT', ignore_missing=True)
        self.calc_image_stats(
            (~nanflag).astype('float') / pix_per_beam,
            header=hdr, stat_func=np.sum, colname='_N_beam')
        self['_N_beam'][self['_N_beam'] < 1] = 1
        self[colname_e] /= np.sqrt(self['_N_beam'])
        self.table.remove_column('_N_beam')


    def calc_surf_dens_atom(
            self, colname='Sigma_atom', unit='Msun pc-2',
            colname_e='e_Sigma_atom', unit_e='Msun pc-2',
            I_HI=None, e_I_HI=None,
            cosi=1., e_sys=None, snr_thresh=None):
        alpha_HI = 0.0197 * u.Unit('Msun pc-2 K-1 km-1 s')
        self[colname] = (alpha_HI * cosi * I_HI).to(unit)
        e_stat = alpha_HI * cosi * e_I_HI
        if e_sys is None:
            e_sys = 0.0
        self[colname_e] = np.sqrt(e_stat**2 + e_sys**2).to(unit_e)
        # mask entries below S/N threshold
        if snr_thresh is None:
            snr_thresh = 3
        low_snr_flag = (I_HI / e_I_HI < snr_thresh)
        self[colname][low_snr_flag] = 0
        # self[colname_e][low_snr_flag] = np.nan

    def calc_effective_line_width(
            self, colname='line_width_effective', unit='km s-1',
            I_HI=None, e_I_HI=None, HI_T_PEAK=None, snr_thresh=None):
        self[colname] = (I_HI/(HI_T_PEAK*np.sqrt(2*np.pi))).to(unit)
        # mask entries below S/N threshold
        if snr_thresh is None:
            snr_thresh = 3
        low_snr_flag = (I_HI / e_I_HI < snr_thresh)
        self[colname][low_snr_flag] = 0

    def calc_mom2_line_width(
            self, colname='line_width_mom2', unit='km s-1',
            HI_MOM2=None, I_HI=None,e_I_HI=None, snr_thresh=None):
        #assumed that sqrt(mom2) has been used in the PHANGS pipe to produce mom2 maps [CHECK]
        #self[colname]= np.sqrt(HI_MOM2).to(unit)
        self[colname]= HI_MOM2.to(unit)
        # mask entries below S/N threshold
        if snr_thresh is None:
            snr_thresh = 3
        low_snr_flag = (I_HI / e_I_HI < snr_thresh)
        self[colname][low_snr_flag] = 0



class PhangsBaseTessellMegaTable(
        TessellMegaTable, PhangsBaseMegaTable):

    """
    TessellMegaTable for PHANGS base data.
    """

    def add_deprojected_coords(
            self, colname_r_gal='r_gal', unit_r_gal='kpc',
            colname_phi_gal='phi_gal', unit_phi_gal='deg',
            dist=None, **kwargs):
        from mega_table.utils import deproject
        r_gal, phi_gal = deproject(**kwargs)
        self[colname_r_gal] = (
            r_gal * u.deg / u.rad * dist).to(unit_r_gal)
        self[colname_phi_gal] = (phi_gal * u.deg).to(unit_phi_gal)


class PhangsBaseRadialMegaTable(
        RadialMegaTable, PhangsBaseMegaTable):

    """
    RadialMegaTable for PHANGS base data.
    """

    def add_linear_r_gal(
            self, colname_r_gal='r_gal', unit_r_gal='kpc',
            r_gal_angl_min=None, r_gal_angl_max=None, dist=None):
        self[colname_r_gal] = (
            (r_gal_angl_min + r_gal_angl_max).to('rad').value / 2 *
            dist).to(unit_r_gal)


# -----------------------------------------------------------------------------


def add_raw_measurements_to_table(
        t, data_paths=None, gal_params=None, verbose=True):

    gal_name = gal_params['name']
    gal_ang2lin = gal_params['dist_Mpc'] * u.Mpc / u.rad

    # MeerKAT HI data
    if verbose:
        print("  Add MeerKAT HI data")
        in_file = data_paths['PHANGS_HI_MeerKAT'].format(
            galaxy=gal_name, product='mom0',
            postfix_resolution='')
        err_file = data_paths['PHANGS_HI_MeerKAT'].format(
            galaxy=gal_name, product='emom0',
            postfix_resolution='')
        in_file_t = data_paths['PHANGS_HI_MeerKAT'].format(
            galaxy=gal_name, product='tpeak',
            postfix_resolution='')
        err_file_t = data_paths['PHANGS_HI_MeerKAT'].format(
                galaxy=gal_name, product='e_tpeak',
                postfix_resolution='')
        #--strict masks for mom2 and ew
        in_file_2 = data_paths['PHANGS_HI_MeerKAT_strict'].format(
            galaxy=gal_name, product='mom2',
            postfix_resolution='')
        err_file_2 = data_paths['PHANGS_HI_MeerKAT_strict'].format(
                    galaxy=gal_name, product='emom2',
                    postfix_resolution='')
        in_file_ew = data_paths['PHANGS_HI_MeerKAT_strict'].format(
            galaxy=gal_name, product='ew',
            postfix_resolution='')
        err_file_ew = data_paths['PHANGS_HI_MeerKAT_strict'].format(
                    galaxy=gal_name, product='eew',
                    postfix_resolution='')

        t.add_area_average_for_image(
            # column to save the output
            colname='I_HI_MeerKAT', unit='K km s-1',
            colname_e="e_I_HI_MeerKAT", unit_e='K km s-1',
            # input parameters
            img_file=in_file, err_file=err_file)
        t.add_area_average_for_image(
            # column to save the output
            colname='HI_MeerKAT_T_PEAK', unit='K',
            colname_e="e_HI_MeerKAT_T_PEAK", unit_e='K',
            # input parameters
            img_file=in_file_t, err_file=err_file_t)
        t.add_area_average_for_image(
            # column to save the output
            colname='HI_MeerKAT_MOM2', unit='km s-1',
            colname_e="e_HI_MeerKAT_MOM2", unit_e='km s-1',
            # input parameters
            img_file=in_file_2, err_file=err_file_2)
        t.add_area_average_for_image(
                    # column to save the output
                    colname='HI_MeerKAT_EW', unit='km s-1',
                    colname_e="e_HI_MeerKAT_EW", unit_e='km s-1',
                    # input parameters
                    img_file=in_file_ew, err_file=err_file_ew)

    #if verbose:
        print("  Add MeerKAT HI data from MHONGOOSE")
        in_file = data_paths['MHONGOOSE_HI_MeerKAT'].format(
            galaxy=gal_name, product='mom0',
            postfix_resolution='')
        err_file = data_paths['MHONGOOSE_HI_MeerKAT'].format(
            galaxy=gal_name, product='emom0',
            postfix_resolution='')
    #Erwin didn't provide the whole fits cube nor the mom2 maps (yet)?
        #in_file_2 = data_paths['MHONGOOSE_HI_MeerKAT'].format(
        #    galaxy=gal_name, product='mom2',
        #    postfix_resolution='')
        #err_file_2 = data_paths['MHONGOOSE_HI_MeerKAT'].format(
        #    galaxy=gal_name, product='emom2',
        #    postfix_resolution='')

        t.add_area_average_for_image(
            # column to save the output
            colname='I_MeerKAT_HI', unit='K km s-1',
            colname_e="e_I_MeerKAT_HI", unit_e='K km s-1',
            # input parameters
            img_file=in_file, err_file=err_file)
        #t.add_area_average_for_image(
        #    # column to save the output
        #    colname='HI_MeerKAT_MOM2', unit='km s-1',
        #    colname_e="e_HI_MeerKAT_MOM2", unit_e='km s-1',
        #    # input parameters
        #    img_file=in_file_2, err_file=err_file_2)



def calc_high_level_params_in_table(
        t, gal_params=None, verbose=True):

    gal_cosi = np.cos((gal_params['incl_deg'] * u.deg).to('rad').value)
    gal_ang2lin = gal_params['dist_Mpc'] * u.Mpc / u.rad
    gal_Rstar = (
        gal_params['Rstar_arcsec'] * u.arcsec * gal_ang2lin).to('kpc')
    gal_Mstar = gal_params['Mstar_Msun'] * u.Msun


    # Atomic gas surface density
    if verbose:
        print("  Calculate atom gas surface density")
    t.calc_surf_dens_atom(
        # columns to save the output
        colname="Sigma_atom_MeerKAT", colname_e="e_Sigma_atom_MeerKAT",
        # input parameters
        I_HI=t['I_MeerKAT_HI'], e_I_HI=t['e_I_MeerKAT_HI'],
        cosi=gal_cosi, snr_thresh=0.0)

    # Atomic effective line width
    if verbose:
        print("  Calculate effective line width")
    t.calc_effective_line_width(
            colname='HI_line_width_effective',
            I_HI=t['I_HI_MeerKAT'], e_I_HI=t['e_I_HI_MeerKAT'],
            HI_T_PEAK=t['HI_MeerKAT_T_PEAK'], snr_thresh=0.0)

    #if verbose:
    #    print(" Calculate mom2 line width")
    #t.calc_mom2_line_width(
    #        colname='HI_line_width_mom2',
    #        I_HI=t['I_HI_MeerKAT'], e_I_HI=t['e_I_HI_MeerKAT'],
    #        HI_MOM2=t['HI_MeerKAT_MOM2'], snr_thresh=0.0)


def build_tessell_base_table(
        tile_shape=None, tile_size_kpc=None, fov_radius_R25=None,
        data_paths=None, gal_params=None,
        notes='', version=0.0, output_format=None, writefile=None,
        verbose=True):

    # paraphrase galaxy parameters
    gal_name = gal_params['name']
    gal_ra = np.round(gal_params['ra_deg'], 6) * u.deg
    gal_dec = np.round(gal_params['dec_deg'], 6) * u.deg
    gal_pa = np.round(gal_params['pa_deg'], 1) * u.deg
    gal_incl = np.round(gal_params['incl_deg'], 1) * u.deg
    gal_dist = np.round(gal_params['dist_Mpc'], 2) * u.Mpc
    gal_ang2lin = gal_dist / u.rad
    gal_R25 = np.round((
        gal_params['R25_arcsec'] * u.arcsec * gal_ang2lin).to('kpc'), 2)
    gal_Rstar = np.round((
        gal_params['Rstar_arcsec'] * u.arcsec * gal_ang2lin).to('kpc'), 2)
    gal_Mstar = 10**np.round(
        np.log10(gal_params['Mstar_Msun']), 2) * u.Msun

    # ------------------------------------------------
    # initialize table
    # ------------------------------------------------

    if verbose:
        print("  Initialize tessell base table")
    fov_radius_arcsec = (
        fov_radius_R25 * (gal_R25 / gal_ang2lin).to('arcsec').value)
    tile_size_arcsec = (
        (tile_size_kpc * u.kpc / gal_ang2lin).to('arcsec').value)
    t = PhangsBaseTessellMegaTable(
        # galaxy center coordinates
        gal_ra.to('deg').value, gal_dec.to('deg').value,
        # full field-of-view radius in arcsec
        fov_radius_arcsec,
        # individual tile size in arcsec
        tile_size_arcsec,
        # individual tile shape
        tile_shape=tile_shape)

    # ------------------------------------------------
    # add deprojected galactocentric coordinates
    # ------------------------------------------------

    if verbose:
        print("  Add deprojected galactocentric coordinates")
    t.add_deprojected_coords(
        # columns to save the output
        colname_r_gal='r_gal', colname_phi_gal='phi_gal',
        # input parameters
        ra=t['RA'], dec=t['DEC'],
        center_ra=gal_ra, center_dec=gal_dec,
        incl=gal_incl, pa=gal_pa, dist=gal_dist)

    # ------------------------------------------------
    # add raw measurements & uncertainties
    # ------------------------------------------------

    add_raw_measurements_to_table(
        t, data_paths=data_paths, gal_params=gal_params,
        verbose=verbose)

    # ------------------------------------------------
    # calculate high-level parameters
    # ------------------------------------------------

    calc_high_level_params_in_table(
        t, gal_params=gal_params, verbose=verbose)

    # ------------------------------------------------
    # clean and format output table
    # ------------------------------------------------

    t.table.sort(['r_gal', 'phi_gal'])
    t['ID'] = np.arange(len(t))
    colnames = (
        ['ID', 'RA', 'DEC', 'r_gal', 'phi_gal'] +
        [str(entry) for entry in output_format['colname']])
    units = (
        ['', 'deg', 'deg', 'kpc', 'deg'] +
        [str(entry) for entry in output_format['unit']])
    formats = (
        ['.0f', '.5f', '.5f', '.2f', '.2f'] +
        [str(entry) for entry in output_format['format']])
    descriptions = (
        ["Aperture ID",
         "Right Ascension of the aperture center",
         "Declination of the aperture center",
         "Deprojected galactocentric radius",
         "Deprojected azimuthal angle (0 = receding major axis)"] +
        [str(entry) for entry in output_format['description']])
    t.format(
        colnames=colnames, units=units, formats=formats,
        descriptions=descriptions, ignore_missing=True)

    # ------------------------------------------------
    # add metadata
    # ------------------------------------------------

    t.meta['GALAXY'] = str(gal_name)
    t.meta['RA_DEG'] = float(np.round(gal_ra.to('deg').value, 6))
    t.meta['DEC_DEG'] = float(np.round(gal_dec.to('deg').value, 6))
    t.meta['INCL_DEG'] = float(np.round(gal_incl.to('deg').value, 1))
    t.meta['PA_DEG'] = float(np.round(gal_pa.to('deg').value, 1))
    t.meta['DIST_MPC'] = float(np.round(gal_dist.to('Mpc').value, 2))
    t.meta['R25_KPC'] = float(np.round(gal_R25.to('kpc').value, 2))
    t.meta['RSCL_KPC'] = float(np.round(gal_Rstar.to('kpc').value, 2))
    t.meta['LOGMSTAR'] = float(
        np.round(np.log10(gal_Mstar.to('Msun').value), 2))
    t.meta['TBLNOTE'] = str(notes)
    t.meta['VERSION'] = float(version)

    # ------------------------------------------------
    # write table
    # ------------------------------------------------

    if writefile is not None:
        if verbose:
            print("  Write table")
        if Path(writefile).suffix == '.ecsv':
            t.write(
                writefile, add_timestamp=True,
                delimiter=',', overwrite=True)
        else:
            t.write(writefile, add_timestamp=True, overwrite=True)
        return writefile
    else:
        return t


def build_radial_base_table(
        annulus_width_kpc=None, fov_radius_R25=None,
        data_paths=None, gal_params=None,
        notes='', version=0.0, output_format=None, writefile=None,
        verbose=True):

    # paraphrase galaxy parameters
    gal_name = gal_params['name']
    gal_ra = np.round(gal_params['ra_deg'], 6) * u.deg
    gal_dec = np.round(gal_params['dec_deg'], 6) * u.deg
    gal_pa = np.round(gal_params['pa_deg'], 1) * u.deg
    gal_incl = np.round(gal_params['incl_deg'], 1) * u.deg
    gal_dist = np.round(gal_params['dist_Mpc'], 2) * u.Mpc
    gal_ang2lin = gal_dist / u.rad
    gal_R25 = np.round((
        gal_params['R25_arcsec'] * u.arcsec * gal_ang2lin).to('kpc'), 2)
    gal_Rstar = np.round((
        gal_params['Rstar_arcsec'] * u.arcsec * gal_ang2lin).to('kpc'), 2)
    gal_Mstar = 10**np.round(
        np.log10(gal_params['Mstar_Msun']), 2) * u.Msun

    # ------------------------------------------------
    # initialize table
    # ------------------------------------------------

    if verbose:
        print("  Initialize table")
    rgal_bin_arcsec = (
        (annulus_width_kpc * u.kpc / gal_ang2lin).to('arcsec').value)
    rgal_max_arcsec = (
        fov_radius_R25 * (gal_R25 / gal_ang2lin).to('arcsec').value)
    t = PhangsBaseRadialMegaTable(
        # galaxy center coordinates
        gal_ra.to('deg').value, gal_dec.to('deg').value,
        # radial bin width in arcsec
        rgal_bin_arcsec,
        # full field-of-view radius in arcsec
        rgal_max_arcsec=rgal_max_arcsec,
        # galaxy orientation parameters
        gal_incl_deg=gal_incl.to('deg').value,
        gal_posang_deg=gal_pa.to('deg').value)

    # ------------------------------------------------
    # add galactocentric radius in linear size
    # ------------------------------------------------

    if verbose:
        print("  Add galactocentric radius in linear size")
    t.add_linear_r_gal(
        # columns to save the output
        colname_r_gal='r_gal',
        # input parameters
        r_gal_angl_min=t['r_gal_angl_min'],
        r_gal_angl_max=t['r_gal_angl_max'],
        dist=gal_dist)

    # ------------------------------------------------
    # add raw measurements & uncertainties
    # ------------------------------------------------

    add_raw_measurements_to_table(
        t, data_paths=data_paths, gal_params=gal_params,
        verbose=verbose)

    # ------------------------------------------------
    # calculate high-level parameters
    # ------------------------------------------------

    calc_high_level_params_in_table(
        t, gal_params=gal_params, verbose=verbose)

    # ------------------------------------------------
    # clean and format output table
    # ------------------------------------------------

    t['ID'] = np.arange(len(t))
    colnames = (
        ['ID', 'r_gal'] +
        [str(entry) for entry in output_format['colname']])
    units = (
        ['', 'kpc'] +
        [str(entry) for entry in output_format['unit']])
    formats = (
        ['.0f', '.2f'] +
        [str(entry) for entry in output_format['format']])
    descriptions = (
        ["Annulus ID",
         "Deprojected galactocentric radius"] +
        [str(entry) for entry in output_format['description']])
    t.format(
        colnames=colnames, units=units, formats=formats,
        descriptions=descriptions, ignore_missing=True)

    # ------------------------------------------------
    # add metadata
    # ------------------------------------------------

    t.meta['GALAXY'] = str(gal_name)
    t.meta['RA_DEG'] = float(np.round(gal_ra.to('deg').value, 6))
    t.meta['DEC_DEG'] = float(np.round(gal_dec.to('deg').value, 6))
    t.meta['INCL_DEG'] = float(np.round(gal_incl.to('deg').value, 1))
    t.meta['PA_DEG'] = float(np.round(gal_pa.to('deg').value, 1))
    t.meta['DIST_MPC'] = float(np.round(gal_dist.to('Mpc').value, 2))
    t.meta['R25_KPC'] = float(np.round(gal_R25.to('kpc').value, 2))
    t.meta['RSCL_KPC'] = float(np.round(gal_Rstar.to('kpc').value, 2))
    t.meta['LOGMSTAR'] = float(
        np.round(np.log10(gal_Mstar.to('Msun').value), 2))
    t.meta['TBLNOTE'] = str(notes)
    t.meta['VERSION'] = float(version)

    # ------------------------------------------------
    # write table
    # ------------------------------------------------

    if writefile is not None:
        if verbose:
            print("  Write table")
        if Path(writefile).suffix == '.ecsv':
            t.write(
                writefile, add_timestamp=True,
                delimiter=',', overwrite=True)
        else:
            t.write(writefile, add_timestamp=True, overwrite=True)
        return writefile
    else:
        return t


# -----------------------------------------------------------------------------


if __name__ == '__main__':

    # warning and logging settings
    warnings.simplefilter('ignore', RuntimeWarning)
    if logging:
        # output log to a file
        orig_stdout = sys.stdout
        log = open(work_dir / (str(Path(__file__).stem) + '.log'), 'w')
        sys.stdout = log
    else:
        orig_stdout = log = None

    # load config files
    with open(config_dir / "config_data_path.json") as f:
        data_paths = json.load(f)
    with open(config_dir / "config_tables.json") as f:
        table_configs = json.load(f)
    t_format = Table.read(config_dir / "format_base.csv")

    # read sample table
    #columns=['name','orient_ra','orient_dec','orient_incl','orient_posang','dist','props_mstar','size_r25','size_r25_unc','size_scalelength','size_scalelength_unc']


    #t_sample = pd.read_csv(data_paths['PHANGS_sample_table'],names=columns, sep=';', keep_default_na=False,comment='#')

    t_sample = Table.read(data_paths['PHANGS_sample_table'])
    #print(t_sample["name"])

    # sub-select sample
    t_sample = t_sample[t_sample['data_has_megatable']]

    # loop through all galaxies
    for row in t_sample:

        # extract galaxy parameters
        gal_params = {
            'name': row['name'].upper(),
            'dist_Mpc': row['dist'],
            'ra_deg': row['orient_ra'],
            'dec_deg': row['orient_dec'],
            'incl_deg': row['orient_incl'],
            'pa_deg': row['orient_posang'],
            'Rstar_arcsec': row['size_scalelength'],
            'R25_arcsec': row['size_r25'],
            'Mstar_Msun': row['props_mstar'],
        }

        print("\n############################################################")
        print(f"# {gal_params['name']}")
        print("############################################################\n")

        # TessellMegaTable
        tile_shape = table_configs['tessell_tile_shape']
        tile_size = (
            table_configs['tessell_tile_size'] *
            u.Unit(table_configs['tessell_tile_size_unit']))
        tile_size_str = (
            str(table_configs['tessell_tile_size']).replace('.', 'p') +
            table_configs['tessell_tile_size_unit'])
        fov_radius_R25 = table_configs['tessell_FoV_radius']

        tessell_base_table_file = (
            work_dir / table_configs['tessell_table_name'].format(
                galaxy=gal_params['name'], content='base',
                tile_shape=tile_shape, tile_size_str=tile_size_str))
        if not tessell_base_table_file.is_file():
            print("Building tessellation statistics table...")
            build_tessell_base_table(
                tile_shape=tile_shape,
                tile_size_kpc=tile_size.to('kpc').value,
                fov_radius_R25=fov_radius_R25,
                data_paths=data_paths,
                gal_params=gal_params,
                version=table_configs['table_version'],
                notes=table_configs['table_notes'],
                output_format=t_format,
                writefile=tessell_base_table_file)
            print("Done\n")

        # RadialMegaTable
        annulus_width = (
            table_configs['radial_annulus_width'] *
            u.Unit(table_configs['radial_annulus_width_unit']))
        annulus_width_str = (
            str(table_configs['radial_annulus_width']).replace('.', 'p') +
            table_configs['radial_annulus_width_unit'])
        fov_radius_R25 = table_configs['radial_FoV_radius']

        radial_base_table_file = (
            work_dir / table_configs['radial_table_name'].format(
                galaxy=gal_params['name'], content='base',
                annulus_width_str=annulus_width_str))
        if not radial_base_table_file.is_file():
            print("Building radial statistics table...")
            build_radial_base_table(
                annulus_width_kpc=annulus_width.to('kpc').value,
                fov_radius_R25=fov_radius_R25,
                data_paths=data_paths,
                gal_params=gal_params,
                version=table_configs['table_version'],
                notes=table_configs['table_notes'],
                output_format=t_format,
                writefile=radial_base_table_file)
            print("Done\n")

    # logging settings
    if logging:
        sys.stdout = orig_stdout
        log.close()
