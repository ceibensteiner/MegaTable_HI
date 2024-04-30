import sys
sys.path.append('/Users/cosimaeibensteiner/Desktop/home/PhD/0-Coding/general_scripts/')
import numpy as np
from astropy import constants as const
from pathlib import Path
from astropy.io import fits
from astropy import units as u, constants as const
from astropy.table import Table, QTable

sys.path.append('/Users/cosimaeibensteiner/Documents/GitHub/MegaTable/')
from mega_table.utils import nanaverage, nanrms#, calc_pixel_per_beam


def calc_pixel_per_beam(header, suppress_no_beam_error=True):
    from astropy.wcs import WCS
    from radio_beam import Beam
    from radio_beam.beam import NoBeamException
    try:

        beam = Beam.from_fits_header(header)
        wcs = WCS(header)


        pixsize = (header['CDELT2'])**2*u.sr


        return (beam.sr / pixsize)

        #return (beam.sr / wcs.proj_plane_pixel_area()).to('').value
    except NoBeamException as e:
        if suppress_no_beam_error:
            return None
        else:
            raise NoBeamException(e)


def add_area_average_for_image(
        table, colname=None, unit=None, colname_e=None,
        unit_e=None, img_file=None, err_file=None):

    if not Path(img_file).is_file():
        print("Input image file not found - Default to NaN")
        table[colname] = np.nan * u.Unit(unit)
        if colname_e is not None:
            table[colname_e] = np.nan * u.Unit(unit_e)
        return

    # sample image
    with fits.open(img_file) as hdul:
        data = hdul[0].data.copy()
        hdr = hdul[0].header.copy()
    table.calc_image_stats(
        data, header=hdr, stat_func=nanaverage, weight=None,
        colname=colname, unit='header')
    table[colname] = table[colname].to(unit)

    # sample error image
    if colname_e is None:
        return
    if not Path(err_file).is_file():
        print("Input error file not found - Default to NaN")
        table[colname_e] = np.nan * u.Unit(unit_e)
        return
    with fits.open(err_file) as hdul:
        data_e = hdul[0].data.copy()
    nanflag = np.isnan(data)

    # calculate the direct RMS of errors in each pixel
    table.calc_image_stats(
        data_e, header=hdr, stat_func=nanrms, weight=~nanflag,
        colname=colname_e, unit='header')
    table[colname_e] = table[colname_e].to(unit_e)

    # account for correlated errors among pixels within each beam
    pix_per_beam = calc_pixel_per_beam(hdr)
    hdr.remove('BUNIT', ignore_missing=True)
    table.calc_image_stats(
        (~nanflag).astype('float') / pix_per_beam,
        header=hdr, stat_func=np.sum, colname='_N_beam')
    table['_N_beam'][table['_N_beam'] < 1] = 1
    table[colname_e] /= np.sqrt(table['_N_beam'])
    #table.remove_column('_N_beam')

    return table


def calc_surf_dens_atom(
            table, colname='Sigma_atom', unit='Msun pc-2',
            colname_e='e_Sigma_atom', unit_e='Msun pc-2',
            I_HI=None, e_I_HI=None,
            cosi=1., e_sys=None, snr_thresh=None):

        alpha_HI = 0.0197 * u.Unit('Msun pc-2 K-1 km-1 s')
        table[colname] = (alpha_HI * cosi * I_HI).to(unit)
        e_stat = alpha_HI * cosi * e_I_HI

        if e_sys is None:
            e_sys = 0.0
        table[colname_e] = np.sqrt(e_stat**2 + e_sys**2).to(unit_e)
        # mask entries below S/N threshold
        if snr_thresh is None:
            snr_thresh = 3
        low_snr_flag = (I_HI / e_I_HI < snr_thresh)
        table[colname][low_snr_flag] = 0
        # self[colname_e][low_snr_flag] = np.nan

        return table


def stellar_mass_volume_density(Sigma_star, e_Sigma_star, radial_scale_lenght):
    ''' Parameters
    Sigma_star: from MegaTable
    radial_scale_lenght: from Leory+21 Sample table
    '''

    #Sigma_star = Sigma_star.to(u.kg/u.m**2)

    radial_scale_lenght = radial_scale_lenght*u.kpc
    radial_scale_lenghtpc = radial_scale_lenght.to(u.pc)

    unit = 'Msun pc-3'
    #we need h_star
    flattening = 7.3 #Kregel+02; Sun+20a
    #rho = Sigma_star/(0.54*radial_scale_lenght).to(unit) # Eq. 13 Sun+20
    h_star = radial_scale_lenghtpc / flattening # flat method
    rho = (Sigma_star / 4 / h_star).to(unit)

    #calculate error
    unit_e = 'Msun pc-3'
    e_stat = e_Sigma_star / 4 / h_star
    e_sys = 0.0
    rho_e = np.sqrt(e_stat**2 + e_sys**2).to(unit_e)

    return rho, rho_e

def calc_dyn_eq_p(Sigma_star=None, radial_scale_lenght=None,
                Sigma_mol=None, Sigma_atom=None, vdisp_atom_z=None, vdisp_mol_z=None,
                e_Sigma_mol=None, e_Sigma_atom=None, e_Sigma_star=None):

    rho,rho_e  = stellar_mass_volume_density(Sigma_star, e_Sigma_star, radial_scale_lenght)

    if vdisp_atom_z is None:
        vdisp_atom_z = 10 * u.Unit('km s-1')  # Leroy+08, Sun+20a
    if vdisp_mol_z is None:
        vdisp_mol_z = vdisp_atom_z  # Leroy+08

    Sigma_gas = Sigma_atom + Sigma_mol
    vdisp_gas_z = (  # Sun+20a Eq.14
        Sigma_mol * vdisp_mol_z +
        Sigma_atom * vdisp_atom_z) / Sigma_gas

    vdisp_gas_z = vdisp_gas_z.to('km s-1')
    print(vdisp_gas_z)
    Sigma_gas = Sigma_gas.to(u.kg/u.m**2)
    first_term = ((np.pi*const.G)/2)*Sigma_gas**2
    second_term = Sigma_gas*np.sqrt(2*const.G*rho)*vdisp_gas_z
    pressure = first_term+second_term

    p = (pressure/const.k_B).to('K cm-3')

    #calculate error
    e_Sigma_gas = np.sqrt(e_Sigma_mol**2 + e_Sigma_atom**2)
    e_Sigma_gas = e_Sigma_gas.to(u.kg/u.m**2)
    e_stat = np.sqrt(
            (np.pi * const.G * Sigma_gas +
             vdisp_gas_z * np.sqrt(2 * const.G * rho))**2 *
             e_Sigma_gas**2 +
            (Sigma_gas * vdisp_gas_z *
             np.sqrt(const.G / 2 / rho))**2 *
            rho_e **2) / const.k_B
    e_sys = 0.0
    p_e = np.sqrt(e_stat**2 + e_sys**2).to('K cm-3')

    return p, p_e








def pde(Sigma_star, radial_scale_lenght, Sigma_mol, Sigma_atom, sigma_hi_z, e_Sigma_mol, e_Sigma_atom, e_Sigma_star):
    '''Parameters
    Sigma_star: from MegaTable
    radial_sclae_lenght: from Sato+2015 (S4G), Table 1, column 6
    Sigma_gas: from MegaTable Sigma_mol+Sigma_HI
    sigma_hi_z: for now fixed to 11km/s
    e_Sigma_mol: from MegaTable
    e_Sigma_atom: from MegaTable
    '''

    rho,rho_e  = stellar_mass_volume_density(Sigma_star, e_Sigma_star, radial_scale_lenght)

    vdisp_atom_z = sigma_hi_z   # Leroy+08, Sun+20a == 10 km/s
    vdisp_mol_z = sigma_hi_z  # Leroy+08
    Sigma_gas = Sigma_atom + Sigma_mol

    vdisp_gas_z = (  # Sun+20a Eq.14
            Sigma_mol * vdisp_mol_z +
            Sigma_atom * vdisp_atom_z) / Sigma_gas

    #sigma_gas_z = sigma_gas_z*(u.km/u.s)
    vdisp_gas_z = vdisp_gas_z.to('km s-1')
    Sigma_gas = Sigma_gas.to(u.kg/u.m**2)
    first_term = ((np.pi*const.G)/2)*Sigma_gas**2
    second_term = Sigma_gas*np.sqrt(2*const.G*rho)*vdisp_gas_z

    pressure = first_term+second_term
    print(pressure)

    #let's get it in units of K/cm-3 (dividing by k_B ==> in Jiayi's code for paper2022; different than in paper2020)
    unit = 'K cm-3'# unit is K cm-3 kB !!
    pressure = (pressure/const.k_B).to(unit)
    print(pressure)


    #calculate error
    unit_e = 'K cm-3'
    e_Sigma_gas = np.sqrt(e_Sigma_mol**2 + e_Sigma_atom**2)
    e_Sigma_gas = e_Sigma_gas.to(u.kg/u.m**2)
    e_stat = np.sqrt(
            (np.pi * const.G * Sigma_gas +
             vdisp_gas_z * np.sqrt(2 * const.G * rho))**2 *
             e_Sigma_gas**2 +
            (Sigma_gas * vdisp_gas_z *
             np.sqrt(const.G / 2 / rho))**2 *
            rho_e **2) / const.k_B
    e_sys = 0.0
    pressure_e = np.sqrt(e_stat**2 + e_sys**2).to(unit_e)

    return pressure, pressure_e



#--- varing vdisp
def pde_varying(Sigma_star, radial_scale_lenght,
                Sigma_mol, Sigma_atom, sigma_hi_z, sigma_mol_z,
                e_Sigma_mol, e_Sigma_atom, e_Sigma_star):

    unit='K cm-3',
    unit_e='K cm-3'
    vdisp_atom_z = sigma_hi_z.to('km s-1')
    vdisp_mol_z = sigma_mol_z.to('km s-1')
    #print(vdisp_atom_z)
    #print(vdisp_mol_z)
    Sigma_gas = Sigma_atom + Sigma_mol
    rho,rho_e  = stellar_mass_volume_density(Sigma_star, e_Sigma_star, radial_scale_lenght)
    vdisp_gas_z = (  # Sun+20a Eq.14
            Sigma_mol * vdisp_mol_z +
            Sigma_atom * vdisp_atom_z) / Sigma_gas
    vdisp_gas_z = vdisp_gas_z.to('km s-1')

    Sigma_gas = Sigma_gas.to(u.kg/u.m**2)

    first_term = ((np.pi*const.G)/2)*Sigma_gas**2
    second_term = Sigma_gas*np.sqrt(2*const.G*rho)*vdisp_gas_z

    pressure = first_term+second_term

    #let's get it in units of K/cm-3 (dividing by k_B ==> in Jiayi's code for paper2022; different than in paper2020)
    unit = 'K cm-3'# unit is K cm-3 kB !!
    pressure_ = (pressure/const.k_B).to(unit)

    e_Sigma_gas = np.sqrt(e_Sigma_mol**2 + e_Sigma_atom**2)
    e_Sigma_gas = e_Sigma_gas.to(u.kg/u.m**2)
    e_stat = np.sqrt(
            (np.pi * const.G * Sigma_gas +
             vdisp_gas_z * np.sqrt(2 * const.G * rho))**2 *
             e_Sigma_gas**2 +
            (Sigma_gas * vdisp_gas_z *
             np.sqrt(const.G / 2 / rho))**2 *
            rho_e **2) / const.k_B
    e_sys = 0.0
    pressure_e_ = np.sqrt(e_stat**2 + e_sys**2).to(unit_e)
    #unit = 'K cm-3'
    #pressure_e = (pressure_e/const.k_B).to(unit)

    return pressure_, pressure_e_
