import sys
sys.path.append('/Users/ceibenst/Desktop/home/1-Projects/general_scripts/')
from tools import *
from astropy.table import Table

'''Function to get new database for megatable based database'''


def replace_zeros(x):
    "Function to replace zeros with nan"
    x_mask = x.copy()
    x_mask[x_mask==0] = 'nan'

    return x_mask

def sn_d_faint(line1, ucs1):
    SN_stack = line1/ucs1
    detect_stack     = np.where(SN_stack>=3)
    return detect_stack

def sn_n_faint(line1, ucs1):
    SN_stack = line1/ucs1
    not_detect_stack = np.where(SN_stack<3)
    return not_detect_stack

def get_new_base(base_dir,base,galaxyID,base_dir_out):
    '''The main fuction'''

    gal = Table.read(base_dir+base)

    #get quantities
    hi = gal["I_21cm_MeerKAT"]
    hi_unc = gal["e_I_21cm_MeerKAT"]
    co = gal["I_CO21"]
    co_unc = gal["e_I_CO21"]
    rgal_kpc = gal["r_gal"]

    hi = replace_zeros(hi)
    co = replace_zeros(co)

    ints_f   = np.array([hi,co])
    ucs_f    = np.array([hi_unc,co_unc])
    lines_f  = np.array(['I_HI','I_CO21'])
    length = len(hi)

    #new database
    database = gal

    for i in range(2):

        #-------------------------------------------------------------------
        # some basic infos first
        #
        detection = sn_d_faint(ints_f[i],ucs_f[i])
        not_significant = sn_n_faint(ints_f[i],ucs_f[i])

        #-------------------------------------------------------------------
        # now getting the upperlimtis
        npoints    = np.arange(0,length,dtype=np.float64)
        upperlim   = 3*ucs_f[i]

        # -- inds for detected and non detection
        detect = npoints[detection]
        notdetect = npoints[not_significant]

        # -- Use the upperlimits:  set to False for no upperlimit (bool_limits)
        notdelta = np.isin(npoints,notdetect)
        database["notdelta_"+str(lines_f[i])] = notdelta

        # -- Not using upperlimits: set to False for upperlimit (bool_not_limits)
        delta    = np.logical_not(notdelta)
        database["delta_"+str(lines_f[i])] = delta

        # -- use calculated upperlims for non detected and fill with nans for detected
        upperlim_nan = upperlim.copy()
        upperlim_nan[delta] = np.nan
        database["upperlim_nan"+str(lines_f[i])] = upperlim_nan

        # -- create mask for: upper limits
        mask_forups  = np.ma.masked_where(database["notdelta_"+str(lines_f[i])],ints_f[i][:])
        database["mask_forups_"+str(lines_f[i])] = mask_forups
        ind_limits         = np.where(database["notdelta_"+str(lines_f[i])]) #ind similar to notdet
        ind_not_limits     = np.where(database["delta_"+str(lines_f[i])])    #ind_not similar to detect

        # -- replace with the numbers of the upperlimits
        int_new = ints_f[i]
        numb    = ind_limits[0]
        upperlim_nan = database["upperlim_nan"+str(lines_f[i])]
        int_new[numb] = upperlim_nan[numb]
        database["ycens"+str(lines_f[i])] = int_new

    fname_dict = base_dir_out+galaxyID+'_base+MeerKAT+P-DE_hexagon_1p5kpc_uplims_ints'+'.ecsv'
    database.write(fname_dict, format='ascii.ecsv',overwrite=True)

    return database


def get_new_base_mhongoose(base_dir,base,galaxyID,base_dir_out):
    '''The main fuction // for MHONGOOSE targets'''

    gal = Table.read(base_dir+base)

    #get quantities
    hi = gal["I_21cm_MHONGOOSE"]
    hi_unc = gal["e_I_21cm_MHONGOOSE"]
    co = gal["I_CO21"]
    co_unc = gal["e_I_CO21"]
    rgal_kpc = gal["r_gal"]

    hi = replace_zeros(hi)
    co = replace_zeros(co)

    ints_f   = np.array([hi,co])
    ucs_f    = np.array([hi_unc,co_unc])
    lines_f  = np.array(['I_HI','I_CO21'])
    length = len(hi)

    #new database
    database = gal

    for i in range(2):

        #-------------------------------------------------------------------
        # some basic infos first
        #
        detection = sn_d_faint(ints_f[i],ucs_f[i])
        not_significant = sn_n_faint(ints_f[i],ucs_f[i])

        #-------------------------------------------------------------------
        # now getting the upperlimtis
        npoints    = np.arange(0,length,dtype=np.float64)
        upperlim   = 3*ucs_f[i]

        # -- ind for detected
        detect = npoints[detection]
        notdetect = npoints[not_significant]

        # -- Use the upperlimits:  set to False for no upperlimit (bool_limits)
        notdelta = np.isin(npoints,notdetect)
        database["notdelta_"+str(lines_f[i])] = notdelta

        # -- Not using upperlimits: set to False for upperlimit (bool_not_limits)
        delta    = np.logical_not(notdelta)
        database["delta_"+str(lines_f[i])] = delta

        # -- use calculated upperlims for non detected and fill with nans for detected
        upperlim_nan = upperlim.copy()
        upperlim_nan[delta] = np.nan
        database["upperlim_nan"+str(lines_f[i])] = upperlim_nan

        # -- create mask for: upper limits
        mask_forups  = np.ma.masked_where(database["notdelta_"+str(lines_f[i])],ints_f[i][:])
        database["mask_forups_"+str(lines_f[i])] = mask_forups
        ind_limits         = np.where(database["notdelta_"+str(lines_f[i])]) #ind similar to notdet
        ind_not_limits     = np.where(database["delta_"+str(lines_f[i])])    #ind_not similar to detect

        # -- replace with the numbers of the upperlimits
        int_new = ints_f[i]
        numb    = ind_limits[0]
        upperlim_nan = database["upperlim_nan"+str(lines_f[i])]
        int_new[numb] = upperlim_nan[numb]
        database["ycens"+str(lines_f[i])] = int_new

    fname_dict = base_dir_out+galaxyID+'_base+MeerKAT+P-DE_hexagon_1p5kpc_uplims_ints'+'.ecsv'
    database.write(fname_dict, format='ascii.ecsv',overwrite=True)

    return database
