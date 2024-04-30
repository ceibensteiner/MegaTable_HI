import sys
import numpy as np
from astropy import units as u, constants as const
from astropy.io import fits
from astropy.table import Table, QTable
from matplotlib.pyplot import cm

from astropy.table import Table
from astropy import constants as const
sys.path.append('/Users/ceibenst/Documents/GitHub/MegaTable_HI/')
from mega_table_functions import *

sys.path.append('/Users/ceibenst/Documents/GitHub/MegaTable/')
from mega_table import RadialMegaTable, TessellMegaTable

#path where databases are with the MeerKAT HI included
p = '/Users/ceibenst/Desktop/home/1-Projects/MEERKAT/scripts/0-data_base/4-Mega-Table_v4p0/'
#path where we store the new databases that will include new calculations of P_DE_new
path_out = '/Users/ceibenst/Desktop/home/1-Projects/MEERKAT/scripts/0-data_base/4-Mega-Table_v4p0/'

n1512 = Table.read(p+'NGC1512_base+MeerKAT_hexagon_1p5kpc.ecsv')
n4535 = Table.read(p+'NGC4535_base+MeerKAT_hexagon_1p5kpc.ecsv')
n7496 = Table.read(p+'NGC7496_base+MeerKAT_hexagon_1p5kpc.ecsv') #kein P_DE
i1954 = Table.read(p+'IC1954_base+MeerKAT_hexagon_1p5kpc.ecsv')  #kein P_DE
n1672 = Table.read(p+'NGC1672_base+MeerKAT_hexagon_1p5kpc.ecsv') #kein P_DE
n1566 = Table.read(p+'NGC1566_base+MeerKAT_hexagon_1p5kpc.ecsv') #kein P_DE
n3511 = Table.read(p+'NGC3511_base+MeerKAT_hexagon_1p5kpc.ecsv')
n5068 = Table.read(p+'NGC5068_base+MeerKAT_hexagon_1p5kpc.ecsv')

hi_n1512 = n1512["Sigma_atom_MeerKAT"]
hi_unc_n1512 = n1512["e_Sigma_atom_MeerKAT"]

co_n1512 = n1512["Sigma_mol"]
co_unc_n1512 = n1512["e_Sigma_mol"]

#------------------------------------
hi_n4535     = n4535["Sigma_atom_MeerKAT"]
hi_unc_n4535 = n4535["e_Sigma_atom_MeerKAT"]

co_n4535     = n4535["Sigma_mol"]
co_unc_n4535 = n4535["e_Sigma_mol"]

#------------------------------------
hi_n7496    = n7496["Sigma_atom_MeerKAT"]
hi_unc_n7496 = n7496["e_Sigma_atom_MeerKAT"]

co_n7496     = n7496["Sigma_mol"]
co_unc_n7496 = n7496["e_Sigma_mol"]

#------------------------------------
hi_i1954     = i1954["Sigma_atom_MHONGOOSE"]
hi_unc_i1954 = i1954["e_Sigma_atom_MHONGOOSE"]

co_i1954     = i1954["Sigma_mol"]
co_unc_i1954 = i1954["e_Sigma_mol"]

#------------------------------------
hi_n1566     = n1566["Sigma_atom_MHONGOOSE"]
hi_unc_n1566 = n1566["e_Sigma_atom_MHONGOOSE"]

co_n1566     = n1566["Sigma_mol"]
co_unc_n1566 = n1566["e_Sigma_mol"]

#------------------------------------
hi_n1672     = n1672["Sigma_atom_MHONGOOSE"]
hi_unc_n1672 = n1672["e_Sigma_atom_MHONGOOSE"]

co_n1672     = n1672["Sigma_mol"]
co_unc_n1672 = n1672["e_Sigma_mol"]

#------------------------------------
hi_n3511     = n3511["Sigma_atom_MHONGOOSE"]
hi_unc_n3511 = n3511["e_Sigma_atom_MHONGOOSE"]

co_n3511     = n3511["Sigma_mol"]
co_unc_n3511 = n3511["e_Sigma_mol"]

#------------------------------------
hi_n5068     = n5068["Sigma_atom_MHONGOOSE"]
hi_unc_n5068 = n5068["e_Sigma_atom_MHONGOOSE"]

co_n5068     = n5068["Sigma_mol"]
co_unc_n5068 = n5068["e_Sigma_mol"]


mstar_1512 = n1512['Sigma_star']
mstar_4535 = n4535['Sigma_star']
mstar_7496 = n7496['Sigma_star']
mstar_1954 = i1954['Sigma_star']
mstar_1566 = n1566['Sigma_star']
mstar_1672 = n1672['Sigma_star']
mstar_3511 = n3511['Sigma_star']
mstar_5068 = n5068['Sigma_star']


mstar_e_1512 = n1512['e_Sigma_star']
mstar_e_4535 = n4535['e_Sigma_star']
mstar_e_7496 = n7496['e_Sigma_star']
mstar_e_1954 = i1954['e_Sigma_star']
mstar_e_1566 = n1566['e_Sigma_star']
mstar_e_1672 = n1672['e_Sigma_star']
mstar_e_3511 = n3511['e_Sigma_star']
mstar_e_5068 = n5068['e_Sigma_star']

scale_length = [6.2, 3.8, 1.5, 3.9, 5.8, 2.4, 1.3, 1.5] #from Leroy2021 - sample table

sigma_hi_z = 11.0 * u.Unit('km s-1')

#Sigma_star, radial_scale_lenght, Sigma_mol, Sigma_atom, sigma_hi_z, e_Sigma_mol, e_Sigma_atom, e_Sigma_star
pde_n1512, pde_e_n1512 = pde(mstar_1512,scale_length[0],co_n1512,hi_n1512,sigma_hi_z, co_unc_n1512, hi_unc_n1512, mstar_e_1512)
pde_n4535, pde_e_n4535 = pde(mstar_4535,scale_length[1],co_n4535,hi_n4535,sigma_hi_z, co_unc_n4535, hi_unc_n4535, mstar_e_4535)
pde_n7496, pde_e_n7496 = pde(mstar_7496,scale_length[2],co_n7496,hi_n7496,sigma_hi_z, co_unc_n7496, hi_unc_n7496, mstar_e_7496)
pde_n1566, pde_e_n1566 = pde(mstar_1566,scale_length[3],co_n1566,hi_n1566,sigma_hi_z, co_unc_n1566, hi_unc_n1566, mstar_e_1566)
pde_n1672, pde_e_n1672 = pde(mstar_1672,scale_length[4],co_n1672,hi_n1672,sigma_hi_z, co_unc_n1672, hi_unc_n1672, mstar_e_1672)
pde_n3511, pde_e_n3511 = pde(mstar_3511,scale_length[5],co_n3511,hi_n3511,sigma_hi_z, co_unc_n3511, hi_unc_n3511, mstar_e_3511)
pde_n5068, pde_e_n5068 = pde(mstar_5068,scale_length[6],co_n5068,hi_n5068,sigma_hi_z, co_unc_n5068, hi_unc_n5068, mstar_e_5068)
pde_i1954, pde_e_i1954 = pde(mstar_1954,scale_length[7],co_i1954,hi_i1954,sigma_hi_z, co_unc_i1954, hi_unc_i1954, mstar_e_1954)


pde_new = [pde_i1954, pde_n1566, pde_n1672, pde_n3511, pde_n5068, pde_n1512, pde_n4535,pde_n7496]
pde_e_new = [pde_e_i1954, pde_e_n1566, pde_e_n1672, pde_e_n3511, pde_e_n5068, pde_e_n1512, pde_e_n4535,pde_e_n7496]

#we only have EW maps for the PHANGS-MeerKAT targests atm; MHONGOOSE are on their way? (asked Erwin)
sigma_hi_z_n1512 = n1512['EW_21cm_MeerKAT']
sigma_hi_z_n4535 = n4535['EW_21cm_MeerKAT']
sigma_hi_z_n7496 = n7496['EW_21cm_MeerKAT']

sigma_co_z_n1512 = n1512['EW_CO']
sigma_co_z_n4535 = n4535['EW_CO']
sigma_co_z_n7496 = n7496['EW_CO']
sigma_co_z_i1954 = i1954['EW_CO']
sigma_co_z_n1566 = n1566['EW_CO']
sigma_co_z_n1672 = n1672['EW_CO']
sigma_co_z_n3511 = n3511['EW_CO']
sigma_co_z_n5068 = n5068['EW_CO']


#P_DE calculation where we vary both (i) sigma mol and (ii) sigma atom !! // atm only for PHANGS-MeerKAT
pde_varying_n1512, pde_e_varying_n1512 = pde_varying(mstar_1512,scale_length[0],co_n1512,hi_n1512,sigma_hi_z_n1512,sigma_co_z_n1512,co_unc_n1512,hi_unc_n1512,mstar_e_1512)
pde_varying_n4535, pde_e_varying_n4535 = pde_varying(mstar_4535,scale_length[1],co_n4535,hi_n4535,sigma_hi_z_n4535,sigma_co_z_n4535,co_unc_n4535,hi_unc_n4535,mstar_e_4535)
pde_varying_n7496, pde_e_varying_n7496 = pde_varying(mstar_7496,scale_length[2],co_n7496,hi_n7496,sigma_hi_z_n7496,sigma_co_z_n7496,co_unc_n7496,hi_unc_n7496,mstar_e_7496)

pde_new_varying = [pde_varying_n1512,pde_varying_n4535,pde_varying_n7496]
e_pde_new_varying = [pde_e_varying_n1512,pde_e_varying_n4535,pde_e_varying_n7496]

#P_DE calculations where we vary only (i) sigma mol and fix sig atom to 10 km/s. 
pde_varying_mol_n1512, pde_e_varying_mol_n1512 = pde_varying(mstar_1512,scale_length[0],co_n1512,hi_n1512,sigma_hi_z,sigma_co_z_n1512,co_unc_n1512,hi_unc_n1512,mstar_e_1512)
pde_varying_mol_n4535, pde_e_varying_mol_n4535 = pde_varying(mstar_4535,scale_length[1],co_n4535,hi_n4535,sigma_hi_z,sigma_co_z_n4535,co_unc_n4535,hi_unc_n4535,mstar_e_4535)
pde_varying_mol_n7496, pde_e_varying_mol_n7496 = pde_varying(mstar_7496,scale_length[2],co_n7496,hi_n7496,sigma_hi_z,sigma_co_z_n7496,co_unc_n7496,hi_unc_n7496,mstar_e_7496)
pde_varying_mol_n1566, pde_e_varying_mol_n1566 = pde_varying(mstar_1566,scale_length[3],co_n1566,hi_n1566,sigma_hi_z,sigma_co_z_n1566,co_unc_n1566,hi_unc_n1566,mstar_e_1566)
pde_varying_mol_n1672, pde_e_varying_mol_n1672 = pde_varying(mstar_1672,scale_length[4],co_n1672,hi_n1672,sigma_hi_z,sigma_co_z_n1672,co_unc_n1672,hi_unc_n1672,mstar_e_1672)
pde_varying_mol_n3511, pde_e_varying_mol_n3511 = pde_varying(mstar_3511,scale_length[5],co_n3511,hi_n3511,sigma_hi_z,sigma_co_z_n3511,co_unc_n3511,hi_unc_n3511,mstar_e_3511)
pde_varying_mol_n5068, pde_e_varying_mol_n5068 = pde_varying(mstar_5068,scale_length[6],co_n5068,hi_n5068,sigma_hi_z,sigma_co_z_n5068,co_unc_n5068,hi_unc_n5068,mstar_e_5068)
pde_varying_mol_i1954, pde_e_varying_mol_i1954 = pde_varying(mstar_1954,scale_length[7],co_i1954,hi_i1954,sigma_hi_z,sigma_co_z_i1954,co_unc_i1954,hi_unc_i1954,mstar_e_1954)

pde_new_varying_mol = [pde_varying_mol_n1512, pde_varying_mol_n4535, pde_varying_mol_n7496]
e_pde_new_varying_mol = [pde_e_varying_mol_n1512, pde_e_varying_mol_n4535, pde_e_varying_mol_n7496]

galaxy = ['IC1954','NGC1566','NGC1672','NGC3511','NGC5068','NGC1512', 'NGC4535', 'NGC7496']
galaxy_pde_var = ['NGC1512', 'NGC4535', 'NGC7496']

pde_new_varying_mol_=[pde_varying_mol_i1954,
                      pde_varying_mol_n1566,
                      pde_varying_mol_n1672,
                      pde_varying_mol_n3511,
                      pde_varying_mol_n5068,
                      pde_varying_mol_n1512,
                      pde_varying_mol_n4535,
                      pde_varying_mol_n7496]

e_pde_new_varying_mol_=[pde_e_varying_mol_i1954,
                      pde_e_varying_mol_n1566,
                      pde_e_varying_mol_n1672,
                      pde_e_varying_mol_n3511,
                      pde_e_varying_mol_n5068,
                      pde_e_varying_mol_n1512,
                      pde_e_varying_mol_n4535,
                      pde_e_varying_mol_n7496]


bases = [i1954,n1566,n1672,n3511,n5068,n1512,n4535,n7496]


for i in range(len(galaxy)):

    t = TessellMegaTable.read(p+galaxy[i]+'_base+MeerKAT_hexagon_1p5kpc.ecsv')

    colname='P_DE_fixed'
    unit='K cm-3'
    colname_e='e_P_DE_fixed'
    unit_e='K cm-3'

    t[colname] = pde_new[i].to(unit)
    t[colname_e] = pde_e_new[i].to(unit_e)

    colname_='P_DE_varying_mol'
    colname_e_='e_P_DE_varying_mol'

    t[colname_] = pde_new_varying_mol_[i].to(unit)
    t[colname_e_] = e_pde_new_varying_mol_[i].to(unit_e)

    #if i > 5:
    #    for k in range(len(galaxy_pde_var)):
    #        _colname='P_DE_varying'
    #        unit='K cm-3'
    #        _colname_e='e_P_DE_varying'
    #        unit_e='K cm-3'

    #        t[_colname] = pde_new_varying[k].to(unit)
    #        t[_colname_e] = e_pde_new_varying[k].to(unit_e)

    #        colname_='P_DE_varying_mol'
    #        colname_e_='e_P_DE_varying_mol'

    #        t[colname_] = pde_new_varying_mol[k].to(unit)
    #        t[colname_e_] = e_pde_new_varying_mol[k].to(unit_e)

    t.write(path_out+galaxy[i]+'_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)


galaxy_pde_var = ['NGC1512', 'NGC4535', 'NGC7496']

for j in range(len(galaxy_pde_var)):

    a = TessellMegaTable.read(p+galaxy_pde_var[j]+'_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv')

    _colname='P_DE_varying'
    unit='K cm-3'
    _colname_e='e_P_DE_varying'
    unit_e='K cm-3'

    a[_colname] = pde_new_varying[j].to(unit)
    a[_colname_e] = e_pde_new_varying[j].to(unit_e)

    #colname_='P_DE_varying_mol'
    #colname_e_='e_P_DE_varying_mol'

    #a[colname_] = pde_new_varying_mol[j].to(unit)
    #a[colname_e_] = e_pde_new_varying_mol[j].to(unit_e)


    a.write(path_out+galaxy_pde_var[j]+'_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv', add_timestamp=True, delimiter=',', overwrite=True)


#p='/Users/ceibenst/Desktop/home/1-Projects/MEERKAT/scripts/0-data_base/4-Mega-Table_v4p0/'
#test = Table.read(p+'NGC1512_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv')
#print(test['P_DE_varying'])
#print(test['P_DE_fixed'])
