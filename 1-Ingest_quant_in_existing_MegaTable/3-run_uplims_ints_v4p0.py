from uplims_functions_intensities import *

p = '/Users/ceibenst/Desktop/home/1-Projects/MEERKAT/scripts/0-data_base/4-Mega-Table_v4p0/'
base_out='/Users/ceibenst/Desktop/home/1-Projects/MEERKAT/scripts/0-data_base/4-Mega-Table_v4p0/'

file = 'NGC1512_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC1512'
n1512 = get_new_base(p,file,galaxyID,base_out)

file = 'NGC4535_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC4535'
n4535 = get_new_base(p,file,galaxyID,base_out)

file = 'NGC7496_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC7496'
n7496 = get_new_base(p,file,galaxyID,base_out)

file = 'IC1954_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'IC1954'
i1954 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC1566_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC1566'
n1566 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC1672_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC1672'
n1672 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC3511_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC3511'
n3511 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC5068_base+MeerKAT+P-DE_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC5068'
n5086 = get_new_base_mhongoose(p,file,galaxyID,base_out)
