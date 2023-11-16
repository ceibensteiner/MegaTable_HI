from uplims_functions import *

p='/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/scripts/0-data_base/1b-incl_new_P_DE/'
base_out='/Users/cosimaeibensteiner/Desktop/home/PhD/1-Project/MEERKAT/scripts/0-data_base/2-MegaTables_with_sn_and_ups/'

file = 'NGC1512_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC1512'
n1512 = get_new_base(p,file,galaxyID,base_out)

file = 'NGC4535_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC4535'
n4535 = get_new_base(p,file,galaxyID,base_out)

file = 'NGC7496_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC7496'
n7496 = get_new_base(p,file,galaxyID,base_out)

file = 'IC1954_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'IC1954'
i1954 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC1566_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC1566'
n1566 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC1672_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC1672'
n1672 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC3511_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC3511'
n3511 = get_new_base_mhongoose(p,file,galaxyID,base_out)

file = 'NGC5068_base+MeerKAT_hexagon_1p5kpc.ecsv'
galaxyID = 'NGC5068'
n5086 = get_new_base_mhongoose(p,file,galaxyID,base_out)
