from __future__ import division

import os
import shutil 
from gunzip_make_hifi import convert_FITS_to_HIFI
from combine_and_average import average_polarizations

list_of_bands = ["1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b", "5a", "6a", "6b", "7a"]

root_directory_of_data = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/")
level_2_5_data = os.path.join(root_directory_of_data, "level_2_5_all_bands")
target_location = os.path.join(root_directory_of_data, "Partially_Reduced_Spectra")

for band in list_of_bands:

    data_location = os.path.join(level_2_5_data, band, "level2_5/myDecon/")
    
    data_location_horizontal = os.path.join(data_location, "myDecon_WBS-H")
    data_location_vertical = os.path.join(data_location, "myDecon_WBS-V")

    convert_FITS_to_HIFI(data_location_horizontal, band+"-horizontal.hifi")
    convert_FITS_to_HIFI(data_location_vertical, band+"-vertical.hifi")

    averaged_file_fullpath = average_polarizations(data_location, band, clobber=True)

    shutil.copy2(averaged_file_fullpath, target_location)

