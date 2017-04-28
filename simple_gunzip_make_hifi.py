"""
This is my proof-of-concept that I could, in principle, do the following even ONCE:

- find a specific .fits.gz file
- gunzip it
- make a simple CLASS script that turns it into a .hifi format
- run that CLASS script.

"""

from __future__ import division
import os
from subprocess import call

convert_fits_scriptname = "convert_fits_to_hifi.class"


def write_class_script_to_convert_FITS_to_HIFI(location, fits_filename, hifi_filename):

    script_text =  ("file out {0} mul\n"
                    "fits read {1}\n"
                    "exit".format(hifi_filename, fits_filename) )

    # write this text to the location:
    path_to_script_with_filename = os.path.join(location, convert_fits_scriptname)
    target = open(path_to_script_with_filename, 'w')
    target.write(script_text)
    target.close()

    return path_to_script_with_filename


destination_folder = "/Users/tsrice/Documents/Data/Herschel_Science_Archive/IRAS16293/level_2_5_all_bands/6b/level2_5/myDecon/myDecon_WBS-H"
fits_gz_filename = "hhifiwbshssb1342191794_25ssv20_1461136819393.fits.gz"
fits_filename = fits_gz_filename.rstrip(".gz")
fits_gz_fullpath = os.path.join(destination_folder, fits_gz_filename)

output_filename = "6b-horizontal.hifi"
output_fullpath = os.path.join(destination_folder, output_filename)


original_directory = os.getcwd()

try: 
    os.chdir(destination_folder)

    # first let's gunzip the file, if needed
    if not os.path.exists(fits_gz_fullpath):
        call(["gunzip", fits_gz_filename])

    # then make a script and execute it to generate the  
    if not os.path.exists(output_fullpath):
        convert_fits_scriptname = write_class_script_to_convert_FITS_to_HIFI(destination_folder, fits_filename, output_filename)
        call(["class", "@"+convert_fits_scriptname])

except Exception as e:
    os.chdir(original_directory)
    raise e
finally:
    os.chdir(original_directory)



