"""
This is my proof-of-concept that I could, in principle, do the following even ONCE:

- find a specific .fits.gz file
- gunzip it
- make a simple CLASS script that turns it into a .hifi format
- run that CLASS script.

"""

from __future__ import division
import os
import glob
from subprocess import call

convert_fits_scriptname = "convert_fits_to_hifi.class"


def write_class_script_to_convert_FITS_to_HIFI(location, fits_filename, hifi_filename, clobber=False):

    if clobber:
        clobber_string = " /OVERWRITE"
    else:
        clobber_string = ""

    script_text =  ("file out {0} mul{2}\n"
                    "fits read {1}\n"
                    "exit".format(hifi_filename, fits_filename, clobber_string) )

    # write this text to the location:
    path_to_script_with_filename = os.path.join(location, convert_fits_scriptname)
    target = open(path_to_script_with_filename, 'w')
    target.write(script_text)
    target.close()

    return path_to_script_with_filename


def convert_FITS_to_HIFI(destination_folder, output_filename, clobber=False):

    original_directory = os.getcwd()

    # this try/except framework is important when changing directories, so we don't get stranded upon error
    try: 
        os.chdir(destination_folder)

        # glob.glob() is like calling "ls" at the command line
        # assume there's only gonna be one .fits file here if any
        if len(glob.glob("*.fits")) == 0:

            # grab the shortest filename among .fits.gz (this helps avoid any backup files)
            fits_gz_filename = sorted(glob.glob("*.fits.gz"), key=len)[0]

            call(["gunzip", fits_gz_filename])

        fits_filename = glob.glob("*.fits")[0]

        # if our desired .hifi file doesn't exist...
        output_fullpath = os.path.join(destination_folder, output_filename)
        if not os.path.exists(output_fullpath) or clobber:
            # ...then make a script and execute it to generate the .hifi file        

            convert_fits_script_filename = write_class_script_to_convert_FITS_to_HIFI(destination_folder, fits_filename, output_filename, clobber=clobber)

            call(["class", "-nw", "@"+convert_fits_script_filename])

    except Exception as e:
        os.chdir(original_directory)
        raise e
    finally:
        os.chdir(original_directory)

    return output_fullpath


if __name__ == "__main__":

    destination_folder = "/Users/tsrice/Documents/Data/Herschel_Science_Archive/IRAS16293/level_2_5_all_bands/7a/level2_5/myDecon/myDecon_WBS-V"
    # fits_gz_filename = "hhifiwbshssb1342191794_25ssv20_1461136819393.fits.gz"
    # fits_filename = fits_gz_filename.rstrip(".gz")
    # fits_gz_fullpath = os.path.join(destination_folder, fits_gz_filename)

    output_filename = "7a-vertical.hifi"
    # output_fullpath = os.path.join(destination_folder, output_filename)

    convert_FITS_to_HIFI(destination_folder, output_filename, clobber=True)

