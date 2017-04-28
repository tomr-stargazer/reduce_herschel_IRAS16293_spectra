""""
"This is my proof-of-concept that I could, in principle, do the following even ONCE\n:

- find a two .hifi files from different polarizations
- make a simple CLASS script that combines & averages them
- run that CLASS script.

"""

from __future__ import division
import os
import shutil 
from subprocess import call

average_polarizations_scriptname = "average_polarizations.class"


def write_class_script_to_average_polarizations(location, band_name):

    if not len(band_name) == 2:
        raise ValueError("`band_name` must be 2 characters (received `{0}`)".format(band_name))

    script_text = ("file in {0}-horizontal.hifi\n"
                   "find\n"
                   "get first\n\n"
                   
                   "file in {0}-vertical.hifi\n"
                   "find\n"
                   "get first\n\n"
                   
                   "set weight equal\n"
                   "accumulate /resample /nocheck position\n"
                   "multiply 0.5"
                   "\n"
                   
                   "file out {0}-averaged.hifi mul /overwrite\n"
                   "write\n"
                   "exit".format(band_name))

    # write this text to the location:
    path_to_script_with_filename = os.path.join(location, average_polarizations_scriptname)
    target = open(path_to_script_with_filename, 'w')
    target.write(script_text)
    target.close()

    return path_to_script_with_filename

destination_folder = "/Users/tsrice/Documents/Data/Herschel_Science_Archive/IRAS16293/level_2_5_all_bands/6b/level2_5/myDecon"

band_name = "6b"

def average_polarizations(destination_folder, band_name, clobber=False):

    original_directory = os.getcwd()

    try: 
        os.chdir(destination_folder)

        horizontal_folder = os.path.join(destination_folder, "myDecon_WBS-H")
        vertical_folder = os.path.join(destination_folder, "myDecon_WBS-V")        

        horizontal_filename = band_name+"-horizontal.hifi"
        vertical_filename = band_name+"-vertical.hifi"

        # first let's copy the .hifi files here, if needed
        if not os.path.exists(horizontal_filename):

            source_file_path_H = os.path.join(horizontal_folder, horizontal_filename)
            shutil.copy2(source_file_path_H, os.curdir)

        if not os.path.exists(vertical_filename):

            source_file_path_V = os.path.join(vertical_folder, vertical_filename)
            shutil.copy2(source_file_path_V, os.curdir)

        # then make a script and execute it to generate the averaged .hifi file
        output_filename = "{0}-averaged.hifi".format(band_name)
        output_fullpath = os.path.join(destination_folder, output_filename)

        if not os.path.exists(output_fullpath) or clobber:

            average_pols_script_filename = write_class_script_to_average_polarizations(destination_folder, band_name)
            call(["class", "-nw", "@"+average_pols_script_filename])

    except Exception as e:
        os.chdir(original_directory)
        raise e
    finally:
        os.chdir(original_directory)

    return output_fullpath

if __name__ == "__main__":

    destination_folder = "/Users/tsrice/Documents/Data/Herschel_Science_Archive/IRAS16293/level_2_5_all_bands/6a/level2_5/myDecon"
    band_name = "6a"

    average_polarizations(destination_folder, band_name, clobber=True)
