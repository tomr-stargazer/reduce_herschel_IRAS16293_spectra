from __future__ import division

import os

list_of_bands = ["1a", "1b", "2a", "2b", "3a", "3b", "4a", "4b", "5a", "6a", "6b", "7a"]

root_directory_of_data = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/")
level_2_5_data = os.path.join(root_directory_of_data, "level_2_5_all_bands")

for band in list_of_bands:

	data_location = os.path.join(level_2_5_data, band, "level2_5/myDecon/")
	filename_horizontal = os.path.join(data_location, )

	print band, "data exists: ", os.path.exists()

	# if they don't exist, we gotta make them exist
	# so the data REALLY come from a file that looks like this:

	# /Users/tsrice/Documents/Data/Herschel_Science_Archive/IRAS16293/level_2_5_all_bands/7a/level2_5/myDecon/myDecon_WBS-V/hhifiwbsvssb1342191762_25ssv20_1461134514602.fits.gz
	# and I think we're supposed to use class to convert that FITS file into 


convert_fits_scriptname = "convert_fits_to_hifi.class"


def write_class_script_to_convert_FITS_to_HIFI(location, fits_filename, hifi_filename):

	script_text = 	"file out {0} mul"
	                "fits read {1}".format(hifi_filename, fits_filename)

	# write this text to the location:
	path_to_script_with_filename = os.path.join(location, convert_fits_scriptname)
	target = open(path_to_script_with_filename, 'w')
	target.write(script_text)
	target.close()

	return path_to_script_with_filename


def call_class_script_to_convert_FITS_to_HIFI(location):

	original_directory = os.getcwd()

	try: 

		from subprocess import call
		os.chdir(location)

		call(["class", "@"+convert_fits_scriptname])

	except Exception, e:
		raise e
	finally:
		os.chdir(original_directory)

	return