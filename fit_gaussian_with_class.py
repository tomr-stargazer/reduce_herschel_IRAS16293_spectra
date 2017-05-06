"""
This is a proof-of-concept.

I wanna be able to take some data and fit


"""

import os
import subprocess
from collections import OrderedDict
import numpy as np

root_directory_of_data = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/")
data_location = os.path.join(root_directory_of_data, "Partially_Reduced_Spectra")
fit_results_path = os.path.join(root_directory_of_data, "Fit_results")

def fit_line_and_return_raw_output(
    filename = "3b-averaged.hifi", freq=863071, 
    line_params = "0 0 0 3.5 0 7", smooth_channels=None,
    save=False, path=fit_results_path, output_file_root=None):

    # fits write fitswrite_test_result.fits /mode spectrum
    # result
    # fits write blah

    if save:
        # construct the two filenames
        if output_file_root is None:
            output_file_root = filename.split('.')[0]

        spectrum_filename = output_file_root+"_spectrum.fits"
        result_filename = output_file_root+"_result.fits"

        spectrum_fullpath = os.path.join(path, spectrum_filename)
        result_fullpath = os.path.join(path, result_filename)

        # clobber the existing files.
        if os.path.exists(spectrum_fullpath):
            os.remove(spectrum_fullpath)
        if os.path.exists(result_fullpath):
            os.remove(result_fullpath)

        save_script_snippet = (
            "fits write \"{0}\" /mode spectrum\n"
            "result\n"
            "fits write \"{1}\" /mode spectrum\n"
            "".format(spectrum_fullpath, result_fullpath)
            )
        # print (save_script_snippet)
    else:
        save_script_snippet = ""

    if smooth_channels is not None and smooth_channels > 1.5:
        smooth_snippet = "smooth box {0}\n".format(smooth_channels)
    else:
        smooth_snippet = ""

    original_directory = os.getcwd()

    try:
        os.chdir(data_location)

        class_script_text = (
            "file in {0}\n" 
            "find\n"
            "get first\n"
            "set mode x auto\n"
            "set unit v v\n"
            "mod freq {1}\n" 
            "set mode x -40 +40\n"
            "set win -5 15\n"
            "base\n"
            "{4}"
            "method gauss\n"
            "lines 1 \"{2}\"\n" 
            "mini\n"
            "{3}"
            "sic\\examine R%HEAD%GAU%RESULT\n"
            "sic\\examine R%HEAD%GAU%ERROR\n"
            "exit\n".format(
                filename, freq, line_params, save_script_snippet, smooth_snippet))

        path_to_script_with_filename = os.path.join(os.curdir, "first_try_linefit.class")
        target = open(path_to_script_with_filename, 'w')
        target.write(class_script_text)
        target.close()

        p = subprocess.Popen(['class', '-nw', "@first_try_linefit"], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()

    except Exception as e:
        os.chdir(original_directory)
        raise e
    finally:
        os.chdir(original_directory)

    return output, error


def extract_line_params_from_raw_output(raw_output):
    """ Will depend on the specifics of what we're passing to class. """

    output_lines = raw_output.decode('utf-8').split('\n')

    # This may change if we add anything to the end of the CLASS script above!
    results_line = output_lines[-8]
    uncertainties_line = output_lines[-4]

    area, v_cen, v_fwhm, a, b = [float(x) for x in results_line.split()]
    e_area, e_v_cen, e_v_fwhm, a, b = [float(x) for x in uncertainties_line.split()]

    results_dict = OrderedDict()
    results_dict['area'] = area
    results_dict['e_area'] = e_area
    results_dict['v_cen'] = v_cen
    results_dict['e_v_cen'] = e_v_cen
    results_dict['v_fwhm'] = v_fwhm
    results_dict['e_v_fwhm'] = e_v_fwhm

    t_peak = area * 2.35482 / ( (2*np.pi)**(1/2) * v_fwhm)

    results_dict['t_peak'] = t_peak

    return results_dict


if __name__ == '__main__':

    filename = "3b-averaged.hifi"
    freq=863071
    line_params = "0 0 0 3.5 0 7"

    output, error = fit_line_and_return_raw_output(filename, freq, line_params)
    results = extract_line_params_from_raw_output(output)

    print(results)


