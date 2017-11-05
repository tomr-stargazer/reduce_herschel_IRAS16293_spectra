"""
This is a proof-of-concept.

I wanna be able to take some data and fit


"""

import os
import subprocess
from collections import OrderedDict
from functools import reduce
import numpy as np
import pdb

root_directory_of_data = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/")
data_location = os.path.join(root_directory_of_data, "Partially_Reduced_Spectra")
fit_results_path = os.path.join(root_directory_of_data, "Fit_results")

def fit_line_and_return_raw_output(
    filename = "3b-averaged.hifi", freq=863071, n_lines=1,
    line_params = "0 0 0 3.5 0 7", custom_window="-5 15", smooth_gauss=None, smooth_channels=None,
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

    if smooth_gauss is not None:
        smooth_snippet_g = "smooth gauss {0}\n".format(smooth_gauss)
    else:
        smooth_snippet_g = ""

    if smooth_channels is not None and smooth_channels > 1.5:
        smooth_snippet_c = "smooth box {0}\n".format(smooth_channels)        
    else:
        smooth_snippet_c = ""

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
            "extract -100 100\n"
            "set win {5}\n"
            "base\n"
            "{4}"
            "{7}"
            "method gauss\n"
            "lines {6:d} \"{2}\"\n" 
            "mini\n"
            "{3}"
            "sic\\examine R%HEAD%GAU%RESULT\n"
            "sic\\examine R%HEAD%GAU%ERROR\n"
            "exit\n".format(
                filename, freq, line_params, save_script_snippet, smooth_snippet_g, custom_window, n_lines, smooth_snippet_c))

        script_base_filename = 'linefit'
        path_to_script_with_filename = os.path.join(os.curdir, script_base_filename+".class")
        target = open(path_to_script_with_filename, 'w')
        target.write(class_script_text)
        target.close()

        p = subprocess.Popen(['class', '-nw', "@"+script_base_filename], stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        output, error = p.communicate()

    except Exception as e:
        os.chdir(original_directory)
        raise e
    finally:
        os.chdir(original_directory)

    return output, error


def values_list_from_result_lines(list_of_str):
    """ Turns 3 lines of numerical output into a list of 15 values. """

    if len(list_of_str) != 3:
        raise ValueError("`values_list_from_result_lines` takes a list of 3 strings, not {0}.".format(len(list_of_str)))

    joined_lines = reduce(lambda x,y: "".join([x,y]), list_of_str)
    values = [float(x) for x in joined_lines.split()]

    return values


def n_lines_from_values_list(values_list): 

    if len(values_list) != 15:
        raise ValueError("Expecting a list of length 15, not {0}".format(len(values_list)))

    for i in range(15//3):

        if (values_list[0+3*i]==0) and (values_list[1+3*i]==0) and (values_list[2+3*i]==0):
            # i is zero-indexed, 
            # and we want to return one less than the number we just checked, 
            # so this functions as `return (i-1) + 1`
            return i


def extract_line_params_from_raw_output(raw_output):
    """ Will depend on the specifics of what we're passing to class. """

    output_lines = raw_output.decode('utf-8').split('\n')

    # # This may change if we add anything to the end of the CLASS script above!
    # results_line = output_lines[-8]
    # uncertainties_line = output_lines[-4]

    # area, v_cen, v_fwhm, a, b = [float(x) for x in results_line.split()]
    # e_area, e_v_cen, e_v_fwhm, a, b = [float(x) for x in uncertainties_line.split()]

    fit_results = values_list_from_result_lines(output_lines[-8:-5])
    fit_uncertainties = values_list_from_result_lines(output_lines[-4:-1])

    n_lines = n_lines_from_values_list(fit_results)

    # pdb.set_trace()

    results_dict = OrderedDict()
    results_dict['n_lines'] = n_lines
    results_dict['area'] = fit_results[0]
    results_dict['e_area'] = fit_uncertainties[0]
    results_dict['v_cen'] = fit_results[1]
    results_dict['e_v_cen'] = fit_uncertainties[1]
    results_dict['v_fwhm'] = fit_results[2]
    results_dict['e_v_fwhm'] = fit_uncertainties[2]

    t_peak = results_dict['area'] * 2.35482 / ( (2*np.pi)**(1/2) * results_dict['v_fwhm'])
    results_dict['t_peak'] = t_peak

    if n_lines > 1:
        results_dict['area_2'] = fit_results[3]
        results_dict['e_area_2'] = fit_uncertainties[3]
        results_dict['v_cen_2'] = fit_results[4]
        results_dict['e_v_cen_2'] = fit_uncertainties[4]
        results_dict['v_fwhm_2'] = fit_results[5]
        results_dict['e_v_fwhm_2'] = fit_uncertainties[5]

        t_peak_2 = results_dict['area_2'] * 2.35482 / ( (2*np.pi)**(1/2) * results_dict['v_fwhm_2'])
        results_dict['t_peak_2'] = t_peak_2

    return results_dict


if __name__ == '__main__':

    filename = "3b-averaged.hifi"
    freq=863071
    line_params = "0 0 0 3.5 0 7"

    output, error = fit_line_and_return_raw_output(filename=filename, freq=freq, line_params=line_params)
    results = extract_line_params_from_raw_output(output)

    print(results)


