"""
This is a script to fit the damn lines.

"""

import glob
import os

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u

from fit_gaussian_with_class import fit_line_and_return_raw_output, extract_line_params_from_raw_output
from read_in_fits_spectrum_from_class import load_a_spectrum

class EmissionLine:
    def __init__(self, freq, molecule_name, line_name):

        self.freq = freq
        self.molecule_name = molecule_name
        self.line_name = line_name

    def line_identification(self):

        return self.molecule_name+" "+self.line_name


def which_hifi_band(freq):

    if 480 < freq < 560:
        return "1a"
    elif 554 < freq < 636:
        return "1b"
    elif 626 < freq < 726:
        return "2a"
    elif 714< freq < 801:
        return "2b"
    elif 799 < freq < 860:
        return "3a"
    elif 858 < freq < 961:
        return "3b"
    elif 949 < freq < 1061:
        return "4a"
    elif 1052 < freq < 1122:
        return "4b"
    elif 1108 < freq < 1238:
        return "5a"
    elif 1481 < freq < 1511:
        return "6a"
    elif 1575 < freq < 1700:
        return "6b"
    elif 1703 < freq < 1799:
        return "7a"
    else:
        return "None"

hcn_line_freqs = [
    531.71639, 
    620.30410, 
    708.87721, 
    797.43366, 
    885.97141, 
    974.48840, 
    1062.98260, 
    1151.45200, 
    1239.89459, 
    1328.30837, 
    1416.69138, 
    1505.04164, 
    1593.35722, 
    1681.63621, 
    1769.87669, 
    1858.07678, 
    1946.23464]

hcn_line_names = [
    "J= 6- 5",  
    "J= 7- 6",  
    "J= 8- 7",  
    "J= 9- 8",  
    "J=10- 9",  
    "J=11-10",  
    "J=12-11",  
    "J=13-12",  
    "J=14-13",  
    "J=15-14",  
    "J=16-15",  
    "J=17-16",  
    "J=18-17",  
    "J=19-18",  
    "J=20-19",  
    "J=21-20",  
    "J=22-21"]

if False:
    for line_freq, line_name in zip(hcn_line_freqs, hcn_line_names):

        band = which_hifi_band(line_freq)

        print(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue

        line_filename = band+"-averaged.hifi"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(filename=line_filename, freq=line_freq*1000, save=True, output_file_root="HCN"+line_name)

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        print(gaussian_fit_results)
        print("")

fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")
list_of_files = glob.glob(fit_results_path+"/HCN*.fits")
list_of_spectra = [x for x in list_of_files if 'spectrum.fits' in x]
list_of_results = [x for x in list_of_files if 'result.fits' in x]

fig = plt.figure()

for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra, list_of_results)):

    ax = fig.add_subplot(3,4,i+1)

    spectrum_tuple = load_a_spectrum(spectrum_fname)
    result_tuple = load_a_spectrum(result_fname)

    ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1.5)
    ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

    ax.set_xlim(-40, 40)
    ax.set_ylim(-0.5, 1)

    line_name = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
    ax.text(-30, 0.8, line_name)

plt.show()