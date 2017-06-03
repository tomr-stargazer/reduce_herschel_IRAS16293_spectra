"""
This is a script to fit the damn lines.

"""

import glob
import os
from collections import OrderedDict

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u

from fit_gaussian_with_class import fit_line_and_return_raw_output, extract_line_params_from_raw_output
from read_in_fits_spectrum_from_class import load_a_spectrum
from c18o_line_frequencies_and_names import *

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
    elif 1052 < freq < 1121:
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


def baseline_and_fit_c18o_lines():
    for line_freq, line_name in zip(c18o_line_freqs, c18o_line_names):

        band = which_hifi_band(line_freq)

        print(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue

        line_filename = band+"-averaged.hifi"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, smooth_gauss=0.62,
            save=True, output_file_root="C18O"+line_name)

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        print(gaussian_fit_results)
        print("")

    return


def baseline_and_fit_c17o_lines():
    for line_freq, line_name in zip(c17o_line_freqs, c17o_line_names):

        band = which_hifi_band(line_freq)

        print(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue

        line_filename = band+"-averaged.hifi"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, smooth_gauss=0.62,
            save=True, output_file_root="C17O"+line_name)

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        print(gaussian_fit_results)
        print("")

    return


def make_c18o_all_lines_fig():
    fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")
    list_of_files = glob.glob(fit_results_path+"/C18O*.fits")
    list_of_spectra = [x for x in list_of_files if 'spectrum.fits' in x]
    list_of_results = [x for x in list_of_files if 'result.fits' in x]

    fig = plt.figure()

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra, list_of_results)):

        ax = fig.add_subplot(3,4,i+1)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-12, 18)
        ax.set_ylim(-0.25, 2)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "C$^{18}$O "+line_name_nospace.lstrip("C18O")
        ax.text(-10, 1.5, line_name)

    plt.show()

    return fig


def make_c17o_all_lines_fig():
    fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")
    list_of_files = glob.glob(fit_results_path+"/C17O*.fits")
    list_of_spectra = [x for x in list_of_files if 'spectrum.fits' in x]
    list_of_results = [x for x in list_of_files if 'result.fits' in x]

    fig = plt.figure()

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra, list_of_results)):

        ax = fig.add_subplot(3,4,i+1)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-12, 18)
        ax.set_ylim(-0.25, 0.75)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "C$^{17}$O "+line_name_nospace.lstrip("C17O")
        ax.text(-10, 0.6, line_name)

    plt.show()

    return fig


def baseline_spectra_and_compute_fits():
    # ok so the REAL loop is gonna look like this:
    c18o_linefits = OrderedDict()

    # first, fit (and save the results of) all the HCN lines.
    for i, (line_freq, line_name) in enumerate(zip(c18o_line_freqs, c18o_line_names)):

        band = which_hifi_band(line_freq)
        print(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue

        line_filename = band+"-averaged.hifi"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, smooth_gauss=0.62,
            save=True, output_file_root="C18O"+line_name)

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        print(gaussian_fit_results)
        print("")
        c18o_linefits[i] = gaussian_fit_results
        c18o_linefits[i]['line_name'] = line_name

    # Keep those fits around in a list or a dict or something.
    # Second, do that with the h13cn lines. 
    c17o_linefits = OrderedDict()

    for i, (line_freq, line_name) in enumerate(zip(c17o_line_freqs, c17o_line_names)):

        band = which_hifi_band(line_freq)
        print(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue
        if i not in c18o_linefits.keys():
            continue

        line_filename = band+"-averaged.hifi"

        # Use the line center and line width as a soft prior for each fit.
        line_params_string = "0 0 0 {0:.2f} 0 {1:.2f}".format(c18o_linefits[i]['v_cen'], c18o_linefits[i]['v_fwhm'] )
        # if line_name == "J=10 - 9":
        #     line_params_string = "0 0 1 {0:.2f} 1 {1:.2f}".format(hcn_linefits[i]['v_cen'], hcn_linefits[i]['v_fwhm'] )

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, line_params=line_params_string,
            smooth_gauss=0.62, save=True, output_file_root="C17O_"+line_name)

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        print(gaussian_fit_results)
        print("")
        c17o_linefits[i] = gaussian_fit_results
        c17o_linefits[i]['line_name'] = line_name

    # Finally, do that with the h15cn lines.
    # Use the line center and the linewidth as firm priors.
    # hc15n_linefits = OrderedDict()

    # for i, (line_freq, line_name) in enumerate(zip(hc15n_line_freqs, hc15n_line_names)):

    #     band = which_hifi_band(line_freq)
    #     print(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
    #     if band == "None":
    #         continue
    #     if i not in hcn_linefits.keys():
    #         continue

    #     line_filename = band+"-averaged.hifi"

    #     # Use the line center and line width as a hard prior for each fit.
    #     n_lines = 1
    #     line_params_string = "0 0 1 {0:.2f} 1 {1:.2f}".format(h13cn_linefits[i]['v_cen'], h13cn_linefits[i]['v_fwhm'] )
    #     if line_name == "J= 6-5":
    #         custom_window = "-40 -29 -5 15"
    #     elif line_name == 'J= 7-6':
    #         custom_window = "-15 15 20 30"
    #         n_lines = 2
    #         line_params_string = "0 0 1 {0:.2f} 1 {1:.2f}\" \"0 1 0 -4.231 0 5.926".format(h13cn_linefits[i]['v_cen'], h13cn_linefits[i]['v_fwhm'])
    #     elif line_name == 'J= 8-7':
    #         custom_window = "-5 15 23 40"
    #     elif line_name == 'J= 9-8':
    #         custom_window = "-30 -18 -5 15"
    #     elif line_name ==  'J=10-9':
    #         custom_window = "-5 15 21 36"

    #     raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
    #         filename=line_filename, freq=line_freq*1000, line_params=line_params_string, n_lines=n_lines, custom_window=custom_window,
    #         smooth_gauss=0.62, save=True, output_file_root="HC15N_"+line_name)

    #     gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

    #     print(gaussian_fit_results)
    #     print("")
    #     hc15n_linefits[i] = gaussian_fit_results
    #     hc15n_linefits[i]['line_name'] = line_name

    return c18o_linefits, c17o_linefits


def make_h13cn_all_lines_fig():
    fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")
    list_of_files_13 = glob.glob(fit_results_path+"/H13CN*.fits")
    list_of_spectra_13 = [x for x in list_of_files_13 if 'spectrum.fits' in x]
    list_of_results_13 = [x for x in list_of_files_13 if 'result.fits' in x]

    fig = plt.figure()

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra_13, list_of_results_13)):

        ax = fig.add_subplot(3,4,i+1)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-40, 40)
        ax.set_ylim(-0.1, 0.2)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "H13CN "+line_name_nospace.lstrip("H13CN_")
        ax.text(-38, 0.14, line_name)

    plt.show()
    return fig


def make_hc15n_all_lines_fig():
    fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")
    list_of_files_15 = glob.glob(fit_results_path+"/HC15N*.fits")
    list_of_spectra_15 = [x for x in list_of_files_15 if 'spectrum.fits' in x]
    list_of_results_15 = [x for x in list_of_files_15 if 'result.fits' in x]

    fig = plt.figure()

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra_15, list_of_results_15)):

        ax = fig.add_subplot(3,4,i+1)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-40, 40)
        ax.set_ylim(-0.05, 0.1)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "HC15N "+line_name_nospace.lstrip("HC15N_")
        ax.text(-38, 0.07, line_name)

    plt.show()
    return fig

def make_hcn_h13cn_hc15n_figure():
    fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")

    list_of_files = glob.glob(fit_results_path+"/HCN*.fits")
    list_of_spectra = [x for x in list_of_files if 'spectrum.fits' in x]
    list_of_results = [x for x in list_of_files if 'result.fits' in x]

    list_of_files_13 = glob.glob(fit_results_path+"/H13CN*.fits")
    list_of_spectra_13 = [x for x in list_of_files_13 if 'spectrum.fits' in x]
    list_of_results_13 = [x for x in list_of_files_13 if 'result.fits' in x]

    list_of_files_15 = glob.glob(fit_results_path+"/HC15N*.fits")
    list_of_spectra_15 = [x for x in list_of_files_15 if 'spectrum.fits' in x]
    list_of_results_15 = [x for x in list_of_files_15 if 'result.fits' in x]

    fig = plt.figure(figsize=(8.6, 4.8))

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra, list_of_results)):

        if i > 4:
            break

        ax = fig.add_subplot(3,5,i+1)
        ax.tick_params(axis='both', labelsize=7)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-30, 30)
        ax.set_ylim(-0.1, 1)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "HCN "+line_name_nospace.lstrip("HCN")
        ax.text(-25, 0.8, line_name, fontsize=9, bbox={'facecolor':'white', 'alpha':0.5, 'edgecolor':'none'})

        ax.tick_params(axis='x', labelbottom='off')
        if i>0:
            ax.tick_params(axis='y', labelleft='off')

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra_13, list_of_results_13)):

        if i > 4:
            break

        ax = fig.add_subplot(3,5,i+6)
        ax.tick_params(axis='both', labelsize=7)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-30, 30)
        ax.set_ylim(-0.09, 0.15)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "H$^{13}$CN "+line_name_nospace.lstrip("H13CN_")
        ax.text(-25, 0.11, line_name, fontsize=9, bbox={'facecolor':'white', 'alpha':0.5, 'edgecolor':'none'})

        ax.tick_params(axis='x', labelbottom='off')
        if i>0:
            ax.tick_params(axis='y', labelleft='off')        

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra_15, list_of_results_15)):

        if i > 4:
            break

        ax = fig.add_subplot(3,5,i+11)
        ax.tick_params(axis='both', labelsize=7)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-30, 30)
        ax.set_ylim(-0.05, 0.1)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "HC$^{15}$N "+line_name_nospace.lstrip("HC15N_")
        ax.text(-25, 0.07, line_name, fontsize=9, bbox={'facecolor':'white', 'alpha':0.5, 'edgecolor':'none'})

        if i==1:
            xs = np.arange(-30, 30, 0.5)
            a = 0.03666304656927235
            b = 3.67
            c = 6.45 / 2.35482
            ys = a * np.exp( - (xs-b)**2 / (2*c**2))

            ax.plot(xs, ys, 'C0', lw=1)

        if i>0:
            ax.tick_params(axis='y', labelleft='off')

    plt.show()
    return fig

if __name__ == "__main__":

    if False:

        fit_tuple = baseline_spectra_and_compute_fits()

        paper_path = os.path.expanduser("~/Documents/Academia/Articles/Nitrogen_Paper/")
        fig_filename = "in_progress_graphics/hifi_hcn_lines.pdf"
        fig_fullpath = os.path.join(paper_path, fig_filename)

        fig = make_hcn_h13cn_hc15n_figure()
        fig.savefig(fig_fullpath, bbox_inches='tight')