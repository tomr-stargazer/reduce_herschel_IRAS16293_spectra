"""
This is a script to fit the damn lines.

"""

import glob
import os
from collections import OrderedDict
import pdb

import numpy as np
import matplotlib.pyplot as plt
import astropy.constants as c
import astropy.units as u

from fit_gaussian_with_class import fit_line_and_return_raw_output, extract_line_params_from_raw_output
from read_in_fits_spectrum_from_class import load_a_spectrum
from line_frequencies_and_names import *

# currently unused.
class EmissionLine:
    def __init__(self, freq, molecule_name, line_name):

        self.freq = freq
        self.molecule_name = molecule_name
        self.line_name = line_name

    def line_identification(self):

        return self.molecule_name+" "+self.line_name


# should re-implement as a dict or something with a simple accompanying function.
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


def make_linestring(plain_molecule_name, Ju):

    return "{0}_Ju={1:02d}".format(plain_molecule_name, Ju)


def parse_linestring(linestring):

    split_string = linestring.split("_Ju=")

    plain_molecule_name = split_string[0]
    Ju = int(split_string[1])

    return plain_molecule_name, Ju


def baseline_and_fit_HCN_lines():
    for line_freq, Ju in zip(hcn_line_freqs, common_Ju_list):

        line_name = "J={0}-{1}".format(Ju, Ju-1)

        band = which_hifi_band(line_freq)

        print(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue

        line_filename = band+"-averaged.hifi"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, smooth_gauss=0.62,
            save=True, output_file_root="HCN"+line_name)

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        print(gaussian_fit_results)
        print("")

    return

def make_hcn_all_lines_fig():
    fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")
    list_of_files = glob.glob(fit_results_path+"/HCN*.fits")
    list_of_spectra = [x for x in list_of_files if 'spectrum.fits' in x]
    list_of_results = [x for x in list_of_files if 'result.fits' in x]

    fig = plt.figure()

    for i, (spectrum_fname, result_fname) in enumerate(zip(list_of_spectra, list_of_results)):

        ax = fig.add_subplot(3,4,i+1)

        spectrum_tuple = load_a_spectrum(spectrum_fname)
        result_tuple = load_a_spectrum(result_fname)

        ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=1, drawstyle='steps-mid')
        ax.plot(result_tuple[2], result_tuple[0], 'r', lw=0.75)

        ax.set_xlim(-40, 40)
        ax.set_ylim(-0.5, 1)

        line_name_nospace = spectrum_fname.split('/')[-1].rstrip('_spectrum.fits')
        line_name = "HCN "+line_name_nospace.lstrip("HCN")
        ax.text(-30, 0.8, line_name)

    plt.show()

    return fig


def baseline_spectra_and_compute_fits(verbose=False):

    def vprint(*args):
        if verbose:
            print(*args)

    # ok so the REAL loop is gonna look like this:
    hcn_linefits = OrderedDict()

    # first, fit (and save the results of) all the HCN lines.
    for i, (line_freq, Ju) in enumerate(zip(hcn_line_freqs, common_Ju_list)):

        line_name = "J={0}-{1}".format(Ju, Ju-1)

        band = which_hifi_band(line_freq)
        vprint(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue

        line_filename = band+"-averaged.hifi"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, smooth_gauss=0.62,
            save=True, output_file_root=make_linestring("HCN", Ju))

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        vprint(gaussian_fit_results)
        vprint("")
        hcn_linefits[i] = gaussian_fit_results
        hcn_linefits[i]['freq'] = line_freq 
        hcn_linefits[i]['Molecule'] = 'HCN'
        hcn_linefits[i]['Ju'] = Ju
        

    # Keep those fits around in a list or a dict or something.
    # Second, do that with the h13cn lines. 
    h13cn_linefits = OrderedDict()

    for i, (line_freq, Ju) in enumerate(zip(h13cn_line_freqs, common_Ju_list)):

        line_name = "J={0}-{1}".format(Ju, Ju-1)

        band = which_hifi_band(line_freq)
        vprint(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue
        if i not in hcn_linefits.keys():
            continue

        line_filename = band+"-averaged.hifi"

        # Use the line center and line width as a soft prior for each fit.
        line_params_string = "0 0 0 {0:.2f} 0 {1:.2f}".format(hcn_linefits[i]['v_cen'], hcn_linefits[i]['v_fwhm'] )
        if line_name == "J=10 - 9":
            line_params_string = "0 0 1 {0:.2f} 1 {1:.2f}".format(hcn_linefits[i]['v_cen'], hcn_linefits[i]['v_fwhm'] )

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, line_params=line_params_string,
            smooth_gauss=0.62, save=True, output_file_root=make_linestring("H13CN", Ju))

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        vprint(gaussian_fit_results)
        vprint("")
        h13cn_linefits[i] = gaussian_fit_results
        h13cn_linefits[i]['freq'] = line_freq 
        h13cn_linefits[i]['Molecule'] = 'H13CN'
        h13cn_linefits[i]['Ju'] = Ju
        

    # Finally, do that with the h15cn lines.
    # Use the line center and the linewidth as firm priors.
    hc15n_linefits = OrderedDict()

    for i, (line_freq, Ju) in enumerate(zip(hc15n_line_freqs, common_Ju_list)):

        line_name = "J={0}-{1}".format(Ju, Ju-1)

        band = which_hifi_band(line_freq)
        vprint(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue
        if i not in hcn_linefits.keys():
            continue

        line_filename = band+"-averaged.hifi"

        # Use the line center and line width as a hard prior for each fit.
        n_lines = 1
        line_params_string = "0 0 1 {0:.2f} 1 {1:.2f}".format(h13cn_linefits[i]['v_cen'], h13cn_linefits[i]['v_fwhm'] )
        if Ju == 6:
            custom_window = "-40 -29 -5 15"
        elif Ju == 7:
            custom_window = "-15 15 20 30"
            n_lines = 2
            line_params_string = "0 0 1 {0:.2f} 1 {1:.2f}\" \"0 1 0 -4.231 0 5.926".format(h13cn_linefits[i]['v_cen'], h13cn_linefits[i]['v_fwhm'])
        elif Ju == 8:
            custom_window = "-5 15 23 40"
        elif Ju == 9:
            custom_window = "-30 -18 -5 15"
        elif Ju == 10:
            custom_window = "-5 15 21 36"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, line_params=line_params_string, n_lines=n_lines, custom_window=custom_window,
            smooth_gauss=0.62, save=True, output_file_root=make_linestring("HC15N", Ju))

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        if n_lines == 2:
            vprint(str(raw_gaussian_result, 'utf-8'), str(raw_gaussian_error, 'utf-8'))
            if verbose:
                pdb.set_trace()

        vprint(gaussian_fit_results)
        vprint("")
        hc15n_linefits[i] = gaussian_fit_results
        hc15n_linefits[i]['freq'] = line_freq 
        hc15n_linefits[i]['Molecule'] = 'HC15N'
        hc15n_linefits[i]['Ju'] = Ju
        

    return hcn_linefits, h13cn_linefits, hc15n_linefits


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


def plottable_latex_string(plain_molecule_name, Ju):
    """ Assumes the molecule name is simple like HCN or H13CN """
    latex_molname_dict = {}
    latex_molname_dict['HCN'] = 'HCN'
    latex_molname_dict['H13CN'] = r'H$^{13}$CN'
    latex_molname_dict['HC15N'] = r'HC$^{15}$N'

    transition_name = "J$={0}-{1}$".format(Ju, Ju-1)

    plot_string = "{0}  {1}".format(latex_molname_dict[plain_molecule_name], transition_name)
    return plot_string


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

    text_params = dict(fontsize=11, family='serif', bbox={'facecolor':'white', 'alpha':0.5, 'edgecolor':'none'})

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

        mol_name, Ju = parse_linestring(spectrum_fname.split('/')[-1].rstrip('_spectrum.fits'))
        plottable_string = plottable_latex_string(mol_name, Ju)
        ax.text(-25, 0.8, plottable_string, text_params)

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

        mol_name, Ju = parse_linestring(spectrum_fname.split('/')[-1].rstrip('_spectrum.fits'))
        plottable_string = plottable_latex_string(mol_name, Ju)
        ax.text(-25, 0.11, plottable_string, text_params)

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

        mol_name, Ju = parse_linestring(spectrum_fname.split('/')[-1].rstrip('_spectrum.fits'))
        plottable_string = plottable_latex_string(mol_name, Ju)
        ax.text(-25, 0.07, plottable_string, text_params)

        # this is to plot the blue partial component of the blended hc15n 7-6 fit (very tacked-on)
        if i==1:
            xs = np.arange(-30, 30, 0.1)
            a = 0.03666304656927235
            b = 3.67
            c = 6.45 / 2.35482
            ys = a * np.exp( - (xs-b)**2 / (2*c**2))

            ax.plot(xs, ys, 'C0', lw=1)

            a2 = 0.35102
            b2 = -4.813
            c2 = 6.062 / 2.35482
            ys2 = a2 * np.exp( - (xs-b2)**2 / (2*c2**2))

            ax.plot(xs, ys2, 'C1', lw=1)

        if i>0:
            ax.tick_params(axis='y', labelleft='off')

    plt.tight_layout(w_pad=0.5, h_pad=0.3)
    plt.show()
    return fig

if __name__ == "__main__":

    if True:

        fit_tuple = baseline_spectra_and_compute_fits(verbose=True)

        paper_path = os.path.expanduser("~/Documents/Academia/Articles/Nitrogen_Paper/")
        fig_filename = "in_progress_graphics/hifi_hcn_lines.pdf"
        fig_fullpath = os.path.join(paper_path, fig_filename)

        fig = make_hcn_h13cn_hc15n_figure()
        fig.savefig(fig_fullpath, bbox_inches='tight')