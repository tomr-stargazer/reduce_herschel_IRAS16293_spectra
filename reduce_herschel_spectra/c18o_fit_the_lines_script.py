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
from c18o_line_frequencies_and_names import *
from fit_the_lines_script import make_linestring, parse_linestring, plottable_latex_string

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
        mol_name, Ju = parse_linestring(spectrum_fname.split('/')[-1].rstrip('_spectrum.fits'))
        plottable_string = plottable_latex_string(mol_name, Ju)

        ax.text(-10, 1.5, plottable_string)

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
        mol_name, Ju = parse_linestring(spectrum_fname.split('/')[-1].rstrip('_spectrum.fits'))
        plottable_string = plottable_latex_string(mol_name, Ju)

        ax.text(-10, 0.6, plottable_string)

    plt.show()

    return fig


def co_baseline_spectra_and_compute_fits(verbose=False):

    def vprint(*args):
        if verbose:
            print(*args)

    # ok so the REAL loop is gonna look like this:
    c18o_linefits = OrderedDict()

    # first, fit (and save the results of) all the c18o lines.
    for i, (line_freq, Ju) in enumerate(zip(c18o_line_freqs, common_Ju_list)):

        line_name = "J={0}-{1}".format(Ju, Ju-1)

        band = which_hifi_band(line_freq)
        vprint(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue

        line_filename = band+"-averaged.hifi"

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, smooth_gauss=0.62/2, smooth_channels=2*line_freq/c18o_line_freqs[0],
            save=True, output_file_root=make_linestring("C18O", Ju))


        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        vprint(gaussian_fit_results)
        vprint("")
        c18o_linefits[i] = gaussian_fit_results
        c18o_linefits[i]['freq'] = line_freq 
        c18o_linefits[i]['Molecule'] = 'C18O'
        c18o_linefits[i]['Ju'] = Ju

    # Keep those fits around in a list or a dict or something.
    # Second, do that with the c17o lines. 
    c17o_linefits = OrderedDict()

    for i, (line_freq, Ju) in enumerate(zip(c17o_line_freqs, common_Ju_list)):

        line_name = "J={0}-{1}".format(Ju, Ju-1)

        band = which_hifi_band(line_freq)
        vprint(line_name+" ({0:.3f} GHz): Band ".format(line_freq)+which_hifi_band(line_freq))
        if band == "None":
            continue
        if i not in c18o_linefits.keys():
            continue

        line_filename = band+"-averaged.hifi"

        # Use the line center and line width as a soft prior for each fit.
        n_lines = 1
        line_params_string = "0 0 0 {0:.2f} 0 {1:.2f}".format(c18o_linefits[i]['v_cen'], c18o_linefits[i]['v_fwhm'] )
        if Ju == 7:
            n_lines = 2
            line_params_string = "0 0 1 {0:.2f} 1 {1:.2f}\" \"0 0 1 {2:.2f} 0 {1:.2f}".format(c18o_linefits[i]['v_cen'], c18o_linefits[i]['v_fwhm'], c18o_linefits[i]['v_cen']-1.525)

        raw_gaussian_result, raw_gaussian_error = fit_line_and_return_raw_output(
            filename=line_filename, freq=line_freq*1000, line_params=line_params_string, n_lines=n_lines, 
            smooth_gauss=0.62/2, smooth_channels=2*line_freq/c18o_line_freqs[0], save=True, output_file_root=make_linestring("C17O", Ju))

        gaussian_fit_results = extract_line_params_from_raw_output(raw_gaussian_result)

        if n_lines == 2:
            vprint(str(raw_gaussian_result, 'utf-8'), str(raw_gaussian_error, 'utf-8'))
            if verbose:
                pdb.set_trace() 

        vprint(gaussian_fit_results)
        vprint("")
        c17o_linefits[i] = gaussian_fit_results
        c17o_linefits[i]['freq'] = line_freq 
        c17o_linefits[i]['Molecule'] = 'C17O'
        c17o_linefits[i]['Ju'] = Ju


    return c18o_linefits, c17o_linefits


def make_c18o_c17o_figure(fit_tuple=None):
    fit_results_path = os.path.expanduser("~/Documents/Data/Herschel_Science_Archive/IRAS16293/Fit_results")

    list_of_files_18 = glob.glob(fit_results_path+"/C18O*.fits")
    list_of_files_17 = glob.glob(fit_results_path+"/C17O*.fits")

    fig = plt.figure(figsize=(8.6*1.2, 4.8*0.75))

    text_params = dict(fontsize=11, family='serif', bbox={'facecolor':'white', 'alpha':0.6, 'edgecolor':'none'})

    file_list_list = [list_of_files_18, list_of_files_17]

    ylims = [(-0.15, 2), (-0.15, 0.65)]
    text_heights = [1.65, 0.525]

    for j, file_list in enumerate(file_list_list):

        spectra_list = [x for x in file_list if 'spectrum.fits' in x]
        results_list = [x for x in file_list if 'result.fits' in x]

        for i, (spectrum_fname, result_fname) in enumerate(zip(spectra_list, results_list)):

            if i > 5:
                break

            ax = fig.add_subplot(2,6,i+1+6*j)
            ax.tick_params(axis='both', labelsize=7)

            spectrum_tuple = load_a_spectrum(spectrum_fname)
            result_tuple = load_a_spectrum(result_fname)
            mol_name, Ju = parse_linestring(spectrum_fname.split('/')[-1].rstrip('_spectrum.fits'))

            ax.plot(spectrum_tuple[2], spectrum_tuple[0], 'k', lw=0.75, drawstyle='steps-mid')
            result_linestyle = 'r-'

            if fit_tuple is not None:
                # retrieve the relevant fit dict from fit_tuple
                this_fit = [x for x in fit_tuple[j].values() if x['Molecule'] == mol_name and x['Ju'] == Ju][0]
                if this_fit['n_lines'] == 2:

                    result_linestyle = 'C0--'

                    # plot the individual line fits
                    xs = np.arange(-30, 30, 0.1)
                    a = this_fit['t_peak'] 
                    b = this_fit['v_cen'] 
                    c = this_fit['v_fwhm'] / 2.35482 
                    ys = a * np.exp( - (xs-b)**2 / (2*c**2))

                    ax.plot(xs, ys, 'r', lw=0.75)

                    a2 = this_fit['t_peak_2'] 
                    b2 = this_fit['v_cen_2'] 
                    c2 = this_fit['v_fwhm_2'] / 2.35482 
                    ys2 = a2 * np.exp( - (xs-b2)**2 / (2*c2**2))

                    ax.plot(xs, ys2, 'C2', lw=0.75)                       

            ax.plot(result_tuple[2], result_tuple[0], result_linestyle, lw=0.75)

            ax.set_xlim(-12, 18)
            ax.set_ylim(ylims[j])

            plottable_string = plottable_latex_string(mol_name, Ju)
            ax.text(-10, text_heights[j], plottable_string, **text_params)

            if j!=1:
                ax.tick_params(axis='x', labelbottom='off')
            if i>0:
                ax.tick_params(axis='y', labelleft='off')
            if i==0 and j==1:
                ax.set_xlabel(r"$V_{\rm{lsr}}$ (km s$^{-1}$)", fontsize=8)
                ax.set_ylabel(r"$T_{\rm{mb}}$ (K)", fontsize=8)

    plt.tight_layout(w_pad=0.5, h_pad=0.3)
    plt.show()
    return fig


if __name__ == "__main__":

    if False:

        co_fit_tuple = co_baseline_spectra_and_compute_fits(verbose=False)

        paper_path = os.path.expanduser("~/Documents/Academia/Articles/Nitrogen_Paper/")
        fig_filename = "in_progress_graphics/hifi_co_lines.pdf"
        fig_fullpath = os.path.join(paper_path, fig_filename)

        fig = make_c18o_c17o_figure(co_fit_tuple)
        fig.savefig(fig_fullpath, bbox_inches='tight')