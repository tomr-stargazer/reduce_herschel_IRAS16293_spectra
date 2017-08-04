""" 
Works with the output of the Herschel-fit lines to make something workable.

"""

import pdb

import numpy as np
import astropy.table
import astropy.units as u

from fit_the_lines_script import plottable_latex_string, which_hifi_band
from make_derived_line_properties_table import herschel_beamsize_from_freq
from groundbased_linefits import ground_beamsize_from_freq

def destructively_preprocess_dict_list(dict_list):
    """ 
    The latest version of astropy breaks my dict-table hacks. 
    So this is a hack around that. 

    """

    allowed_keys = dict_list[0].keys()
    disallowed_keys = []

    for dict_ in dict_list:
        for key in dict_.keys():
            if key not in allowed_keys:
                disallowed_keys.append(key)

    for key in disallowed_keys:
        for i in range(len(dict_list)):
            try:
                del dict_list[i][key]
            except KeyError:
                pass

    # it's in-place destructive, so:
    return None 


def make_linefit_table(fit_tuple):
    """ Takes the output of `baseline_spectra_and_compute_fits`, returns a table """

    hcn_fits, h13cn_fits, hc15n_fits = fit_tuple

    hcn_table = astropy.table.Table([x for x in hcn_fits.values()], names=hcn_fits[0].keys())
    # hcn_table['Molecule'] = np.array(['HCN']*len(hcn_table))
    hcn_table_re = hcn_table[('Molecule', 'Ju', *hcn_table.colnames[:-2])]

    h13cn_table = astropy.table.Table([x for x in h13cn_fits.values()], names=h13cn_fits[0].keys())
    # h13cn_table['Molecule'] = np.array(['H13CN']*len(h13cn_table))
    h13cn_table_re = h13cn_table[('Molecule', 'Ju', *h13cn_table.colnames[:-2])]

    hc15n_dict_list = [x for x in hc15n_fits.values()]
    destructively_preprocess_dict_list(hc15n_dict_list)
    hc15n_table = astropy.table.Table(hc15n_dict_list, names=hc15n_fits[0].keys())
    # hc15n_table['Molecule'] = np.array(['HC15N']*len(hc15n_table))
    hc15n_table_re = hc15n_table[('Molecule', 'Ju', *hc15n_table.colnames[:-2])]

    hcn_meta_table = astropy.table.vstack([hcn_table_re[:5], h13cn_table_re[:5], hc15n_table_re[:5]])
    hcn_meta_table['freq'].unit = u.GHz

    return hcn_meta_table


def turn_data_row_into_tex_row(table_row):

    linestring = plottable_latex_string(table_row['Molecule'], table_row['Ju'])

    band_name = "HIFI "+which_hifi_band(table_row['freq'])
    beamsize = herschel_beamsize_from_freq(table_row['freq']*u.GHz).to(u.arcsec).value

    # species & transition. frequency. telescope. theta main beam. V_0. FWHM. Peak flux. Area.
    tex_line = r"{0}{1} & {9:.4f} & {10} & {11:.1f} & {2:.2f} $\pm$ {3:.2f} & {4:.2f} $\pm$ {5:.2f} & {6:.3f} & {7:.2f} $\pm$ {8:.2f} \\".format(
        '', linestring, table_row['v_cen'], table_row['e_v_cen'], 
        table_row['v_fwhm'], table_row['e_v_fwhm'], table_row['t_peak'], table_row['area'], table_row['e_area'], table_row['freq'], band_name, beamsize)

    return tex_line


def prepare_linefit_table_for_latex(meta_table):

    list_of_latex_lines = [turn_data_row_into_tex_row(row) for row in meta_table]

    return list_of_latex_lines


def make_co_linefit_table(fit_tuple):
    """ Takes the output of `baseline_spectra_and_compute_fits`, returns a table """

    c18o_fits, c17o_fits = fit_tuple

    c18o_table = astropy.table.Table([x for x in c18o_fits.values()], names=c18o_fits[0].keys())
    c18o_table_re = c18o_table[('Molecule', 'Ju', *c18o_table.colnames[:-2])]

    c17o_dict_list = [x for x in c17o_fits.values()]
    destructively_preprocess_dict_list(c17o_dict_list)
    c17o_table = astropy.table.Table(c17o_dict_list, names=c17o_fits[0].keys())
    c17o_table_re = c17o_table[('Molecule', 'Ju', *c17o_table.colnames[:-2])]

    co_meta_table = astropy.table.vstack([c18o_table_re[:6], c17o_table_re[:6]])
    co_meta_table['freq'].unit = u.GHz

    return co_meta_table


def turn_enhanced_ground_row_into_tex_row(table_row):

    linestring = plottable_latex_string(table_row['Molecule'], table_row['Ju'])

    if table_row['Frequency'] < 300 * u.GHz:
        band_name = "IRAM-30m"
    elif table_row['Frequency'] > 300 * u.GHz:
        band_name = "JCMT 15m"
    beamsize = ground_beamsize_from_freq(table_row['Frequency']).to(u.arcsec).value

    # species & transition. frequency. telescope. theta main beam. V_0. FWHM. Peak flux. Area.
    tex_line = r"{0}{1} & {9:.4f} & {10} & {11:.1f} & {2:.2f} $\pm$ {3:.2f} & {4:.2f} $\pm$ {5:.2f} & {6:.3f} & {7:.2f} $\pm$ {8:.2f} \\".format(
        '', linestring, table_row['Vo'].value, table_row['δVo'].value, 
        table_row['FWHM'].value, table_row['δFWHM'].value, table_row['Int'].value, table_row['Flux'].value, table_row['δFlux'].value, table_row['Frequency'].to(u.GHz).value, band_name, beamsize)

    return tex_line
