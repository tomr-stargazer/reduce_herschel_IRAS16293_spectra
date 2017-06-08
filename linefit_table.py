""" 
Works with the output of the Herschel-fit lines to make something workable.

"""

import pdb

import numpy as np
import astropy.table
import astropy.units as u

from fit_the_lines_script import plottable_latex_string

def make_linefit_table(fit_tuple):
    """ Takes the output of `baseline_spectra_and_compute_fits`, returns a table """

    hcn_fits, h13cn_fits, hc15n_fits = fit_tuple

    hcn_table = astropy.table.Table([x for x in hcn_fits.values()], names=hcn_fits[0].keys())
    # hcn_table['Molecule'] = np.array(['HCN']*len(hcn_table))
    hcn_table_re = hcn_table[('Molecule', 'Ju', *hcn_table.colnames[:-2])]

    h13cn_table = astropy.table.Table([x for x in h13cn_fits.values()], names=h13cn_fits[0].keys())
    # h13cn_table['Molecule'] = np.array(['H13CN']*len(h13cn_table))
    h13cn_table_re = h13cn_table[('Molecule', 'Ju', *h13cn_table.colnames[:-2])]

    hc15n_table = astropy.table.Table([x for x in hc15n_fits.values()], names=hc15n_fits[0].keys())
    # hc15n_table['Molecule'] = np.array(['HC15N']*len(hc15n_table))
    hc15n_table_re = hc15n_table[('Molecule', 'Ju', *hc15n_table.colnames[:-2])]

    hcn_meta_table = astropy.table.vstack([hcn_table_re[:5], h13cn_table_re[:5], hc15n_table_re[:5]])
    hcn_meta_table['freq'].unit = u.GHz

    return hcn_meta_table


def turn_data_row_into_tex_row(table_row):

    linestring = plottable_latex_string(table_row['Molecule'], table_row['Ju'])

    # species & transition. frequency. telescope. theta main beam. V_0. FWHM. Peak flux. Area.
    tex_line = r"{0}{1} & freq & band & beamsize & {2:.2f} $\pm$ {3:.2f} & {4:.2f} $\pm$ {5:.2f} & {6:.3f} & {7:.2f} $\pm$ {8:.2f} \\".format(
        '', linestring, table_row['v_cen'], table_row['e_v_cen'], 
        table_row['v_fwhm'], table_row['e_v_fwhm'], table_row['t_peak'], table_row['area'], table_row['e_area'])

    return tex_line


def prepare_linefit_table_for_latex(hcn_meta_table):

    list_of_latex_lines = [turn_data_row_into_tex_row(row) for row in hcn_meta_table]

    return list_of_latex_lines



def make_co_linefit_table(fit_tuple):
    """ Takes the output of `baseline_spectra_and_compute_fits`, returns a table """

    c18o_fits, c17o_fits = fit_tuple

    c18o_table = astropy.table.Table([x for x in c18o_fits.values()], names=c18o_fits[0].keys())
    c18o_table_re = c18o_table[('Molecule', 'Ju', *c18o_table.colnames[:-2])]

    c17o_table = astropy.table.Table([x for x in c17o_fits.values()], names=c17o_fits[0].keys())
    c17o_table_re = c17o_table[('Molecule', 'Ju', *c17o_table.colnames[:-2])]

    co_meta_table = astropy.table.vstack([c18o_table_re[:6], c17o_table_re[:6]])
    co_meta_table['freq'].unit = u.GHz

    return co_meta_table

