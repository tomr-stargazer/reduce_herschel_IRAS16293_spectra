""" 
Works with 

"""

import numpy as np
import astropy.table

def make_linefit_table(fit_tuple):
	""" Takes the output of `baseline_spectra_and_compute_fits`, returns a table """

	hcn_fits, h13cn_fits, hc15n_fits = fit_tuple

	hcn_table = astropy.table.Table([x for x in hcn_fits.values()], names=hcn_fits[0].keys())
	hcn_table['Molecule'] = np.array(['HCN']*len(hcn_table))
	hcn_table_re = hcn_table[('Molecule', 'line_name', *hcn_table.colnames[:-2])]

	h13cn_table = astropy.table.Table([x for x in h13cn_fits.values()], names=h13cn_fits[0].keys())
	h13cn_table['Molecule'] = np.array(['H13CN']*len(h13cn_table))
	h13cn_table_re = h13cn_table[('Molecule', 'line_name', *h13cn_table.colnames[:-2])]

	hc15n_table = astropy.table.Table([x for x in hc15n_fits.values()], names=hc15n_fits[0].keys())
	hc15n_table['Molecule'] = np.array(['HC15N']*len(hc15n_table))
	hc15n_table_re = hc15n_table[('Molecule', 'line_name', *hc15n_table.colnames[:-2])]

	hcn_meta_table = astropy.table.vstack([hcn_table_re[:5], h13cn_table_re[:5], hc15n_table_re[:5]])

	return hcn_meta_table


def turn_data_row_into_tex_row(table_row):

	# species & transition. frequency. telescope. theta main beam. V_0. FWHM. Peak flux. Area.
	tex_line = r"{0} {1} & freq & band & beamsize & {2:.2f} $\pm$ {3:.2f} & {4:.2f} $\pm$ {5:.2f} & {6:.3f} & {7:.2f} $\pm$ {8:.2f}".format(
		table_row['Molecule'], table_row['line_name'], table_row['v_cen'], table_row['e_v_cen'], 
		table_row['v_fwhm'], table_row['e_v_fwhm'], table_row['t_peak'], table_row['area'], table_row['e_area'])

	return tex_line


def prepare_linefit_table_for_latex(hcn_meta_table):

	list_of_latex_lines = [turn_data_row_into_tex_row(row) for row in hcn_meta_table]

	return list_of_latex_lines





