"""
This is really gonna pull everything together, I hope!

"""

import numpy as np
import astropy.units as u

from fit_the_lines_script import baseline_spectra_and_compute_fits
from linefit_table import make_linefit_table, make_co_linefit_table

from load_timasss_table import timasss_table
from groundbased_linefits import make_derived_props_table_ground, select_molecules, enhance_groundtable_with_mol_Ju_cols, make_co_derived_props_table_ground
from make_derived_line_properties_table import make_derived_props_table

def prepare_linefit_table_for_latex(table):

    list_of_latex_lines = [turn_data_row_into_tex_row(row) for row in table]

    return list_of_latex_lines


def turn_data_row_into_tex_row(table_row, ratio_12C_13C=69):

    tr = table_row

    # Ju= , tau_h13cn, N_u_h13cn, tau_HCN, N_u_hcn, | tau_hc15n, N_u_hc15n, h13cn / hc15n, 14N / 15N

    tex_line = (r"{0} & {1:.2f} & {2:.2e} & {3:.2f} & {4:.2e} & {5:.2f} & {6:.2e} & {7:.1f} & {8:.1f}".
                format(tr['J_upper'], tr['tau_h13cn'], tr['N(h13cn)_upper'], tr['tau_h13cn']*ratio_12C_13C, tr['N(h13cn)_upper']*ratio_12C_13C,
                       tr['tau_hc15n'], tr['N(hc15n)_upper'], tr['N(h13cn)_upper']/tr['N(hc15n)_upper'], tr['N(h13cn)_upper']/tr['N(hc15n)_upper']*ratio_12C_13C))

    # species & transition. frequency. telescope. theta main beam. V_0. FWHM. Peak flux. Area.
    # tex_line = r"{0}{1} & {9:.4f} & {10} & {11:.1f} & {2:.2f} $\pm$ {3:.2f} & {4:.2f} $\pm$ {5:.2f} & {6:.3f} & {7:.2f} $\pm$ {8:.2f} \\".format(
    #     '', linestring, table_row['Vo'].value, table_row['δVo'].value, 
    #     table_row['FWHM'].value, table_row['δFWHM'].value, table_row['Int'].value, table_row['Flux'].value, table_row['δFlux'].value, table_row['Frequency'].to(u.GHz).value, band_name, beamsize)

    return tex_line


def texify_derived_props_table(herschel_lineprops, ground_lineprops):

    [print(x) for x in prepare_linefit_table_for_latex(ground_lineprops)]
    [print(x) for x in prepare_linefit_table_for_latex(herschel_lineprops)]

    return


def get_N_HCN():

    herschel_fit_tuple = baseline_spectra_and_compute_fits(verbose=False)
    herschel_linefit_table = make_linefit_table(herschel_fit_tuple)
    herschel_lineprops = make_derived_props_table(herschel_linefit_table)

    ground_hcn_table = select_molecules(timasss_table, ['HCN', 'H13CN', 'HC15N'])
    enhance_groundtable_with_mol_Ju_cols(ground_hcn_table)
    ground_hcn_higher_lines= ground_hcn_table[ground_hcn_table['Ju']>2]
    ground_lineprops = make_derived_props_table_ground(ground_hcn_higher_lines)

    texify_derived_props_table(herschel_lineprops, ground_lineprops)

    # in principle we could spit out N_h13cn here too. It's more direct.
    # We will also want to condense things into a big table along the way, somehow.
    #  One for just the linefits, and then one for the derived properties.
    #  In principle it'd probably be more appropriate to first make that table,
    #  and *then* pull the N_HCN out of that, rather than duplicating the math
    #  in this function proper.

    f_c = 1.75
    carbon_12c_13c_ratio = 69
    N_HCN = f_c*carbon_12c_13c_ratio * (
                np.sum(herschel_lineprops['N(h13cn)_upper']) + 
                np.sum(ground_lineprops['N(h13cn)_upper']) )*u.cm**-2

    print("N_HCN = {:.3e}".format(N_HCN))

    return N_HCN

from c18o_fit_the_lines_script import co_baseline_spectra_and_compute_fits
from c18o_derived_line_properties_table import co_make_derived_props_table

# ok now for the other shoe, metaphorically
def get_N_C18O():

    herschel_fit_tuple = co_baseline_spectra_and_compute_fits(verbose=False)
    herschel_linefit_table = make_co_linefit_table(herschel_fit_tuple)
    herschel_lineprops = co_make_derived_props_table(herschel_linefit_table)

    ground_co_table = select_molecules(timasss_table, ['C18O', 'C17O'])
    enhance_groundtable_with_mol_Ju_cols(ground_co_table)
    ground_co_higher_lines= ground_co_table[ground_co_table['Ju']>1]
    ground_lineprops = make_co_derived_props_table_ground(ground_co_higher_lines)

    f_c = 1.45
    # carbon_12c_13c_ratio = 70
    N_C18O = f_c * (
                np.sum(herschel_lineprops['N(c18o)_upper']) + 
                np.sum(ground_lineprops['N(c18o)_upper']) )*u.cm**-2

    print("N_C18O = {:.3e}".format(N_C18O))

    return N_C18O


def get_N_H2(N_C18O):

    canonical_ratio = 500e4 # from Plume 2012
    N_H2 = canonical_ratio * N_C18O

    print("N_H2 = {:.3e}".format(N_H2))

    return N_H2



if __name__ == "__main__":

    N_HCN = get_N_HCN()

    N_H2 = get_N_H2(get_N_C18O())

    X_HCN = N_HCN/N_H2

    print("X_HCN = {:.3e}".format(X_HCN))

    HCN_over_Si = X_HCN * (1/2) * 3.16e4

    print("HCN/Si = {:.3e}".format(HCN_over_Si))
