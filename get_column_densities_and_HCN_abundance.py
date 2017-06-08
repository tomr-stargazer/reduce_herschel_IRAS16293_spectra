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

def get_N_HCN():

    herschel_fit_tuple = baseline_spectra_and_compute_fits(verbose=False)
    herschel_linefit_table = make_linefit_table(herschel_fit_tuple)
    herschel_lineprops = make_derived_props_table(herschel_linefit_table)

    ground_hcn_table = select_molecules(timasss_table, ['HCN', 'H13CN', 'HC15N'])
    enhance_groundtable_with_mol_Ju_cols(ground_hcn_table)
    ground_hcn_higher_lines= ground_hcn_table[ground_hcn_table['Ju']>2]
    ground_lineprops = make_derived_props_table_ground(ground_hcn_higher_lines)

    # in principle we could spit out N_h13cn here too. It's more direct.
    # We will also want to condense things into a big table along the way, somehow.
    #  One for just the linefits, and then one for the derived properties.
    #  In principle it'd probably be more appropriate to first make that table,
    #  and *then* pull the N_HCN out of that, rather than duplicating the math
    #  in this function proper.

    f_c = 2
    carbon_12c_13c_ratio = 70
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

    f_c = 1.6
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