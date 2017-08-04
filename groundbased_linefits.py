"""
This takes the ground-based line fits from Caux et al. 2011 and puts them in a format 
amenable to combining with my Herschel-derived linefit data.

"""

import pdb
import numpy as np
import astropy.table
import astropy.units as u
import astropy.constants as c

from load_timasss_table import timasss_table
from make_derived_line_properties_table import iso_Nupper
from optical_depth_handling import tau_iso, tau_main, isotopic_ratio_from_taus


def select_molecule(table, molecule_name):

    # first several chars of the "Species & Transition" entry must equal the molname exactly
    padded_molname = molecule_name+" "
    relevant_rows = [name[:len(padded_molname)] == padded_molname for name in table['Species & Transition']]

    table_subset = table[relevant_rows]

    return table_subset


def select_molecules(table, list_of_molnames):

    # first several chars of the "Species & Transition" entry must equal the molname exactly
    relevant_rows = np.zeros(len(table), dtype=bool)
    for molecule_name in list_of_molnames:
        padded_molname = molecule_name+" "
        relevant_rows += np.array([name[:len(padded_molname)] == padded_molname for name in table['Species & Transition']])

    table_subset = table[relevant_rows]

    return table_subset


def select_line(table, line_name):

    relevant_rows = [line_name in name for name in table['Species & Transition']]
    table_subset = table[relevant_rows]

    return table_subset


def select_relevant_hcn_lines_ground(table):

    # hcn = select_molecule(table, 'HCN')
    # h13cn = select_molecule(table, 'H13CN')
    # hc15n = select_molecule(table, 'HC15N')
    # stacked_subset_of_tables = astropy.table.vstack([hcn, h13cn, hc15n])

    return select_molecules(timasss_table, ['HCN', 'H13CN', 'HC15N'])  


def which_dish(freq):
    JCMT_diameter = 15*u.m
    IRAM_30m_diameter = 30*u.m

    if u.Quantity(freq, u.MHz) < u.Quantity(300, u.GHz):
        return IRAM_30m_diameter
    else:
        return JCMT_diameter


def ground_beamsize_from_freq(freq):
    dish_diameter = which_dish(freq)
    wavelength = c.c / freq
    
    theta = 1.22 * (wavelength / dish_diameter) * u.rad
    return theta.to(u.arcsec)


def make_derived_props_table_ground(linefit_table):

    # so ultimately we're deriving N_upper.
    hcn_subtable = linefit_table[linefit_table['Molecule'] == 'HCN']
    h13cn_subtable = linefit_table[linefit_table['Molecule'] == 'H13CN']
    hc15n_subtable = linefit_table[linefit_table['Molecule'] == 'HC15N']

    Ju_column = hcn_subtable['Ju']
    h13cn_tau = tau_iso(h13cn_subtable['Int'], hcn_subtable['Int'])
    hc15n_tau = tau_iso(hc15n_subtable['Int'], hcn_subtable['Int'])

    theta_source = 1.29*u.arcsec # from ALMA image of h13cn 8-7 emission
    beamsizes = [ ground_beamsize_from_freq(row['Frequency']) for row in h13cn_subtable ]
    beamsizes_array = u.Quantity(beamsizes)
    eta_bf = theta_source**2 / (theta_source**2 + beamsizes_array**2)


    N_h13cn_upper_column = iso_Nupper(h13cn_subtable['Flux'], h13cn_tau, h13cn_subtable['Aij'], h13cn_subtable['Frequency'], eta_bf)
    N_hc15n_upper_column = iso_Nupper(hc15n_subtable['Flux'], hc15n_tau, hc15n_subtable['Aij'], hc15n_subtable['Frequency'], eta_bf)
    fractionation_column_ism = isotopic_ratio_from_taus(tau_main(h13cn_tau, 69), hc15n_tau)
    fractionation_column_solar = isotopic_ratio_from_taus(tau_main(h13cn_tau, 89), hc15n_tau)

    # return N_h13cn_upper_column
    new_table = astropy.table.Table([Ju_column, h13cn_tau, N_h13cn_upper_column, hc15n_tau, N_hc15n_upper_column, fractionation_column_ism, fractionation_column_solar], 
                                    names=['J_upper', 'tau_h13cn', "N(h13cn)_upper", "tau_hc15n", "N(hc15n)_upper", "14N/15N ratio (ISM)", "14N/15N ratio (solar)"])


    return new_table


def make_co_derived_props_table_ground(linefit_table):

    # so ultimately we're deriving N_upper.
    c18o_subtable = linefit_table[linefit_table['Molecule'] == 'C18O']
    c17o_subtable = linefit_table[linefit_table['Molecule'] == 'C17O']

    Ju_column = c18o_subtable['Ju']
    c18o_tau = 0.001 * np.ones_like(Ju_column)
    # h13cn_tau = tau_iso(h13cn_subtable['Int'], hcn_subtable['Int'])
    # hc15n_tau = tau_iso(hc15n_subtable['Int'], hcn_subtable['Int'])

    theta_source = 15.5*u.arcsec # from ALMA image of h13cn 8-7 emission
    beamsizes = [ ground_beamsize_from_freq(row['Frequency']) for row in c18o_subtable ]
    beamsizes_array = u.Quantity(beamsizes)
    eta_bf = theta_source**2 / (theta_source**2 + beamsizes_array**2)


    N_upper_column = iso_Nupper(c18o_subtable['Flux'], c18o_tau, c18o_subtable['Aij'], c18o_subtable['Frequency'], eta_bf)
    # fractionation_column_ism = isotopic_ratio_from_taus(tau_main(h13cn_tau, 69), hc15n_tau)
    # fractionation_column_solar = isotopic_ratio_from_taus(tau_main(h13cn_tau, 89), hc15n_tau)

    # return N_upper_column
    new_table = astropy.table.Table([Ju_column, N_upper_column], 
                                    names=['J_upper', "N(c18o)_upper"])

    return new_table

# let's try this again.

# WHAT WE WANT is to explicitly put a molecule and Ju columns
# so... like, given a 'Species & Transition' string, 
# under certain assumptions, rip out the species and Ju.

def molname_Ju_from_string(species_transition_string):

    mol_name, line_name = species_transition_string.split(' ')

    # Assumes transition strings that look like:
    # (4-3) 
    # or 
    # (11-01)
    Ju, Jl = line_name.strip("()").split("-")
    if int(Ju)-int(Jl) != 1:
        Ju = Ju[0]

    return mol_name, int(Ju)


# Don't put the WHOLE TIMASS table in here dude
def enhance_groundtable_with_mol_Ju_cols(table):

    mol_list = []
    Ju_list = []

    for name in table['Species & Transition']:
        mol, Ju = molname_Ju_from_string(name)
        mol_list.append(mol)
        Ju_list.append(Ju)

    table['Molecule'] = mol_list
    table['Ju'] = Ju_list

    return 

