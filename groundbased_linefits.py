"""
This takes the ground-based line fits from Caux et al. 2011 and puts them in a format 
amenable to combining with my Herschel-derived linefit data.

"""

import numpy as np
import astropy.table
import astropy.units as u
import astropy.constants as c

from load_timasss_table import timasss_table
from make_derived_line_properties_table import h13cn_Nupper
from optical_depth_handling import tau_iso, tau_main, isotopic_ratio_from_taus


def select_molecule(table, molecule_name):

    # first several chars of the "Species & Transition" entry must equal the molname exactly
    padded_molname = molecule_name+" "
    relevant_rows = [name[:len(padded_molname)] == padded_molname for name in table['Species & Transition']]

    table_subset = table[relevant_rows]

    return table_subset


def select_line(table, line_name):

    relevant_rows = [line_name in name for name in table['Species & Transition']]
    table_subset = table[relevant_rows]

    return table_subset


# actually - this is wrong. 
# I'm not doing it by row, I'm doing it by species.
def Nupper_from_row(table, transition_name, source_size):
    pass


def select_relevant_hcn_lines_ground(table):

    hcn = select_molecule(table, 'HCN')
    h13cn = select_molecule(table, 'H13CN')
    hc15n = select_molecule(table, 'HC15N')

    stacked_subset_of_tables = astropy.table.vstack([hcn, h13cn, hc15n])

    return stacked_subset_of_tables


def which_dish(freq):
    JCMT_diameter = 15*u.m
    IRAM_30m_diameter = 30*u.m

    if freq < 300 * u.GHz:
        return IRAM_30m_diameter
    else:
        return JCMT_diameter


def ground_beamsize_from_freq(freq):
    dish_diameter = which_dish(freq)
    wavelength = c.c / freq
    
    theta = 1.22 * (wavelength / dish_diameter) * u.rad
    return theta.to(u.arcsec)


def make_derived_props_table_ground(stacked_table):

    # so ultimately we're deriving N_upper.
    lines = ['(3-2)', '(4-3)']
    hcn_table_34 = astropy.table.vstack([ select_line(select_molecule(stacked_table, 'HCN'), line) for line in lines ])
    h13cn_table_34 = astropy.table.vstack([ select_line(select_molecule(stacked_table, 'H13CN'), line) for line in lines ])
    hc15n_table_34 = astropy.table.vstack([ select_line(select_molecule(stacked_table, 'HC15N'), line) for line in lines ])
    
    Ju_column = astropy.table.Column(np.array([3,4]), name='Ju')
    hcn_table_34.add_column(Ju_column)
    h13cn_table_34.add_column(Ju_column)
    hc15n_table_34.add_column(Ju_column)

    theta_source = 1.29*u.arcsec # from ALMA image of h13cn 8-7 emission
    eta_bf = theta_source**2 / (theta_source**2 + ground_beamsize_from_freq(h13cn_table_34['Frequency']*u.GHz)**2)


    pass



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

    return mol_name, Ju


