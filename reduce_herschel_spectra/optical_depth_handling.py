"""
Some functions to deal with getting optical depths out of this.

"""

import numpy as np

def tau_iso(T_iso, T_main):
    """ See Eq. 1 from Crockett+14a. """
    return -np.log(1 - T_iso/T_main)


def tau_main(tau_iso, isotopic_ratio):
    """ See Eq. 2 from Crockett+14a. """
    return tau_iso * isotopic_ratio


def isotopic_ratio_from_taus(tau_main, tau_iso):

    ratio = tau_main / tau_iso
    return ratio


# let's do a rough approach. We wanna compile the following things:

# for the following situations
# - every line
# - each 12c/13c ratio under consideration

# given an input table of line fits, which might possibly include 

# ok so now I'm expanding this a bit. Our underlying table really should have the following unambiguous columns:
# J_u (int)
# Molecule name (str, in a clear format -- probably just like H13CN)
# so that we can iterate over both "J_u" and "Molecule name" and match things appropriately.
# This may require editing our linelists.

def compute_taus_and_nitrogen_fractions(linefit_table):

    # takes in an astropy table

    # returns two or three things I guess:
    # a. h13cn optical depths for each line
    # b. 
    pass

def compute_taus(linefit_table):
    """
    Given a table of hcn and h13cn lines, computes the h13cn tau for each line.

    Returns it as a column/array thing.
    """

    hcn_subset = linefit_table[linefit_table['Molecule'] == 'HCN']
    h13cn_subset = linefit_table[linefit_table['Molecule'] == 'H13CN']

    h13cn_tau = tau_iso(h13cn_subset['t_peak'], hcn_subset['t_peak'])

    # we return the data to strip out the erroneous column metadata. It's unitless.
    return h13cn_tau.data


# def make_