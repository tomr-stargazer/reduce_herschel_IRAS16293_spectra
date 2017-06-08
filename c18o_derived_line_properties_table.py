"""
Now it's c18o.

We wanna make a table with these derived properties:

transition (Ju-Jl)
tau_c18o
N_upper, hcn (maybe multiple values for different assumed 12c/13c ratios)
14N/15N ratio

This is a place where the assumed source size will have to come into play!

"""

import pdb

import numpy as np
import astropy.table
import astropy.units as u
import astropy.constants as c

from optical_depth_handling import tau_iso, tau_main, isotopic_ratio_from_taus

def c18o_Nupper(flux, tau_c18o, A_ul, freq, eta_bf):

    flux = u.Quantity(flux, u.K * u.km/u.s)
    freq = u.Quantity(freq, u.GHz)
    A_ul = u.Quantity(A_ul, u.s**-1)

    first_term = 8*np.pi*c.k_B* freq**2 / (A_ul * c.h * c.c**3)
    second_term = flux / eta_bf
    third_term = tau_c18o / (1 - np.exp(-tau_c18o))

    Nupper = first_term * second_term * third_term

    return Nupper.to(u.cm**-2)


def herschel_beamsize_from_freq(freq):
        herschel_diameter = 3.5 * u.m
        wavelength = c.c / freq
        
        theta = 1.22 * (wavelength / herschel_diameter) * u.rad
        return theta.to(u.arcsec)

"""
from http://home.strw.leidenuniv.nl/~moldata/datafiles/c18o.dat
!TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)
    1     2     1   6.266e-08      109.7821734       5.27
    2     3     2   6.011e-07      219.5603541      15.81
    3     4     3   2.172e-06      329.3305525      31.61
    4     5     4   5.330e-06      439.0887658      52.68
    5     6     5   1.062e-05      548.8310055      79.02
    6     7     6   1.860e-05      658.5532782     110.63
    7     8     7   2.978e-05      768.2515933     147.50
    8     9     8   4.468e-05      877.9219553     189.63
    9    10     9   6.380e-05      987.5603822     237.03
   10    11    10   8.762e-05     1097.1628753     289.68
   11    12    11   1.166e-04     1206.7254487     347.60
"""


A_10 = 6.266e-08
A_21 = 6.011e-07
A_32 = 2.172e-06
A_43 = 5.330e-06
A_54 = 1.062e-05
A_65 = 1.860e-05
A_76 = 2.978e-05
A_87 = 4.468e-05
A_98 = 6.380e-05
A_ten9 = 8.762e-05

# c18o_Auls = [A_10, A_32, A_43, A_65, A_76, A_87, A_98, A_ten9]
c18o_Auls = [A_54, A_65, A_76, A_87, A_98, A_ten9]
c18o_Auls = [x / u.s for x in c18o_Auls]

# currently only handles the Herschel data.
def co_make_derived_props_table(linefit_table):

    c18o_subtable = linefit_table[linefit_table['Molecule'] == 'C18O']
    c17o_subtable = linefit_table[linefit_table['Molecule'] == 'C17O']

    Ju_column = c18o_subtable['Ju']
    # skipping taus. assuming it's all optically thin.
    c18o_tau = 0.0001 * np.ones_like(Ju_column)
    # c18o_tau = tau_iso(c18o_subtable['t_peak'], hcn_subtable['t_peak'])

    theta_source = 15.5*u.arcsec # from SMA image of c18o 2-1 emission, Jørgensen+2011
    eta_bf = theta_source**2 / (theta_source**2 + herschel_beamsize_from_freq(c18o_subtable['freq'])**2)

    N_upper_column = c18o_Nupper(c18o_subtable['area'], c18o_tau, c18o_Auls, c18o_subtable['freq'], eta_bf)
    # not doing fractionation...
    # fractionation_column_ism = isotopic_ratio_from_taus(tau_main(c18o_tau, 69), hc15n_tau)
    # fractionation_column_solar = isotopic_ratio_from_taus(tau_main(c18o_tau, 89), hc15n_tau)
    # fractionation_column_ism[-1] = np.nan
    # fractionation_column_solar[-1] = np.nan

    # return N_upper_column
    new_table = astropy.table.Table([Ju_column, N_upper_column], 
                                    names=['J_upper', "N(c18o)_upper"])

    pdb.set_trace()

    return new_table

