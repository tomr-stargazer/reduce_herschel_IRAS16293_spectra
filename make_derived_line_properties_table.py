"""
We wanna make a table with these derived properties:

transition (Ju-Jl)
tau_h13cn
N_upper, hcn (maybe multiple values for different assumed 12c/13c ratios)
14N/15N ratio

This is a place where the assumed source size will have to come into play!

"""

import numpy as np
import astropy.table
import astropy.units as u
import astropy.constants as c

from optical_depth_handling import tau_iso, tau_main, isotopic_ratio_from_taus

def iso_Nupper(flux, tau_iso, A_ul, freq, eta_bf):
    """ 
    Calculates the N_upper for a chosen isotopologue i.e. probably not the main species

    """

    flux = u.Quantity(flux, u.K * u.km/u.s)
    freq = u.Quantity(freq, u.GHz)
    A_ul = u.Quantity(A_ul, u.s**-1)

    first_term = 8*np.pi*c.k_B* freq**2 / (A_ul * c.h * c.c**3)
    second_term = flux / eta_bf
    third_term = tau_iso / (1 - np.exp(-tau_iso))

    Nupper = first_term * second_term * third_term

    return Nupper.to(u.cm**-2)

def herschel_beamsize_from_freq(freq):
    herschel_diameter = 3.5 * u.m
    wavelength = c.c / freq
    
    theta = 1.22 * (wavelength / herschel_diameter) * u.rad
    return theta.to(u.arcsec)


A_10 = 2.2256e-05
A_21 = 2.1360e-04
A_32 = 7.7240e-04
A_43 = 1.8984e-03
A_54 = 3.7923e-03
A_65 = 6.6521e-03
A_76 = 1.0680e-02
A_87 = 1.6074e-02
A_98 = 2.3033e-02
A_ten9 = 3.1755e-02

# h13cn_Auls = [A_10, A_32, A_43, A_65, A_76, A_87, A_98, A_ten9]
h13cn_Auls = [A_65, A_76, A_87, A_98, A_ten9]
h13cn_Auls = [x / u.s for x in h13cn_Auls]

# !TRANS + UP + LOW + EINSTEINA(s^-1) + FREQ(GHz) + E_u(K)
#     1     2     1   2.2029e-05      86.0549610      4.13
#     2     3     2   2.1146e-04     172.1079560     12.39
#     3     4     3   7.6473e-04     258.1571000     24.78
#     4     5     4   1.8793e-03     344.2003199     41.30
#     5     6     5   3.7541e-03     430.2357675     61.95
#     6     7     6   6.5851e-03     516.2614602     86.72
#     7     8     7   1.0572e-02     602.2754469    115.63
#     8     9     8   1.5913e-02     688.2757768    148.66
#     9    10     9   2.2799e-02     774.2604988    185.82
#    10    11    10   3.1433e-02     860.2276619    227.10
#    11    12    11   4.2016e-02     946.1753151    272.51

hc15n_A_10 = 2.2029e-05
hc15n_A_21 = 2.1146e-04
hc15n_A_32 = 7.6473e-04
hc15n_A_43 = 1.8793e-03
hc15n_A_54 = 3.7541e-03
hc15n_A_65 = 6.5851e-03
hc15n_A_76 = 1.0572e-02
hc15n_A_87 = 1.5913e-02
hc15n_A_98 = 2.2799e-02
hc15n_A_ten9 = 3.1433e-02

# hc15n_Auls = [A_10, A_32, A_43, A_65, A_76, A_87, A_98, A_ten9]
hc15n_Auls = [hc15n_A_65, hc15n_A_76, hc15n_A_87, hc15n_A_98, hc15n_A_ten9]
hc15n_Auls = [x / u.s for x in hc15n_Auls]

# currently only handles the Herschel data.
def make_derived_props_table(linefit_table):

    hcn_subtable = linefit_table[linefit_table['Molecule'] == 'HCN']
    h13cn_subtable = linefit_table[linefit_table['Molecule'] == 'H13CN']
    hc15n_subtable = linefit_table[linefit_table['Molecule'] == 'HC15N']

    Ju_column = hcn_subtable['Ju']
    h13cn_tau = tau_iso(h13cn_subtable['t_peak'], hcn_subtable['t_peak'])
    hc15n_tau = tau_iso(hc15n_subtable['t_peak'], hcn_subtable['t_peak'])

    theta_source = 1.29*u.arcsec # from ALMA image of h13cn 8-7 emission
    eta_bf = theta_source**2 / (theta_source**2 + herschel_beamsize_from_freq(h13cn_subtable['freq'])**2)

    N_h13cn_upper_column = iso_Nupper(h13cn_subtable['area'], h13cn_tau, h13cn_Auls, h13cn_subtable['freq'], eta_bf)
    N_hc15n_upper_column = iso_Nupper(hc15n_subtable['area'], hc15n_tau, hc15n_Auls, hc15n_subtable['freq'], eta_bf)
    fractionation_column_ism = isotopic_ratio_from_taus(tau_main(h13cn_tau, 69), hc15n_tau)
    fractionation_column_solar = isotopic_ratio_from_taus(tau_main(h13cn_tau, 89), hc15n_tau)
    fractionation_column_ism[-1] = np.nan
    fractionation_column_solar[-1] = np.nan

    # return N_h13cn_upper_column
    new_table = astropy.table.Table([Ju_column, h13cn_tau, N_h13cn_upper_column, hc15n_tau, N_hc15n_upper_column, fractionation_column_ism, fractionation_column_solar], 
                                    names=['J_upper', 'tau_h13cn', "N(h13cn)_upper", "tau_hc15n", "N(hc15n)_upper", "14N/15N ratio (ISM)", "14N/15N ratio (solar)"])

    return new_table

