"""
Literally just loads the table so that I get the correct output.

Three S's. 
The Iras16293-2422 Millimeter And Submillimeter Spectral Survey.

"""

import os
import astropy.table
import astropy.units as u

path_to_timasss_table = os.path.expanduser("~/Documents/Code/reduce_herschel_IRAS16293_spectra/timasss_table")

timasss_table = astropy.table.QTable.read(os.path.join(path_to_timasss_table, "just_the_table.html"), data_start=2, header_start=0)

column_unit_dict = {
    'Frequency': u.MHz,
    'Eup': u.K,
    'Aij': u.s**-1,
    'Vo': u.km * u.s**-1,
    'FWHM': u.km * u.s**-1,
    'Int': u.K, 
    'Flux': u.K * u.km * u.s**-1}

for colname, unit in column_unit_dict.items():
    timasss_table[colname].unit = unit
    try:
        timasss_table["δ"+colname].unit = unit
    except:
        pass


if __name__ == "__main__":

    timasss_table.pprint()