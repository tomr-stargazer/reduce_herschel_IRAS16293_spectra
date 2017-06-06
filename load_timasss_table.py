"""
Literally just loads the table so that I get the correct outcome.

Three S's. 
The Iras16293-2422 Millimeter And Submillimeter Spectral Survey.

"""

import os
import astropy.table

path_to_timasss_table = os.path.expanduser("~/Documents/Code/reduce_herschel_IRAS16293_spectra/timass_table")

timasss_table = astropy.table.Table.read(os.path.join(path_to_timasss_table, "just_the_table.html"), data_start=2, header_start=0)

if __name__ == "__main__":

    timasss_table.pprint()