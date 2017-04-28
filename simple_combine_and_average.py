"""
This is my proof-of-concept that I could, in principle, do the following even ONCE:

- find a two .hifi files from different polarizations
- make a simple CLASS script that combines & averages them
- run that CLASS script.

"""

from __future__ import division
import os
from subprocess import call

convert_fits_scriptname = "convert_fits_to_hifi.class"

def 