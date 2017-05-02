"""
I've got some spectrum files that are output from CLASS
(via "fits write blah.fits /mode spectrum")
and I'd like to get them into python in a useful way.

"""

import numpy as np
from astropy.io.fits import getdata 
import astropy.io.fits

def load_a_spectrum(filename):

    data, header = getdata(filename, header=True)

    if type(data) == astropy.io.fits.fitsrec.FITS_rec:
        raise ValueError("This only accepts FITS spectra written in `/mode SPECTRUM`, not `/mode INDEX`.")

    central_freq = header['RESTFREQ']
    reference_pixel = header['CRPIX1']
    delta_freq = header['CDELT1']
    delta_v = header['DELTAV']

    spectrum = data[0][0][0]
    freq_array = (
        (np.arange(0, len(spectrum)) * delta_freq) + 
        (central_freq - (reference_pixel*delta_freq)))
    # downshift to km/s.
    vel_array = (0.001 * (
        (np.arange(0, len(spectrum)) * delta_v) + 
        (0 - reference_pixel*delta_v) ) )

    # for now we'll return these things, but eventually we'll maybe return an Object?
    # or at least the raw data/header thingies too?
    return spectrum, freq_array, vel_array