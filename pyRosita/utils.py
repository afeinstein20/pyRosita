# Series of useful eROSITA helper functions
import numpy as np

__all__ = ['telescope_ID', 'filter_wheel_num', 'energy_band']


# These functions disentangle information about the observing strategy from
#    a file's name. More information on how the naming conventions were created
#    can be found here: https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_ProductsDescription/file_naming_scheme_dr1.html

def telescope_ID(k):
    """Uses information from a given file name to determine what Telescope ID
    was used for that observation."""
    if k == 0:
        return('All cameras/telescopes')
    elif k == 8:
        return('Cameras with on-chip filter')
    elif k == 9:
        return('Cameras without on-chip filter')
    else:
        return('Individual cameras/telescopes')

def filter_wheel_num(l):
    """Uses information from a given file name to determine what filter wheel
    was used for that observation."""
    if l == 0:
        return('Closed')
    elif l == 1:
        return('Open')
    elif l == 2:
        return('Filter')
    elif l == 3:
        return('Calibration source')
    else:
        return('Undefined')

def energy_band(m):
    """Uses information from a given file name to determine what energy band the
    observations were taken in. Energy bands are given in units of keV."""
    if m == 0:
        return([0.2, 10.0])
    elif m == 1:
        return([0.1, 0.6])
    elif m == 2:
        return([0.6, 2.3])
    elif m == 3:
        return([2.3, 5.0])
    elif m == 4:
        return([0.2, 2.3])
    elif m == 5:
        return([0.2, 0.5])
    elif m == 6:
        return([0.5, 1.0])
    elif m == 7:
        return([1.0, 2.0])
    else:
        return('Invalid energy band key number.')
