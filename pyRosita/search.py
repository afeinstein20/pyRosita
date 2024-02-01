import os
import numpy as np
from tqdm import tqdm
from astropy import units
from astropy.table import Table
from astroquery.simbad import Simbad
import astropy.coordinates as coordinates

__all__ = ["SearchImages"]

class SearchImages(object):
    """
    Determines which SkyTile any target or list of targets fall on.

    Parameters
    ----------
    name : str, list, optional
       The name of the target or list of target names.
    coords : list, optional
       A set of coordinate(s) for a target or list of targets. Should be of form
       `[[ra, dec], [ra, dec]]`. Should be given in units of degrees.
    """

    def __init__(self, name=None, coords=None):
        self.name = name

        if coords is not None:
            coords = coordinates.SkyCoord(coords[:,0],
                                          coords[:,1], unit=units.deg)
            self.coords = coords

        # sets an individual name as a list
        if name is not None:
            if type(name) == str:
                self.name = [name]
            self.get_coordinates()

        self.fn_dir = os.path.join(os.path.expanduser('~'), '.pyRosita')
        self.catalog_fn = 'eRASS1_Main.v1.1.fits'
        self.skytile_fn = 'SKYMAPS_052022_MPE.fits'

        if os.path.isfile(os.path.join(self.fn_dir, self.catalog_fn)):
            self.open_source_catalog()
        else:
            self.download_source_catalog()
            #self.open_source_catalog()

        if os.path.isfile(os.path.join(self.fn_dir, self.skytile_fn)):
            self.open_skytile_catalog()
        else:
            self.download_skytile_catalog()
            self.open_skytile_catalog()


    def download_source_catalog(self):
        """
        Downloads the source catalog from eRosita DR1.

        Attributes
        ----------
        catalog_fn : str
           Filename of the Source catalog.
        """
        url = 'https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Main.tar.gz'

        if not os.path.exists(self.fn_dir):
            try:
                os.mkdir(self.fn_dir)
            except OSError:
                self.fn_dir = '.'
                warnings.warn('Warning: unable to create {}. '
                              'Downloading to the current '
                              'working directory instead.'.format(self.fn_dir))
        os.system('curl -L {0} -o {1}'.format(url,
                                              os.path.join(self.fn_dir, 'eRASS1_Main.tar.gz')))
        os.system('tar -xvzf {0}'.format(os.path.join(self.fn_dir,
                                                             'eRASS1_Main.tar.gz')))
        os.system('mv {0} {1}'.format(self.catalog_fn,
                                      os.path.join(self.fn_dir, self.catalog_fn)))


    def download_skytile_catalog(self):
        """
        Download the FITS file with the center (RA, Dec) and edges (RA, Dec) of
        each SkyTile available.

        Attributes
        ----------
        skytile_fn : str
           Filename of the FITS catalog.
        """
        url = "https://erosita.mpe.mpg.de/dr1/eSASS4DR1/eSASS4DR1_ProductsDescription/SKYMAPS_052022_MPE.fits"

        if not os.path.exists(self.fn_dir):
            try:
                os.mkdir(self.fn_dir)
            except OSError:
                self.fn_dir = '.'
                warnings.warn('Warning: unable to create {}. '
                              'Downloading to the current '
                              'working directory instead.'.format(self.fn_dir))
        os.system('curl -L {0} -o {1}'.format(url,
                                              os.path.join(self.fn_dir,
                                                           self.skytile_fn)))

    def open_skytile_catalog(self):
        """
        Opens the SkyTile catalog.

        Attributes
        ----------
        skytile_table : astropy.table.Table
           Astropy table of the center locations for each SkyTile.
        skytile_centers : astropy.coordinates.SkyCoord
           Coordinates object for the center (RA, Dec) for each SkyTile.
        """
        self.skytile_table = Table.read(os.path.join(self.fn_dir,
                                                     self.skytile_fn),
                                        format='fits')
        self.skytile_centers = coordinates.SkyCoord(self.skytile_table['RA_CEN'],
                                                    self.skytile_table['DE_CEN'],
                                                    unit=units.deg)
        return

    def get_coordinates(self):
        """
        Loops through and gets the coordinates for objects if names were
        passed into the function.
        """
        coords = np.zeros((len(self.name), 2), dtype='U30')

        for i in tqdm(range(len(self.name)), desc='Retrieving coordinates'):
            res = Simbad.query_object(self.name[i])
            coords[i][0] = res['RA'][0] + ' hours'
            coords[i][1] = res['DEC'][0] + ' degrees'

        self.coords = coordinates.SkyCoord(coordinates.Angle(coords[:,0]),
                                           coordinates.Angle(coords[:,1]))

    def find_closest_skytile(self, sep_val=1.9):
        """
        Crossmatches the coordinates of the targets in question with the centers
        of all available SkyTiles.

        Parameters
        ----------
        sep_val : float
           Largest separation from the center location of a SkyTile that a
           target would have observations.

        Attributes
        ----------
        tile_numbers : np.array
           Array of SkyTile numbers the targets would have observations on.
        separation : np.array
           Separation between the target and the center location of a SkyTile
           in units of degrees.
        match : np.array
           Array of 0s or 1s depending on if the target has a matching SkyTile.
           0 = not a match; 1 = a match.
        """
        tile_numbers = np.zeros(len(self.coords), dtype='U6')
        separation   = np.zeros(len(self.coords))
        match = np.zeros(len(self.coords), dtype=bool)

        for i in range(len(self.coords)):
            sep = self.coords[i].separation(self.skytile_centers).deg
            argmin = np.argmin(sep)

            tile_numbers[i] = self.skytile_table['SRVMAP'][argmin]
            separation[i]   = sep[argmin]

            if sep[argmin] <= sep_val:
                match[i] = True

        self.tile_numbers = tile_numbers
        self.separation   = separation
        self.match = match

        return

    def download_skytile(self, which_sources='all'):
        """
        Downloads the full SkyTile for a given target

        Parameters
        ----------
        which_sources : np.array, optional
           Sets which targets to download the data products for. Default is all
           targets under `self.coords` or `self.name`. If you do not data for
           all of the targets, pass in a `np.array` of integers for which indices
           the code should download.
        """
        return('This is a WIP.')
