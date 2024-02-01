import os
import tarfile
import numpy as np
from tqdm import tqdm
from astropy import units
from astropy.table import Table
from astroquery.simbad import Simbad
import astropy.coordinates as coordinates

__all__ = ["SearchDatabase"]

class SearchDatabase(object):
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
        self.product_dir = os.path.join(self.fn_dir, 'data_products')
        self.tile_dir = os.path.join(self.fn_dir, 'skytiles')

        self.catalog_fn = 'eRASS1_Main.v1.1.fits'
        self.skytile_fn = 'SKYMAPS_052022_MPE.fits'

        if os.path.isfile(os.path.join(self.fn_dir, self.catalog_fn)):
            self.open_source_catalog()
        else:
            self.download_source_catalog()
            self.open_source_catalog()

        if os.path.isfile(os.path.join(self.fn_dir, self.skytile_fn)):
            self.open_skytile_catalog()
        else:
            self.download_skytile_catalog()
            self.open_skytile_catalog()

    def create_directories(self, path):
        """
        Creates a network of directories to store the data to.

        Parameters
        ----------
        path : str
           The path of the directory to create.
        """
        if not os.path.exists(path):
            try:
                os.mkdir(path)
            except OSError:
                path = '.'
                warnings.warn('Warning: unable to create {}. '
                              'Downloading to the current '
                              'working directory instead.'.format(path))
        return

    def download_source_catalog(self):
        """
        Downloads the source catalog from eRosita DR1.

        Attributes
        ----------
        catalog_fn : str
           Filename of the Source catalog.
        """
        url = 'https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/MerloniA_DR1/eRASS1_Main.tar.gz'

        self.create_directories(self.fn_dir) # creates directory if needed

        os.system('curl -L {0} -o {1}'.format(url,
                                              os.path.join(self.fn_dir, 'eRASS1_Main.tar.gz')))
        os.system('tar -xvzf {0}'.format(os.path.join(self.fn_dir,
                                                             'eRASS1_Main.tar.gz')))
        os.system('mv {0} {1}'.format(self.catalog_fn,
                                      os.path.join(self.fn_dir, self.catalog_fn)))
        return

    def open_source_catalog(self):
        """
        Opens the source catalog.

        Attributes
        ----------
        source_cat : astropy.table.Table
           Source catalog.
        source_cat_coords : astropy.coordinates.SkyCoord
           (RA, Dec) coordinates for all targets in the source catalog.
        """
        self.source_cat = Table.read(os.path.join(self.fn_dir, self.catalog_fn),
                                     format='fits')
        self.source_cat_coords = coordinates.SkyCoord(self.source_cat['RA'],
                                                      self.source_cat['DEC'],
                                                      unit=units.deg)
        return

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

        self.create_directories(self.fn_dir) # creates directory if needed

        os.system('curl -L {0} -o {1}'.format(url,
                                              os.path.join(self.fn_dir,
                                                           self.skytile_fn)))
        return

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
        return

    def find_closest_source(self, sep_val=1.0*units.arcmin):
        """
        Crossmatches the coordinates of the targets in question with the source
        catalog.

        Parameters
        ----------
        sep_val : float, optional
           Largest separation from the target that could be considered a match.
           Default is 1 arcminute.
        closest_sources : astropy.table.Table
           Table of the closest matched sources. Sources which do not have a match
           within the value of `sep_val` are assigned a row of 0s, to preserve
           order.
        """
        closest_sources = Table(names=self.source_cat.colnames,
                                dtype=self.source_cat.dtype)

        seps = np.zeros(len(self.coords))

        for i in tqdm(range(len(self.coords))):
            sep = self.coords[i].separation(self.source_cat_coords).deg
            argmin = np.argmin(sep)
            if (sep[argmin]*units.deg).to(sep_val.unit) <= sep_val:
                closest_sources.add_row(self.source_cat[argmin])
                seps[i] = sep[argmin]
            else:
                blank = np.ma.array(np.zeros(len(self.source_cat.colnames)),
                                    mask = np.ones(len(self.source_cat.colnames)))
                closest_sources.add_row(blank)

        self.closest_sources = closest_sources
        return

    def download_data_products(self, idx='all'):
        """
        Downloads the data products for a given eROSITA source ID.

        Parameters
        ----------
        idx : np.array, optional
           An array of indices corresponding to which sources the user wants
           downloaded. Default is 'all' (downloads data for all sources).
        """
        url = "https://erosita.mpe.mpg.de/dr1/erodat/data/download_source/{0}/"

        self.create_directories(self.product_dir) # creates directory if needed

        if idx == 'all':
            idx = np.arange(0, len(self.coords), 1, dtype=int)

        for i in idx:
            detuid = self.closest_sources['DETUID'][i]

            if detuid != '0.0':
                download_url = url.format(detuid)
                tar_path = os.path.join(self.product_dir, detuid+'.tar.xz')
                file_path = os.path.join(self.product_dir, detuid)

                # if the directory of files does not exist
                if not os.path.exists(file_path):

                    # if the tar file does not exist, download it
                    if not os.path.exists(tar_path):
                        os.system('curl -L {0} -o {1}'.format(download_url,
                                                              tar_path))

                    # open the tar file
                    with tarfile.open(tar_path, 'r:xz') as tar:
                        tar.extractall(path=file_path)

                    os.system('rm {0}'.format(tar_path))

                else:
                    print('Files for {} are already downloaded.'.format(self.closest_sources['IAUNAME'][i]))
        return


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
        tile_sep : np.array
           Separation between the target and the center location of a SkyTile
           in units of degrees.
        tile_match : np.array
           Array of 0s or 1s depending on if the target has a matching SkyTile.
           0 = not a match; 1 = a match.
        """
        tile_numbers = np.zeros(len(self.coords), dtype='U6')
        tile_sep     = np.zeros(len(self.coords))
        tile_match   = np.zeros(len(self.coords), dtype=bool)

        for i in range(len(self.coords)):
            sep = self.coords[i].separation(self.skytile_centers).deg
            argmin = np.argmin(sep)

            tile_numbers[i] = self.skytile_table['SRVMAP'][argmin]
            tile_sep[i]     = sep[argmin]

            if sep[argmin] <= sep_val:
                tile_match[i] = True

        self.tile_numbers = tile_numbers
        self.tile_sep     = tile_sep
        self.tile_match   = tile_match

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
        self.create_directories(self.tile_dir) # creates directory if needed

        return('This is a WIP.')
