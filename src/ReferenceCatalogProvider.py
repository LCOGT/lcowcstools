import abc
import numpy as np
import sqlite3
import math
from astropy.table import Table
import logging
import os.path
import astropy.units as u
from astropy.coordinates import SkyCoord
from astroquery.gaia import Gaia

__author__ = 'drharbeck@gmail.com'

log = logging.getLogger(__name__)


class ReferenceCatalogProvider(metaclass=abc.ABCMeta):
    """ Base class to provide a Reference on sky catalog with RA/Dec coordinates

        The purpose is to provide an abstract interface to read a reference catalog.
    """

    @abc.abstractmethod
    def get_reference_catalog(self, ra, dec, radius) -> Table:
        """

        :param ra:
        :param dec:
        :param radius:
        :return: astropy Table object that has at least the columns 'RA' and 'Dec'
        """
        pass


class gaiaonline(ReferenceCatalogProvider):
    ''' Read the GAIA catalog via astroquery

        Requires a network connection for the online query.
    '''

    def get_reference_catalog(self, ra, dec, radius):
        coord = SkyCoord(ra=ra, dec=dec, unit=(u.degree, u.degree), frame='icrs')
        radius = u.Quantity(radius, u.deg)
        j = Gaia.cone_search_async(coord, radius)
        r = j.get_results()
        # astroquery already assigns units to the resulting columns, but we do not want that.
        retTable = Table([r['ra'] / u.degree, r['dec'] / u.degree], names=["RA", 'Dec'])
        log.debug("Gaia reference catalog has %d items" % len(retTable))
        return retTable


class refcat2(ReferenceCatalogProvider):
    """Interface to query consolidated  sqlite3 db file that was generated from Tonry (2018) refcat2.

        Requires access to the offline database somewhere in the file system
    """

    FILTERMAPPING = {}
    FILTERMAPPING['gp'] = {'refMag': 'g', 'colorTerm': 0.0, 'airmassTerm': 0.20, 'defaultZP': 0.0}
    FILTERMAPPING['rp'] = {'refMag': 'r', 'colorTerm': 0.0, 'airmassTerm': 0.12, 'defaultZP': 0.0}
    FILTERMAPPING['ip'] = {'refMag': 'i', 'colorTerm': 0.0, 'airmassTerm': 0.08, 'defaultZP': 0.0}
    FILTERMAPPING['zp'] = {'refMag': 'z', 'colorTerm': 0.0, 'airmassTerm': 0.05, 'defaultZP': 0.0}

    ###  PS to SDSS color transformations according to  Finkbeiner 2016
    ###  http://iopscience.iop.org/article/10.3847/0004-637X/822/2/66/meta#apj522061s2-4 Table 2
    ###  Note that this transformation is valid for stars only. For the purpose of photometric
    ###  calibration, it is desirable to select point sources only from the input catalog.

    ## Why reverse the order of the color term entries? Data are entered in the order as they are
    ## shown in paper. Reverse after the fact to avoid confusion when looking at paper

    ps1colorterms = {}
    ps1colorterms['g'] = [-0.01808, -0.13595, +0.01941, -0.00183][::-1]
    ps1colorterms['r'] = [-0.01836, -0.03577, +0.02612, -0.00558][::-1]
    ps1colorterms['i'] = [+0.01170, -0.00400, +0.00066, -0.00058][::-1]
    ps1colorterms['z'] = [-0.01062, +0.07529, -0.03592, +0.00890][::-1]

    def __init__(self, dbfile):
        """

        :param dbfile: location of the consolidated refcat2 slite database
        """

        if (dbfile is None) or (not os.path.isfile(dbfile)):
            log.error("Unable to find reference catalog: %s" % (str(dbfile)))
            self.dbfile = None
            return
        self.dbfile = dbfile

    def isInCatalogFootprint(self, ra, dec):
        return True


    def PStoSDSS(self, table):
        """
        Modify table in situ from PS1 to SDSS, requires column names compatible with ps1colorterms definition.

        :param table:
        :return: modified table.
        """
        if table is not None:
            pscolor = table['g'] - table['i']
            for filter in self.ps1colorterms:
                colorcorrection = np.polyval(self.ps1colorterms[filter], pscolor)
                table[filter] -= colorcorrection

        return table

    def get_reference_catalog(self, ra, dec, radius):

        rows = None
        try:
            connection = sqlite3.connect(self.dbfile)
            cursor = connection.cursor()
            if (not radius == None) and (radius > 0):
                min_dec = dec - radius
                max_dec = dec + radius
                min_ra = ra - radius / math.cos(math.radians(dec))
                max_ra = ra + radius / math.cos(math.radians(dec))

            sql_command = 'select sources.RA, sources.Dec, sources.g,sources.r,sources.i, sources.z from sources, positions ' \
                          'where positions.ramin >= {ramin} and positions.ramax <= {ramax} ' \
                          'and positions.decmin >= {decmin} and positions.decmax <= {decmax} ' \
                          'and positions.objid = sources.objid'

            sql_command = sql_command.format(ramin=min_ra, ramax=max_ra, decmin=min_dec, decmax=max_dec)
            cursor.execute(sql_command)
            rows = np.asarray(cursor.fetchall())
            table = Table(rows, names=['RA', 'Dec', 'g', 'r', 'i', 'z'])
            cursor.close()
        except:
            log.exception("While trying to read from database:")
            return None

        table = self.PStoSDSS(table)
        log.debug("Reference Catalog has  %d entries" % len(table))
        return table


if __name__ == '__main__':
    pass
    # gaia = gaiaonline()
    # print (gaia.get_reference_catalog(10,10,0.2))
