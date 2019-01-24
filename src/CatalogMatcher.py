import logging
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import os
import sqlite3
import numpy as np
import math
import matplotlib.pyplot as plt



log = logging.getLogger(__name__)
logging.basicConfig(level=getattr(logging, 'DEBUG'),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')


class CatalogMatcher:
    '''
        Class to match two input catalogs:
        sourcecatalog is a catalog of sources extracted from an image, in coordinates of pixels (x,y)
        referencecatalog is a catalog of on-sky objects based on existing surveys, in coordinates of (RA, Dec)
        WCS is a astropy world coordiante system.

        the source catalog shall be a astropy Table with the columns 'x', 'y'
        the reference catalog shall be a astropy Table with the columns 'RA', 'Dec'

    '''
    source = None
    reference = None
    matchedCatalog = None

    @staticmethod
    def createMatchedCatalogForLCOe91 (imagepath, refcat2db, matchradius=5):

        e91image = fits.open(imagepath)
        ra = e91image['SCI'].header['CRVAL1']
        dec = e91image['SCI'].header['CRVAL2']

        try:
            sourceCatalog = e91image['CAT'].data
            log.debug ("Source Catalog has %d entries" % len(sourceCatalog))

        except:
            log.warning("%s - No extension \'CAT\' available, skipping." % (e91image))
            e91image
            return None
        # instanciate the initial guess WCS from the image header

        image_wcs = WCS (e91image['SCI'].header)

        # fetch a reference catalog:
        referenceCatalogProvider = refcat2(refcat2db)
        referenceCatalog = referenceCatalogProvider.get_reference_catalog(ra,dec,0.25)

        matchedCatalog = CatalogMatcher()
        matchedCatalog.matchCatalogs(sourceCatalog, referenceCatalog, image_wcs, matchradius)
        return matchedCatalog


    def matchCatalogs (self, source, reference, wcs, matchradius = 5):

        retCatalog = None
        self.wcs = wcs
        if WCS is None:
            pass

        # transform source catalog to RADEC
        try:

            sourcera,sourcedec = wcs.all_pix2world(source['x'], source['y'], 1)
            sourceSkyCoords = SkyCoord (ra =  sourcera * u.degree, dec = sourcedec * u.degree)

            referenceSkyCoords = SkyCoord (ra=reference['RA'] * u.degree, dec = reference['Dec'] * u.degree)

            idx, d2d,d3d = referenceSkyCoords.match_to_catalog_sky(sourceSkyCoords)
            distance = referenceSkyCoords.separation(sourceSkyCoords[idx]).arcsecond

            matchcondition = (distance < matchradius)
            retCatalog = Table ( [source['x'][idx][matchcondition],
                                  source['y'][idx][matchcondition],
                                  reference['RA'][matchcondition],
                                  reference['Dec'][matchcondition],
                                  distance[matchcondition]],
                                 names=['x','y','RA','Dec', 'distarcsec']
            )

            self.matchedCatalog = retCatalog

        except:
            log.exception("Error while transforming and matching")

        return retCatalog


    def diagnosticPlots (self, basename, wcs=None):

        sourcera,sourcedec = self.wcs.all_pix2world(self.matchedCatalog['x'], self.matchedCatalog['y'], 1)

        plt.subplot (projection=self.wcs)
        plt.plot (sourcera,sourcedec, '.')
        plt.plot (self.matchedCatalog['RA'], self.matchedCatalog['Dec'],'.')
        plt.savefig ("%s_RADEC.png" % basename)

        plt.clf()
        plt.plot (self.matchedCatalog['x']-self.wcs.wcs.crpix[0], self.matchedCatalog['distarcsec'],'.')
        plt.xlabel ("X [pixels]")
        plt.ylabel ("Distance [\'\']")
        plt.savefig ("%s_RAdist.png" % basename)


        plt.clf()
        plt.plot (self.matchedCatalog['y']-self.wcs.wcs.crpix[1], self.matchedCatalog['distarcsec'],'.')
        plt.xlabel ("Y [pixels]")
        plt.ylabel ("Distance [\'\']")
        plt.savefig ("%s_Decdist.png" % basename)


        plt.clf()
        plt.plot (np.sqrt ( (self.matchedCatalog['y']-self.wcs.wcs.crpix[1])**2 + (self.matchedCatalog['x']-self.wcs.wcs.crpix[0])**2) ,
              self.matchedCatalog['distarcsec'],'.')
        plt.xlabel ("Y [pixels]")
        plt.ylabel ("Distance [\'\']")
        plt.savefig ("%s_radialdist.png" % basename)



class refcat2:
    # Interface to query consolidat3ed sqlite3 db file that was geernated from Tonry (2018) refcat2


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


    def get_reference_catalog(self, ra, dec, radius, overwrite_select=False):

        rows = None
        try:
            connection = sqlite3.connect(self.dbfile)
            cursor = connection.cursor()
            if (not radius == None and radius > 0):
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
            table = Table (rows, names=['RA','Dec','g','r','i','z'])
            cursor.close()
        except:
            log.exception("While trying to read from database:")
            return None

        table = self.PStoSDSS(table)
        log.debug("Reference Catalog has  %d entries" % len(table))
        return table

if __name__ == '__main__':
    matchedCatalog = CatalogMatcher.createMatchedCatalogForLCOe91('/archive/engineering/lsc/fa15/20190122/processed/lsc1m005-fa15-20190122-0323-e91.fits.fz',
                                                 '/nfs/AstroCatalogs/Atlas-refcat2/refcat2.db', 3)

    matchedCatalog.diagnosticPlots('test')