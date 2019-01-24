import logging
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from random import random

from src.ReferenceCatalogProvider import refcat2, gaiaonline

#from src.ReferenceCatalogProvider import gaiaonline

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

    @staticmethod
    def createMatchedCatalogForLCOe91(imagepath, referenceCatalogProvider, matchradius=5):
        ''' Automatically load source catalog from an LCO e91 processed file, fetch a reference catalog, and return
         a matchedcatalog object.'''

        e91image = fits.open(imagepath)
        ra = e91image['SCI'].header['CRVAL1']
        dec = e91image['SCI'].header['CRVAL2']
        log.debug("Image is at RA / Dec %f %f " % (ra, dec))
        try:
            sourceCatalog = e91image['CAT'].data
            log.debug("Source Catalog has %d entries" % len(sourceCatalog))

        except:
            log.warning("%s - No extension \'CAT\' available, skipping." % (e91image))
            e91image.close()
            return None
        # instanciate the initial guess WCS from the image header

        image_wcs = WCS(e91image['SCI'].header)
        e91image.close()

        # fetch a reference catalog:
        referenceCatalog = referenceCatalogProvider.get_reference_catalog(ra, dec, 0.25)

        matchedCatalog = CatalogMatcher()
        matchedCatalog.matchCatalogs(sourceCatalog, referenceCatalog, image_wcs, matchradius)
        return matchedCatalog

    def matchCatalogs(self, source=None, reference=None, wcs=None, matchradius=5):
        ''' match input catalogs.
        If no new catalogs are given, the match will be done on the chached catalogs of the class.
        '''

        retCatalog = None

        # Cache management
        if wcs is not None:
            self.wcs = wcs

        if source is not None:
            self.source = source

        if reference is not None:
            self.reference = reference

        # transform source catalog to RADEC
        try:
            sourcera, sourcedec = self.wcs.all_pix2world(self.source['x'], self.source['y'], 1)
            sourceSkyCoords = SkyCoord(ra=sourcera * u.degree, dec=sourcedec * u.degree)

            referenceSkyCoords = SkyCoord(ra=self.reference['RA'] * u.degree, dec=self.reference['Dec'] * u.degree)

            idx, d2d, d3d = referenceSkyCoords.match_to_catalog_sky(sourceSkyCoords)
            distance = referenceSkyCoords.separation(sourceSkyCoords[idx]).arcsecond

            matchcondition = (distance < matchradius)
            self.matchedCatalog = Table([self.source['x'][idx][matchcondition],
                                self.source['y'][idx][matchcondition],
                                self.reference['RA'][matchcondition],
                                self.reference['Dec'][matchcondition],
                                distance[matchcondition]],
                               names=['x', 'y', 'RA', 'Dec', 'distarcsec']
                               )

            log.debug("Catalog matched %d entries" % len(self.matchedCatalog))

        except:
            log.exception("Error while transforming and matching")

        return self.matchedCatalog

    def updateWCSandUpdateRMS(self, usewcs=None):
        ''' transform the pixel list with a new wcs and get the distance based merrit function of that sollution.
        Note that when this is called, there should be already a matched  catalog avaiable. '''

        if usewcs is not None:
            self.wcs = usewcs
            #log.debug ("WCS updated for MatchedCatalog")
        else:
            pass
            #log.info ("WCS not updated")

        sourcera, sourcedec = self.wcs.all_pix2world(self.matchedCatalog['x'], self.matchedCatalog['y'], 1)

        sourceSkyCoords = SkyCoord(ra=sourcera * u.degree, dec=sourcedec * u.degree)
        referenceSkyCoords = SkyCoord(ra=self.matchedCatalog['RA'] * u.degree,
                                      dec=self.matchedCatalog['Dec'] * u.degree)

        self.matchedCatalog['distarcsec'] = referenceSkyCoords.separation(sourceSkyCoords).arcsecond

        result = np.sqrt (np.sum (self.matchedCatalog['distarcsec']**2)  / len( self.matchedCatalog['distarcsec']))
        #log.info ("WCS CRVAL % 12.9f % 12.9f , Source RA / Dec [0] %f %f  Merrit %f" % (self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], sourcera[0], sourcedec[0],  result))

        return result

    def diagnosticPlots(self, basename):
        ''' Generate some helpful diagnostics for the distortion.
        '''

        sourcera, sourcedec = self.wcs.all_pix2world(self.matchedCatalog['x'], self.matchedCatalog['y'], 1)

        #plt.subplot(projection=self.wcs)
        plt.plot(sourcera, sourcedec, '.')
        plt.plot(self.matchedCatalog['RA'], self.matchedCatalog['Dec'], '.')
        plt.xlabel ("RA")
        plt.ylabel ("DEC")
        plt.savefig("%s_RADEC.png" % basename)
        plt.close()

        plt.clf()
        plt.subplot (4,1,1)
        plt.plot(self.matchedCatalog['x'] - self.wcs.wcs.crpix[0], (self.matchedCatalog['RA'] - sourcera)*3600., '.')
        plt.xlabel("X [pixels]")
        plt.ylabel("residual RA [\'\']")
        plt.ylim([-0.75,0.75])


        plt.subplot (4,1,2)
        plt.plot(self.matchedCatalog['x'] - self.wcs.wcs.crpix[0], (self.matchedCatalog['Dec'] - sourcedec)*3600., '.')
        plt.xlabel("X [pixels]")
        plt.ylabel("resiudal Dec [\'\']")
        plt.ylim([-0.75,0.75])


        plt.subplot (4,1,3)
        plt.plot(self.matchedCatalog['y'] - self.wcs.wcs.crpix[1], (self.matchedCatalog['RA'] - sourcera)*3600., '.')
        plt.xlabel("Y [pixels]")
        plt.ylabel("residual ra [\'\']")
        plt.ylim([-0.75,0.75])


        plt.subplot (4,1,4)
        plt.plot(self.matchedCatalog['y'] - self.wcs.wcs.crpix[1], (self.matchedCatalog['Dec'] - sourcedec)*3600., '.')
        plt.xlabel("Y [pixels]")
        plt.ylabel("residual dec [\'\']")
        plt.ylim([-0.75,0.75])


        plt.savefig("%s_residuals.png" % basename, dpi=200)
        plt.close()
        # plt.clf()
        # plt.plot(np.sqrt((self.matchedCatalog['y'] - self.wcs.wcs.crpix[1]) ** 2 + (
        #             self.matchedCatalog['x'] - self.wcs.wcs.crpix[0]) ** 2),
        #          self.matchedCatalog['distarcsec'], '.')
        # plt.xlabel("radius [pixels]")
        # plt.ylabel("Distance [\'\']")
        # plt.savefig("%s_radialdist.png" % basename)


class SIPOptimizer:

    def __init__(self, newMatchedCatalog, maxorder):
        self.matchedCatalog = newMatchedCatalog
        if self.matchedCatalog.wcs is None:
            log.error("Cannot proceed without a wcs in matched catalog. aborting.")
            return None

        # bootstrap the initial SIP wcs
        self.maxorder = maxorder
        log.info (self.matchedCatalog.wcs.wcs.ctype)
        crval = self.matchedCatalog.wcs.wcs.crval
        cd = self.matchedCatalog.wcs.wcs.cd

        self.wcsfitparams = [crval[0] , crval[1], cd[0][0], cd[0][1], cd[1][0], cd[1][1]]
        merrit = SIPOptimizer.merritFunction(self.wcsfitparams, self.matchedCatalog)
        log.info("SIPOptimizer init: merrit function is %12.7f" % (merrit))


    @staticmethod
    def merritFunction(sipcoefficients, matchedCatalog):

        matchedCatalog.wcs.wcs.crval[0] = sipcoefficients[0]
        matchedCatalog.wcs.wcs.crval[1] = sipcoefficients[1]
        matchedCatalog.wcs.wcs.cd[0][0] = sipcoefficients[2]
        matchedCatalog.wcs.wcs.cd[0][1] = sipcoefficients[3]
        matchedCatalog.wcs.wcs.cd[1][0] = sipcoefficients[4]
        matchedCatalog.wcs.wcs.cd[1][1] = sipcoefficients[5]
        matchedCatalog.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        merrit = matchedCatalog.updateWCSandUpdateRMS(matchedCatalog.wcs)
        #log.debug("% 12.9f % 12.9f f % 12.7f" % (sipcoefficients[0], sipcoefficients[1], merrit))
        return merrit

    def improveSIP(self):

        bestfit = optimize.minimize(SIPOptimizer.merritFunction, self.wcsfitparams, args=(self.matchedCatalog),
                                     )
        log.info("Optimizer return        %s" % bestfit)

        #


if __name__ == '__main__':
    refcat = refcat2('/nfs/AstroCatalogs/Atlas-refcat2/refcat2.db')
    #refcat = gaiaonline()
    matchedCatalog = CatalogMatcher.createMatchedCatalogForLCOe91(
        '/archive/engineering/lsc/fa15/20190122/processed/lsc1m005-fa15-20190122-0323-e91.fits.fz',
        refcat, 1)

    opt = SIPOptimizer(matchedCatalog, 10)
    matchedCatalog.diagnosticPlots('test_prefit')
    opt.improveSIP()


    matchedCatalog.matchCatalogs(matchradius=0.3)
    opt = SIPOptimizer(matchedCatalog, 10)
    opt.improveSIP()
    matchedCatalog.diagnosticPlots('test_postiteration1')


