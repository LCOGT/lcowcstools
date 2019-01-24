import logging
from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS
from astropy.coordinates import SkyCoord
from astropy import units as u
import numpy as np
import math
import matplotlib.pyplot as plt

from src.ReferenceCatalogProvider import refcat2
from src.ReferenceCatalogProvider import gaiaonline


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
    def createMatchedCatalogForLCOe91 (imagepath, referenceCatalogProvider, matchradius=5):
        ''' Automatically load source catalog from an LCO e91 processed file, fetch a reference catalog, and return
         a matchedcatalog object.'''

        e91image = fits.open(imagepath)
        ra = e91image['SCI'].header['CRVAL1']
        dec = e91image['SCI'].header['CRVAL2']
        log.info ("Image is at RA / Dec %f %f " % (ra,dec))
        try:
            sourceCatalog = e91image['CAT'].data
            log.debug ("Source Catalog has %d entries" % len(sourceCatalog))

        except:
            log.warning("%s - No extension \'CAT\' available, skipping." % (e91image))
            e91image.close()
            return None
        # instanciate the initial guess WCS from the image header

        image_wcs = WCS (e91image['SCI'].header)
        e91image.close()

        # fetch a reference catalog:
        #referenceCatalogProvider = refcat2(refcat2db)
        referenceCatalog = referenceCatalogProvider.get_reference_catalog(ra,dec,0.25)

        matchedCatalog = CatalogMatcher()
        matchedCatalog.matchCatalogs(sourceCatalog, referenceCatalog, image_wcs, matchradius)
        return matchedCatalog


    def matchCatalogs (self, source=None, reference=None, wcs=None, matchradius = 5):
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
            sourcera,sourcedec = self.wcs.all_pix2world(self.source['x'], self.source['y'], 1)
            sourceSkyCoords = SkyCoord (ra =  sourcera * u.degree, dec = sourcedec * u.degree)

            referenceSkyCoords = SkyCoord (ra=self.reference['RA'] * u.degree, dec = self.reference['Dec'] * u.degree)

            idx, d2d,d3d = referenceSkyCoords.match_to_catalog_sky(sourceSkyCoords)
            distance = referenceSkyCoords.separation(sourceSkyCoords[idx]).arcsecond

            matchcondition = (distance < matchradius)
            retCatalog = Table ( [self.source['x'][idx][matchcondition],
                                  self.source['y'][idx][matchcondition],
                                  self.reference['RA'][matchcondition],
                                  self.reference['Dec'][matchcondition],
                                  distance[matchcondition]],
                                 names=['x','y','RA','Dec', 'distarcsec']
            )

            self.matchedCatalog = retCatalog
            log.debug ("Catalog matched %d entries" % len(retCatalog))

        except:
            log.exception("Error while transforming and matching")

        return retCatalog



    def getUpdatedRMS (self, usewcs = None):
        ''' transform the pixel list with a new wcs and get the distance based merrit function of that sollution.
        Note that when this is called, there should be already a matched  catalog avaiable. '''

        if usewcs is None:
            usewcs = self.wcs
        sourcera,sourcedec = usewcs.all_pix2world(self.matchedCatalog['x'], self.matchedCatalog['y'], 1)

        sourceSkyCoords = SkyCoord (ra =  sourcera * u.degree, dec = sourcedec * u.degree)
        referenceSkyCoords = SkyCoord (ra=self.matchedCatalog['RA'] * u.degree, dec = self.matchedCatalog['Dec'] * u.degree)
        distance = referenceSkyCoords.separation(sourceSkyCoords).arcsecond

        return  math.sqrt ( np.sum (distance**2) / len(distance) )



    def diagnosticPlots (self, basename):
        ''' Generate some helpful diagnostics for the distortion.
        '''

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
        plt.xlabel ("radius [pixels]")
        plt.ylabel ("Distance [\'\']")
        plt.savefig ("%s_radialdist.png" % basename)






class SIPOptimizer:

    def __init__(self, MatchedCatalog, maxorder):
        self.MatchedCatalog = MatchedCatalog
        self.maxorder = maxorder


    def merritFunction (self, matchedCatalog, sipcoefficients):

        matchedCatalog.wcs.update (sipcoefficients)
        matchedCatalog.matchCatalogs()
        return matchedCatalog.gerrms()

    def improveSIP (self):
        initialGuess = np.zeros[self.maxorder]
        scipy.optimize.minimize (self.merritFunction, initialGuess)




if __name__ == '__main__':
    refcat =  refcat2('/nfs/AstroCatalogs/Atlas-refcat2/refcat2.db')
    #refcat = gaiaonline()
    matchedCatalog = CatalogMatcher.createMatchedCatalogForLCOe91('/archive/engineering/lsc/fa15/20190122/processed/lsc1m005-fa15-20190122-0323-e91.fits.fz',
                                                refcat , 3)

    log.info ("Residual error of matched catalog: % 7.3f" % matchedCatalog.getUpdatedRMS())
    matchedCatalog.diagnosticPlots('test')


    #opt = SIPOptimizer (matchedCatalog,10)
    #opt.improveSIP()