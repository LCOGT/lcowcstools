import copy
import logging
import os

from astropy.io import fits
from astropy.table import Table
from astropy.wcs import WCS, Sip
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import astropy
import numpy as np
import math
import matplotlib.pyplot as plt
import scipy.optimize as optimize
from random import random

import argparse

from LCOWCSLookupProvider import getWCSForcamera, transformList, astrometryServiceRefinceWCS
from ReferenceCatalogProvider import refcat2, gaiaonline
from SourceCatalogProvider import e91SourceCatalogProvider, SEPSourceCatalogProvider
from wcsfitsdatabase import wcsfitdatabase
log = logging.getLogger(__name__)

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
    def createMatchedCatalogForLCOe91(imagepath, referenceCatalogProvider, matchradius=5, minobjects=1e20, undistort = False):
        ''' Automatically load source catalog from an LCO e91 processed file, fetch a reference catalog, and return
         a matchedcatalog object.'''

        if ('e91.fits' in imagepath):
            sourceCatalogProvider = e91SourceCatalogProvider()
        else:
            sourceCatalogProvider = SEPSourceCatalogProvider()

        sourceCatalog, image_wcs = sourceCatalogProvider.get_source_catalog(imagepath)
        if (sourceCatalog is None) or (image_wcs is None):
            return None

        if len(sourceCatalog['x']) < minobjects:
            log.info ("Not enough stars found in source catalog (%d). %d are required. Skipping this one." % (len(sourceCatalog['x']), minobjects))
            return None

        ra = image_wcs.wcs.crval[0]
        dec = image_wcs.wcs.crval[1]

        # TODO: get camera identifier, date obs, etc
        hdu = fits.open (imagepath)
        exptime=None
        filter = None
        camera = None
        dateobs = None
        azimuth = None
        altitude = None
        if 'EXPTIME' in hdu[0].header:
            exptime = hdu[0].header['EXPTIME']
        if 'EXPTIME' in hdu[1].header:
            exptime = hdu[1].header['EXPTIME']
        if ('FILTER') in hdu[0].header:
            filter = hdu[0].header['FILTER']
        if ('FILTER') in hdu[1].header:
            filter = hdu[1].header['FILTER']
        if 'DATE-OBS' in hdu[0].header:
            dateobs = hdu[0].header['DATE-OBS']
        if 'DATE-OBS' in hdu[1].header:
            dateobs = hdu[1].header['DATE-OBS']
        if 'INSTRUME' in hdu[0].header:
            camera = hdu[0].header['INSTRUME']
        if 'INSTRUME' in hdu[1].header:
            camera = hdu[1].header['INSTRUME']
        if 'AZIMUTH' in hdu[0].header:
            azimuth = hdu[0].header['AZIMUTH']
        if 'AZIMUTH' in hdu[1].header:
            azimuth = hdu[1].header['AZIMUTH']
        if 'ALTITUDE' in hdu[0].header:
            altitude = hdu[0].header['ALTITUDE']
        if 'ALTITUDE' in hdu[1].header:
            altitude = hdu[1].header['ALTITUDE']

        hdu.close()

        # remove the distortion from the input catalog.
        if undistort:
            sip =getWCSForcamera (camera, image_wcs.wcs.crpix[0],image_wcs.wcs.crpix[1])
            if sip is not None:
                log.info ("undistorting image")
                u,v = transformList (sourceCatalog['x'], sourceCatalog['y'], sip)
                sourceCatalog['x'] = u
                sourceCatalog['y'] = v
                image_wcs = astrometryServiceRefinceWCS (sourceCatalog, image_wcs)

        # fetch a reference catalog:
        referenceCatalog = referenceCatalogProvider.get_reference_catalog(ra, dec, 0.25)

        matchedCatalog = CatalogMatcher()
        matchedCatalog.matchCatalogs(sourceCatalog, referenceCatalog, image_wcs, matchradius)
        matchedCatalog.exptime = exptime
        matchedCatalog.filter = filter
        matchedCatalog.dateobs = dateobs
        matchedCatalog.camera = camera
        matchedCatalog.altitude = altitude
        matchedCatalog.azimuth = azimuth
        matchedCatalog.azimuth = azimuth
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
        log.info ("MatchCatalogs found % 10i pairs at search radius % 6.3f" % (len(self.matchedCatalog['x']), matchradius))
        return self.matchedCatalog

    def updateWCSandUpdateRMS(self, usewcs=None):
        ''' transform the pixel list with a new wcs and get the distance based merrit function of that sollution.
        Note that when this is called, there should be already a matched  catalog avaiable. '''

        if usewcs is not None:
            self.wcs = usewcs
            # log.debug ("WCS updated for MatchedCatalog")
        else:
            pass
            # log.info ("WCS not updated")

        sourcera, sourcedec = self.wcs.all_pix2world(self.matchedCatalog['x'], self.matchedCatalog['y'], 1)

        sourceSkyCoords = SkyCoord(ra=sourcera * u.degree, dec=sourcedec * u.degree)
        referenceSkyCoords = SkyCoord(ra=self.matchedCatalog['RA'] * u.degree,
                                      dec=self.matchedCatalog['Dec'] * u.degree)

        self.matchedCatalog['distarcsec'] = referenceSkyCoords.separation(sourceSkyCoords).arcsecond

        result = math.sqrt(np.sum(self.matchedCatalog['distarcsec'] ** 2) / len(self.matchedCatalog['distarcsec']))
        # log.info ("WCS CRVAL % 12.9f % 12.9f , Source RA / Dec [0] %f %f  Merrit %f" % (self.wcs.wcs.crval[0], self.wcs.wcs.crval[1], sourcera[0], sourcedec[0],  result))

        return result

    def diagnosticPlots(self, basename):
        ''' Generate some helpful diagnostics for the distortion.
        '''

        sourcera, sourcedec = self.wcs.all_pix2world(self.matchedCatalog['x'], self.matchedCatalog['y'], 1)

        deccor = math.cos(self.wcs.wcs.crval[1] * math.pi / 180)

        plt.subplot(projection=self.wcs)
        plt.plot(sourcera, sourcedec, '.')
        plt.plot(self.matchedCatalog['RA'], self.matchedCatalog['Dec'], '.')
        plt.xlabel("RA")
        plt.ylabel("DEC")
        plt.title (basename)
        plt.savefig("%s_RADEC.png" % basename)
        plt.close()

        plt.clf()
        plt.subplot(4, 1, 1)
        plt.title (basename)

        plt.plot(self.matchedCatalog['x'] - self.wcs.wcs.crpix[0],
                 (self.matchedCatalog['RA'] - sourcera) * 3600. / deccor, '.')
        plt.xlabel("X [pixels]")
        plt.ylabel("residual RA [\'\']")
        plt.ylim([-1.75, 1.75])

        plt.subplot(4, 1, 2)
        plt.plot(self.matchedCatalog['x'] - self.wcs.wcs.crpix[0], (self.matchedCatalog['Dec'] - sourcedec) * 3600.,
                 '.')
        plt.xlabel("X [pixels]")
        plt.ylabel("resiudal Dec [\'\']")
        plt.ylim([-1.75, 1.75])

        plt.subplot(4, 1, 3)
        plt.plot(self.matchedCatalog['y'] - self.wcs.wcs.crpix[1],
                 (self.matchedCatalog['RA'] - sourcera) * 3600. / deccor, '.')
        plt.xlabel("Y [pixels]")
        plt.ylabel("residual ra [\'\']")
        plt.ylim([-1.75, 1.75])

        plt.subplot(4, 1, 4)
        plt.plot(self.matchedCatalog['y'] - self.wcs.wcs.crpix[1], (self.matchedCatalog['Dec'] - sourcedec) * 3600.,
                 '.')
        plt.xlabel("Y [pixels]")
        plt.ylabel("residual dec [\'\']")
        plt.ylim([-1.75, 1.75])
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

    def __init__(self, newMatchedCatalog, maxorder=2):
        self.matchedCatalog = newMatchedCatalog
        if self.matchedCatalog.wcs is None:
            log.error("Cannot proceed without a wcs in matched catalog. aborting.")
            return None

        # bootstrap the initial SIP wcs
        self.maxorder = maxorder
        crval = self.matchedCatalog.wcs.wcs.crval
        cd = self.matchedCatalog.wcs.wcs.cd

        self.wcsfitparams = [crval[0], crval[1], cd[0][0], cd[0][1], cd[1][0], cd[1][1]]

        if maxorder > 1:
            self.wcsfitparams.extend([0, 0, 0, 0, 0, 0])

        if maxorder > 2:
            self.wcsfitparams.extend([0, 0, 0, 0])

        merrit = SIPOptimizer.merritFunction(self.wcsfitparams, self.matchedCatalog)
        #log.info("SIPOptimizer init: merrit function is %12.7f" % (merrit))

    @staticmethod
    def merritFunction(sipcoefficients, matchedCatalog):

        matchedCatalog.wcs.wcs.crval[0] = sipcoefficients[0]
        matchedCatalog.wcs.wcs.crval[1] = sipcoefficients[1]
        matchedCatalog.wcs.wcs.cd[0][0] = sipcoefficients[2]
        matchedCatalog.wcs.wcs.cd[0][1] = sipcoefficients[3]
        matchedCatalog.wcs.wcs.cd[1][0] = sipcoefficients[4]
        matchedCatalog.wcs.wcs.cd[1][1] = sipcoefficients[5]

        m = 3
        sip_a = np.zeros((m + 1, m + 1), np.double)
        sip_b = np.zeros((m + 1, m + 1), np.double)

        if len(sipcoefficients) > 6:
            sip_a[1][1] = sipcoefficients[6]
            sip_a[2][0] = sipcoefficients[7]
            sip_a[0][2] = sipcoefficients[8]
            sip_b[1][1] = sipcoefficients[9]
            sip_b[2][0] = sipcoefficients[10]
            sip_b[0][2] = sipcoefficients[11]

        if len(sipcoefficients) > 12:
            sip_a[3][0] = sipcoefficients[12]
            sip_a[0][3] = sipcoefficients[13]
            sip_b[3][0] = sipcoefficients[14]
            sip_b[0][3] = sipcoefficients[15]

        sip = Sip(sip_a, sip_b,
                  None, None, matchedCatalog.wcs.wcs.crpix)
        matchedCatalog.wcs.sip = sip

        matchedCatalog.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        matchedCatalog.wcs.wcs.set()
        merrit = matchedCatalog.updateWCSandUpdateRMS(matchedCatalog.wcs)
        # log.debug("% 12.9f % 12.9f f % 12.7f" % (sipcoefficients[0], sipcoefficients[1], merrit))
        return merrit

    def improveSIP(self):

        deltas = [
            0.001, 0.001,  # CRVAL
            1e-5, 1e-5, 1e-5, 1e-5,  # CD Matrix
            1e-8, 1e-8, 1e-8,  # SIP order 2 A
            1e-8, 1e-8, 1e-8,  # SIP order 2
            1e-10, 1e-10, 1e-10, 1e-10,  # SIP ORDER 3
        ]

        deltas = deltas[:len(self.wcsfitparams)]
        bestfit = optimize.minimize(SIPOptimizer.merritFunction, self.wcsfitparams, args=(self.matchedCatalog))
        merrit = SIPOptimizer.merritFunction(bestfit.x, self.matchedCatalog)
        log.debug("Optimizer return        % 10.4f" % merrit)


def iterativelyFitWCSmany(images, args,  refcat=None):
    if refcat is None:
        refcat = refcat2(args.refcat2)

    if len(images) > 0:
        for image in images:
            iterativelyFitWCSsingle(image, args, searchradii=args.searchradii, refcat=refcat)


def iterativelyFitWCSsingle(image, args, searchradii=[10, 10, 2, 1.5, 1], refcat=None):

    log.info ("Starting to process {}".format (image))
    if refcat is None:
        refcat = refcat2(args.refcat2)

    pngbasename = os.path.basename(image)
    if args.database:
        wcsdb = wcsfitdatabase(args.database)
        if wcsdb.checkifalreadyused(pngbasename) and  not args.reprocess:
            log.info("File already measured. Not doing the wame work twice; skipping")
            wcsdb.close()
            return

    matchedCatalog = CatalogMatcher.createMatchedCatalogForLCOe91(
        image,
        refcat, searchradii[0], minobjects=args.minmatched, undistort=args.undistort)

    if matchedCatalog is None:
        log.info("returned empty catalog, not continuing.")
        return

    if args.makepng:
        matchedCatalog.diagnosticPlots('%s_prefit' % pngbasename)

    # Preserve the
    initialPointing = copy.deepcopy (matchedCatalog.wcs.wcs.crval)

    # do a full fit
    if len(matchedCatalog.matchedCatalog['x']) < args.minmatched:
        log.warning("Not enough stars in input catalog: %d found, %d are required to start. Giving up" % (
        len(matchedCatalog.matchedCatalog['x']), args.minmatched))
        return

    opt = SIPOptimizer(matchedCatalog, maxorder=args.fitorder)
    opt.improveSIP()

    for searchradius in searchradii[1:]:
        matchedCatalog.matchCatalogs(matchradius=searchradius)
        opt = SIPOptimizer(matchedCatalog, maxorder=args.fitorder)
        opt.improveSIP()

    if args.makepng:
        matchedCatalog.diagnosticPlots('%s_postfits' % pngbasename)

    if (args.database):
        wcsdb.addmeasurement(pngbasename, matchedCatalog.dateobs, matchedCatalog.camera, matchedCatalog.filter, None, None, matchedCatalog.azimuth, matchedCatalog.altitude,
                            wcsdb.wcstojson(matchedCatalog.wcs))
    log.info(matchedCatalog.wcs)

    # Final report
    finalPointing = matchedCatalog.wcs.wcs.crval
    log.info ('Fitting updated pointing at CRPIX by {}"'.format ( (finalPointing-initialPointing)*3600))
    if (args.database):
        wcsdb.close()




def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='LCO WCS Tool')

    parser.add_argument('--inputfiles', type=str, nargs='+', help="FITS file for which to derive the WCS function.")
    parser.add_argument('--refcat2', type=str, default='/nfs/AstroCatalogs/Atlas-refcat2/refcat2.db',
                        help='Location of Atlas refcat2 catalog in slite forrmat')
    parser.add_argument('--minmatched', type=int, default=50,
                        help='Minimum number of matched stars to accept solution or even proceed to fit.')
    parser.add_argument('--fitorder', type=int, default=2)
    parser.add_argument('--makepng', action='store_true', help="Create a png output of wcs before and after fit.")
    parser.add_argument('--undistort', action='store_true', help="Undistort input coordiante catalog before WCS fitting.")
    parser.add_argument('--reprocess', action='store_true', help="Reprocess even though file may have been processed already.")
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    parser.add_argument('--searchradii', type=float, nargs='+', default=[10,10,5,3,2])
    parser.add_argument('--database', default="wcsfits.sqlite")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


if __name__ == '__main__':
    args = parseCommandLine()
    log.info("Processing %d input files" % len(args.inputfiles))

    iterativelyFitWCSmany(args.inputfiles, args)
