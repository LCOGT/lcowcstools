import logging
from astropy.table import Table
from astropy.coordinates import SkyCoord
from astropy import units as u
from astropy.io import fits
import numpy as np
import math
import matplotlib.pyplot as plt
from lcowcstools.LCOWCSLookupProvider import getWCSForcamera, transformList
from lcowcstools.gaiaastrometryservicetools import astrometryServiceRefineWCSFromCatalog
from lcowcstools.SourceCatalogProvider import e91SourceCatalogProvider, SEPSourceCatalogProvider

__author__ = 'drharbeck@gmail.com'

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
    def createMatchedCatalogForLCO(imagepath, referenceCatalogProvider, matchradius=5, minobjects=1e20,
                                   undistort=False):
        ''' Automatically load source catalog from an LCO e91 processed file, fetch a reference catalog, and return
         a matchedcatalog object.'''

        if ('e91.fits' in imagepath) or ('s91.fits' in imagepath):
            sourceCatalogProvider = e91SourceCatalogProvider()
        else:
            sourceCatalogProvider = SEPSourceCatalogProvider(refineWCSViaLCO=True)

        sourceCatalog, image_wcs = sourceCatalogProvider.get_source_catalog(imagepath)
        if (sourceCatalog is None) or (image_wcs is None):
            return None

        if len(sourceCatalog['x']) < minobjects:
            log.info("Not enough stars found in source catalog (%d). %d are required. Skipping this one." % (
                len(sourceCatalog['x']), minobjects))
            return None

        ra = image_wcs.wcs.crval[0]
        dec = image_wcs.wcs.crval[1]

        # TODO: get camera identifier, date obs, etc

        exptime = None
        filter = None
        camera = None
        dateobs = None
        azimuth = None
        altitude = None

        hdu = fits.open(imagepath)
        # TODO: We are opening and closing fits files quite a lot here, might be not most efficient.
        # Go searching for meta data, in multiple extension ssince we might have a .fz compressed file :-(
        for extension in [0, 1]:
            if 'EXPTIME' in hdu[extension].header:
                exptime = hdu[extension].header['EXPTIME']

            if ('FILTER') in hdu[extension].header:
                filter = hdu[extension].header['FILTER']

            if 'DATE-OBS' in hdu[extension].header:
                dateobs = hdu[extension].header['DATE-OBS']

            if 'INSTRUME' in hdu[extension].header:
                camera = hdu[extension].header['INSTRUME']

            if 'AZIMUTH' in hdu[extension].header:
                azimuth = hdu[extension].header['AZIMUTH']

            if 'ALTITUDE' in hdu[extension].header:
                altitude = hdu[extension].header['ALTITUDE']
        hdu.close()

        # remove the distortion from the input catalog if requested and refine the WCS.
        if undistort:
            sip = getWCSForcamera(camera, image_wcs.wcs.crpix[0], image_wcs.wcs.crpix[1])
            if sip is not None:
                log.info("undistorting image")
                u, v = transformList(sourceCatalog['x'], sourceCatalog['y'], sip)
                sourceCatalog['x'] = u
                sourceCatalog['y'] = v
                dedistortedwcs = astrometryServiceRefineWCSFromCatalog(sourceCatalog, image_wcs)
                if dedistortedwcs is not None:
                    image_wcs = dedistortedwcs
                else:
                    log.warning("astrometry.net did not find a solution on the undistorted image. Using original wcs")

        # fetch a reference catalog:
        referenceCatalog = referenceCatalogProvider.get_reference_catalog(ra, dec, 0.30)

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

        self.matchedCatalog = None

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
                                         distance[matchcondition]
                                         ],
                                        names=['x', 'y', 'RA', 'Dec', 'distarcsec']
                                        )
        except:
            log.exception("Error while transforming and matching")

        nummatched = len(self.matchedCatalog) if self.matchedCatalog is not None else 0
        log.info("MatchCatalogs found {: 10d} pairs at search radius {: 6.3f}".format(nummatched, matchradius))
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

    def diagnosticPlots(self, basename, resrange=20):
        ''' Generate some helpful diagnostics for the distortion.
        '''
        if not self.matchedCatalog:
            return

        sourcera, sourcedec = self.wcs.all_pix2world(self.matchedCatalog['x'], self.matchedCatalog['y'], 1)

        deccor = math.cos(self.wcs.wcs.crval[1] * math.pi / 180)

        plt.subplot(projection=self.wcs)
        plt.plot(sourcera, sourcedec, '.', markersize=1)
        plt.plot(self.matchedCatalog['RA'], self.matchedCatalog['Dec'], '.', markersize=1)
        plt.xlabel("RA")
        plt.ylabel("DEC")
        plt.title(basename)
        plt.savefig("%s_RADEC.png" % basename)
        plt.close()



        plt.clf()
        fig = plt.gcf()
        fig.set_size_inches(5, 10)


        ax1 = plt.subplot(4, 1, 1)
        plt.title(basename)
        plt.plot(self.matchedCatalog['x'] - self.wcs.wcs.crpix[0],
             (self.matchedCatalog['RA'] - sourcera) * 3600. * deccor, '.', markersize=1)
        plt.xlabel("X [pixels]")
        plt.ylabel("res RA [\'\']")
        plt.ylim([-1*resrange,resrange])


        ax2=plt.subplot(4, 1, 2)
        plt.plot(self.matchedCatalog['x'] - self.wcs.wcs.crpix[0],
                 (self.matchedCatalog['Dec'] - sourcedec) * 3600.,
                 '.', markersize=1)
        plt.xlabel("X [pixels]")
        plt.ylabel("res Dec [\'\']")
        plt.ylim([-1*resrange,resrange])








        plt.subplot(4, 1, 3)
        plt.plot(self.matchedCatalog['y'] - self.wcs.wcs.crpix[1],
                 (self.matchedCatalog['RA'] - sourcera) * 3600. * deccor, '.', markersize=1)
        plt.xlabel("Y [pixels]")
        plt.ylabel("res RA [\'\']")
        plt.ylim([-1*resrange,resrange])

        plt.subplot(4, 1, 4)
        plt.plot(self.matchedCatalog['y'] - self.wcs.wcs.crpix[1],
                (self.matchedCatalog['Dec'] - sourcedec) * 3600., '.', markersize=1)
        plt.xlabel("Y [pixels]")
        plt.ylabel("res Dec [\'\']")
        plt.ylim([-1*resrange,resrange])

        plt.tight_layout()
        plt.savefig("%s_residuals.png" % basename, dpi=300)
        plt.close()
        # plt.clf()
        # plt.plot(np.sqrt((self.matchedCatalog['y'] - self.wcs.wcs.crpix[1]) ** 2 + (
        #             self.matchedCatalog['x'] - self.wcs.wcs.crpix[0]) ** 2),
        #          self.matchedCatalog['distarcsec'], '.')
        # plt.xlabel("radius [pixels]")
        # plt.ylabel("Distance [\'\']")
        # plt.savefig("%s_radialdist.png" % basename)



