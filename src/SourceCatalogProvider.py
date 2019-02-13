import abc
import requests
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import sep
import logging

from LCOWCSLookupProvider import astrometryServiceRefinceWCS

log = logging.getLogger(__name__)


class SourceCatalogProvider(metaclass=abc.ABCMeta):
    '''Interface to get an source catalog in pixel x/y coordinates out of a FITS image
     '''

    @abc.abstractmethod
    def get_source_catalog(self, imagename) -> (Table, WCS):
        '''

        :param imagename:
        :return: (sourcecatalog, WCS)
        '''
        pass


class e91SourceCatalogProvider(SourceCatalogProvider):
    ''' Read a source catalog and WCS from a LCO level 91 reduced image frame.

    '''

    def get_source_catalog(self, imagename):
        e91image = fits.open(imagename)
        # ra = e91image['SCI'].header['CRVAL1']
        # dec = e91image['SCI'].header['CRVAL2']
        # log.debug("Image is at RA / Dec %f %f " % (ra, dec))
        try:
            sourceCatalog = e91image['CAT'].data
            sourceCatalog['x'] = sourceCatalog['xwin']
            sourceCatalog['y'] = sourceCatalog['ywin']
            log.debug("Source Catalog has %d entries" % len(sourceCatalog))

        except:
            log.warning("%s - No extension \'CAT\' available, skipping." % (e91image))
            e91image.close()
            return (None, None)

        # instantiate the initial guess WCS from the image header
        image_wcs = WCS(e91image['SCI'].header)
        e91image.close()
        return (sourceCatalog, image_wcs)


class SEPSourceCatalogProvider(SourceCatalogProvider):
    ''' Submit fits image to LCO GAIA astrometry service for WCS refinenement and run SEP source finder on image data.
    '''

    def get_source_catalog(self, imagename):

        image_wcs = None

        fitsimage = fits.open(imagename)

        for hdu in fitsimage:
            # search the first valid WCS. Search since we might have a fz compressed file which complicates stuff.
            try:
                original_wcs = WCS(hdu.header)
                continue
            except:
                log.warning("NO RA/DEC found yet, trying next extension")
                original_wcs = None

        log.debug ("FITS file provided WCS is:\n{}".format (original_wcs))

        # Create a source catalog
        # TODO:    Better job of identifying the correct fits extension
        image_data = fitsimage[1].data
        image_data = image_data.astype(float)
        backGround = sep.Background(image_data)
        image_data = image_data - backGround
        backGround = sep.Background(image_data)

        objects = sep.extract(image_data, 5, backGround.globalrms)
        flux_radii, flag = sep.flux_radius(image_data, objects['x'], objects['y'],
                                           6.0 * objects['a'], [0.25, 0.5, 0.75],
                                           normflux=objects['flux'], subpix=5)
        sig = 2.0 / 2.35 * flux_radii[:, 1]
        xwin, ywin, flag = sep.winpos(image_data, objects['x'], objects['y'], sig)
        sourcecatalog = Table([xwin, ywin, objects['flux']], names=['x', 'y', 'flux'])

        log.debug ("Sep found {} sources in image".format (len(sourcecatalog['x'])))

        # Lets refine the WCS solution.
        # TODO: Define condition when we want to refine the WCS
        submitImageInsteadofCatalog = False
        if not 0:
            if submitImageInsteadofCatalog:
                pass
            else:
                log.info ("Sending raw source catalog to astrometry.net service")
                image_wcs = astrometryServiceRefinceWCS (sourcecatalog, original_wcs)
                if image_wcs is None:
                    image_wcs = original_wcs
        else:
            image_wcs = original_wcs

        return sourcecatalog, image_wcs


if __name__ == '__main__':
    # TODO: Make this a test code
    logging.basicConfig(level=getattr(logging, 'DEBUG'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    sourcecatalogProvider = SEPSourceCatalogProvider()

    sourcecatalogProvider.get_source_catalog(
        '/archive/engineering/lsc/ak01/20190107/raw/lsc1m009-ak01-20190107-0484-e00.fits.fz')
