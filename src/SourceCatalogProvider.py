import abc
import requests
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import sep
import logging

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
        ra = e91image['SCI'].header['CRVAL1']
        dec = e91image['SCI'].header['CRVAL2']
        log.debug("Image is at RA / Dec %f %f " % (ra, dec))
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


class blindGaiaAstrometrySourceCatalogProvider(SourceCatalogProvider):
    ''' Submit fits image to LCO GAIA astrometry service for WCS refinenement and run SEP source finder on image data.
    '''

    # lco astrometry service URL
    url = 'http://astrometry.lco.gtn/'

    def get_source_catalog(self, imagename):
        ra = None
        image_wcs = None
        fitsimage = fits.open(imagename)
        pixscale= None
        naxis1 = None
        naxis2 = None


        for hdu in fitsimage:
            # search the first valid WCS. Search since we might have a fz compressed file which complicates stuff.
            try:
                ra = hdu.header['CRVAL1']
                dec = hdu.header['CRVAL2']
                crpix1=hdu.header['CRPIX1']
                crpix2=hdu.header['CRPIX2']
                if ('ZNAXIS') in hdu.header:
                    naxis1=hdu.header['ZNAXIS1']
                    naxis2=hdu.header['ZNAXIS2']
                else:
                    naxis1=hdu.header['NAXIS1']
                    naxis2=hdu.header['NAXIS2']
                pixscale = hdu.header['PIXSCALE']
                image_wcs = WCS(hdu.header)
                continue
            except:
                log.warning("NO RA/DEC found yet, trying next extension")

        log.info("RA/Dec of image is: %s %s" % (ra, dec))

        # Get a source catalog
        # TODO:    Better job of identifying the correct fits extension
        image_data = fitsimage[1].data
        image_data = image_data.astype(float)
        backGround = sep.Background(image_data)
        image_data = image_data - backGround
        backGround = sep.Background(image_data)
        log.debug("Background: %f \pm %f " % (backGround.globalback, backGround.globalrms))

        objects = sep.extract(image_data, 5, backGround.globalrms)
        flux_radii, flag = sep.flux_radius(image_data, objects['x'], objects['y'],
                                           6.0 * objects['a'], [0.25, 0.5, 0.75],
                                           normflux=objects['flux'], subpix=5)
        sig = 2.0 / 2.35 * flux_radii[:, 1]
        xwin, ywin, flag = sep.winpos(image_data, objects['x'], objects['y'], sig)
        sourcecatalog = Table([xwin, ywin], names=['x', 'y'])

        # Lets refine the WCS solution.
        # TODO: Define condition when we want to refine the WCS
        submitImageInsteadofCatalog = True
        image_wcs = None
        if not 0:
            if submitImageInsteadofCatalog:
                payload = {'ra': ra, 'dec': dec, 'image_path': imagename}
                response = requests.post("{}/image/".format(self.url), json=payload)
            else:

                payload = {'ra': ra, 'dec': dec,
                           'crpix1': crpix1, 'crpix2': crpix2,
                           'pixel_scale': pixscale,
                           'naxis1': naxis1, 'naxis2': naxis2, 'naxis':2,
                           'X': list (i for i in sourcecatalog['x']),
                           'Y': list (i for i in sourcecatalog['y']),
                           'FLUX': list (i for i in objects['flux'])
                    }

                response = requests.post("{}/catalog/".format(self.url), json=payload)
            try:
                response = response.json()
                log.info (response)
            except:
                log.error ("Error while executing astrometry.net service %s" % payload)
                return None, None

            if 'CD1_1' in response:
                # Build a WCS, reverse-transform the source list.
                image_wcs = WCS(naxis=2)
                image_wcs.wcs.crpix = [response['CRPIX1'], response['CRPIX2']]
                image_wcs.wcs.crval = [response['CRVAL1'], response['CRVAL2']]
                image_wcs.wcs.cd = np.zeros((2, 2))
                image_wcs.wcs.cd[0][0] = response['CD1_1']
                image_wcs.wcs.cd[0][1] = response['CD1_2']
                image_wcs.wcs.cd[1][0] = response['CD2_1']
                image_wcs.wcs.cd[1][1] = response['CD2_2']
                image_wcs.wcs.ctype = [response['CTYPE1'], response['CTYPE2']]

        # We now have a good first order wcs, let's find all the sources in the image.

        log.debug (image_wcs)
        return sourcecatalog, image_wcs


if __name__ == '__main__':
    # TODO: Make this a test code
    logging.basicConfig(level=getattr(logging, 'DEBUG'),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    sourcecatalogProvider = blindGaiaAstrometrySourceCatalogProvider()

    sourcecatalogProvider.get_source_catalog(
        '/archive/engineering/lsc/ak01/20190107/raw/lsc1m009-ak01-20190107-0484-e00.fits.fz')
