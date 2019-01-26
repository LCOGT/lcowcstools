import abc

import requests
import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.table import Table
import sep

import logging

log = logging.getLogger(__name__)
logging.basicConfig(level=getattr(logging, 'DEBUG'),
                    format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')


class SourceCatalogProvider(metaclass=abc.ABCMeta):

    @abc.abstractmethod
    def get_source_catalog(self, imagename):
        '''

        :param imagename:
        :return: (sourcecatalog, WCS)
        '''
        pass


class e91SourceCatalogProvider(SourceCatalogProvider):

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
        # instanciate the initial guess WCS from the image header

        image_wcs = WCS(e91image['SCI'].header)
        e91image.close()

        return (sourceCatalog, image_wcs)


class blindGaiaAstrometrySourceCatalogProvider(SourceCatalogProvider):
    url = 'http://astrometry.lco.gtn/image/'

    def get_source_catalog(self, imagename):
        fitsimage = fits.open(imagename)
        ra = None
        image_wcs = None
        for hdu in fitsimage:

            try:
                ra = hdu.header['CRVAL1']
                dec = hdu.header['CRVAL2']
                image_wcs = WCS(hdu.header)

                continue
            except:
                log.debug("NO RA/DEC found yet, trying next extension")

        log.debug("RA/Dec of image is: %s %s" % (ra, dec))

        # No WCS? Get us one! Or always refine it :-)
        if not 0:
            payload = {'ra': ra, 'dec': dec, 'image_path': imagename}
            response = requests.post(self.url, json=payload)
            response = response.json()
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

        # We now have a good first order wcs, let's find all thje sources in the image.

        # sourceRA = response['matched_sources']['ra']
        # sourceDec = response['matched_sources']['dec']
        # sourcex, sourcey =  image_wcs.all_world2pix(sourceRA, sourceDec, 1)

        image_data = fitsimage[1].data
        image_data = image_data.astype(float)

        backGround = sep.Background(image_data)
        image_data = image_data - backGround
        log.debug("Background: %f \pm %f " % (backGround.globalback, backGround.globalrms))

        objects = sep.extract(image_data, 10, backGround.globalrms)
        flux_radii, flag = sep.flux_radius(image_data, objects['x'], objects['y'],
                                           6.0 * objects['a'], [0.25, 0.5, 0.75],
                                           normflux=objects['flux'], subpix=5)


        sig = 2.0 / 2.35 * flux_radii[:, 1]
        xwin, ywin, flag = sep.winpos(image_data, objects['x'], objects['y'], sig)

        retTable = Table([xwin,ywin], names=['x', 'y'])
        return retTable, image_wcs


if __name__ == '__main__':
    sourcecatalogProvider = blindGaiaAstrometrySourceCatalogProvider()

    sourcecatalogProvider.get_source_catalog(
        '/archive/engineering/tlv/ak10/20190124/raw/tlv1m0XX-ak10-20190124-0012-e00.fits.fz')
