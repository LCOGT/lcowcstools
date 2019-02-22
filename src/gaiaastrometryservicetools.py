import argparse
import logging
import numpy as np
import requests
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
__author__ = 'drharbeck@gmail.com'



log = logging.getLogger(__name__)

LCO_GAIA_ASTROMETRY_URL = 'http://astrometry.lco.gtn/'


def astrometryServiceRefineWCSFromCatalog(sourceCatalog, wcs):
    originalPointing = np.asarray([wcs.wcs.crval[0], wcs.wcs.crval[1]])

    payload = {'ra': wcs.wcs.crval[0], 'dec': wcs.wcs.crval[1],
               'crpix1': wcs.wcs.crpix[0], 'crpix2': wcs.wcs.crpix[1],
               'pixel_scale': np.mean(proj_plane_pixel_scales(wcs)) * 3600.,
               'naxis1': int(np.max(sourceCatalog['x'])), 'naxis2': int(np.max(sourceCatalog['y'])), 'naxis': 2,
               'X': list(i for i in sourceCatalog['x']),
               'Y': list(i for i in sourceCatalog['y']),
               'FLUX': list(i for i in sourceCatalog['flux'])
               }

    try:
        response = requests.post("{}/catalog/".format(LCO_GAIA_ASTROMETRY_URL), json=payload)
        response = response.json()
        log.debug(response)
    except:
        log.error("Error while executing astrometry.net service with payload:\n %s\n\nGot Response %s" % payload,
                  response)
        return wcs

    image_wcs = astrometryServicereadWCSFromResponse(response, originalPointing)
    return image_wcs


def astrometryServicereadWCSFromResponse(response, originalPointing=None):
    if response is None:
        return None

    image_wcs = None
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
        log.debug("Updated Gaia-service wcs:\n {}".format(image_wcs))
        newPointing = np.asarray([image_wcs.wcs.crval[0], image_wcs.wcs.crval[1]])
        if originalPointing is not None:
            log.info(
                'Gaia service updated image pointing at CRPIX by {}"'.format((newPointing - originalPointing) * 3600))
    else:
        log.warning("Astrometry.net could not find a solution!")
    return image_wcs


def astrometryServiceRefineWCSFromImage(imagepath, wcs=None):
    log.debug(imagepath)
    if wcs is None:
        fitsimage = fits.open(imagepath)
        for hdu in fitsimage:

            try:
                wcs = WCS(hdu.header)
                fitsimage.close()
                continue
            except:
                log.warning("NO RA/DEC found yet, trying next extension")
                wcs = None
    if wcs is None:
        log.error("No WCS found in image header. skipping")
        return None

    originalPointing = None if wcs is None else np.asarray([wcs.wcs.crval[0], wcs.wcs.crval[1]])

    payload = {'ra': wcs.wcs.crval[0], 'dec': wcs.wcs.crval[1],
               'image_path': imagepath,
               }
    try:
        response = requests.post("{}/image/".format(LCO_GAIA_ASTROMETRY_URL), json=payload)
        response = response.json()
        log.debug(response)
    except:
        log.error("Error while executing astrometry.net service with payload:\n %s\n\nGot Response %s" % payload,
                  response)
        return wcs

    image_wcs = astrometryServicereadWCSFromResponse(response, originalPointing)
    return image_wcs


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='LCO WCS Tools - gaiaastrometryservicetools')

    parser.add_argument('inputfile', type=str, nargs=1, help="FITS file for which to derive the WCS function.")

    parser.add_argument('--loglevel', dest='log_level', default='DEBUG', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')

    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


if __name__ == '__main__':
    args = parseCommandLine()
    for file in args.inputfile:
        astrometryServiceRefineWCSFromImage(file, None)
