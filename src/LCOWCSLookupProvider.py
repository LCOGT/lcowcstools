import json
import logging
import numpy as np
import requests
from astropy.wcs import Sip, WCS
from astropy.wcs.utils import proj_plane_pixel_scales

log = logging.getLogger(__name__)

akwcslookup = {

    'ak01': {'SIPA_1_1': 2.8875384573560257e-06, 'SIPA_0_2': -1.2776642259520679e-05, 'SIPA_2_0': 6.873210426347869e-06,
             'SIPB_1_1': 1.8322056773537455e-05, 'SIPB_0_2': 4.49844648740455e-06, 'SIPB_2_0': 2.076956178814459e-06},

    'ak06': {'SIPA_1_1': -9.220745938359626e-07, 'SIPA_0_2': -1.085016354858831e-05, 'SIPA_2_0': 1.3092279633304746e-05,
             'SIPB_1_1': 2.311880941128194e-05, 'SIPB_0_2': -1.5497196904810146e-06, 'SIPB_2_0': -3.348831288004136e-07},

    'ak10': {'SIPA_1_1': 8.681404689671723e-07, 'SIPA_0_2': -1.0897815933209444e-05, 'SIPA_2_0': 1.2154989981759681e-05,
             'SIPB_1_1': 2.246038112054045e-05, 'SIPB_0_2': 2.698442176640065e-06, 'SIPB_2_0': 9.98348136299767e-07},

}


def getWCSForcamera(cameraname, crpix1, crpix2):
    """ Return  SIP non-linear coordiante correction object intialized for a camera from a lookup table.

    If the camera is not in the lookup table, an identify transformation is returned.

    TODO: variable order, so far limit ouselves to second order
    TODO: Time-constraint lookup.

    :param cameraname: Name of camera, e.g., ak01
    :param crpix1: CRPIX1 for camera, as that my have changed over time
    :param crpix2: CRPIX2 for camera, as that my have changed over time
    :return:
    """
    m = 2
    sip_a = np.zeros((m + 1, m + 1), np.double)
    sip_b = np.zeros((m + 1, m + 1), np.double)

    if cameraname in akwcslookup:
        sip_a[1][1] = akwcslookup[cameraname]['SIPA_1_1']
        sip_a[2][0] = akwcslookup[cameraname]['SIPA_2_0']
        sip_a[0][2] = akwcslookup[cameraname]['SIPA_0_2']
        sip_b[1][1] = akwcslookup[cameraname]['SIPB_1_1']
        sip_b[2][0] = akwcslookup[cameraname]['SIPB_2_0']
        sip_b[0][2] = akwcslookup[cameraname]['SIPB_0_2']

    sip = Sip(sip_a, sip_b, None, None, [crpix1, crpix2])
    return sip


def transformList (x,y, sip):
    """

    :param sip:
    :param x:
    :param y:
    :return: (u,v) transformed pixels, but with CRPIX reapplied
    """

    log.debug ("undistorting image with sip %s" %  sip.crpix)
    uv = sip.pix2foc (np.asarray([x,y]).T,1)
    u = uv[:,0] + sip.crpix[0]
    v = uv[:,1] + sip.crpix[1]
    return u,v



def astrometryServiceRefinceWCS (sourceCatalog, wcs):
    url = 'http://astrometry.lco.gtn/'

    originalPointing = np.asarray([wcs.wcs.crval[0],  wcs.wcs.crval[1]])

    payload = {'ra': wcs.wcs.crval[0], 'dec': wcs.wcs.crval[1],
               'crpix1': wcs.wcs.crpix[0], 'crpix2': wcs.wcs.crpix[1],
               'pixel_scale': np.mean (proj_plane_pixel_scales(wcs)) * 3600.,
               'naxis1': int(np.max(sourceCatalog['x'])), 'naxis2': int(np.max (sourceCatalog['y'])), 'naxis':2,
               'X': list (i for i in sourceCatalog['x']),
               'Y': list (i for i in sourceCatalog['y']),
               'FLUX': list (i for i in sourceCatalog['flux'])
               }

    try:
        response = requests.post("{}/catalog/".format(url), json=payload)
        response = response.json()
        log.debug (response)
    except:
        log.error ("Error while executing astrometry.net service with payload:\n %s\n\nGot Response %s" % payload, response)
        return wcs

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

    log.debug ("Updated Gaia-service wcs:\n {}".format(image_wcs))

    newPointing = np.asarray([image_wcs.wcs.crval[0],  image_wcs.wcs.crval[1]])
    log.info ('Gaia service updated image pointing at CRPIX by {}"'.format ( (newPointing-originalPointing)*3600))
    return image_wcs
