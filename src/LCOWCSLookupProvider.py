import logging
import numpy as np
from astropy.wcs import Sip

__author__ = 'drharbeck@gmail.com'

log = logging.getLogger(__name__)

akwcslookup = {

    'ak01': {'SIPA_1_1': 2.8875384573560257e-06, 'SIPA_0_2': -1.2776642259520679e-05, 'SIPA_2_0': 6.873210426347869e-06,
             'SIPB_1_1': 1.8322056773537455e-05, 'SIPB_0_2': 4.49844648740455e-06, 'SIPB_2_0': 2.076956178814459e-06},

    'ak05': {'SIPA_1_1': -6.034337103080695e-06, 'SIPA_0_2': -1.1276384587835203e-05,
             'SIPA_2_0': 1.2802411758347013e-05,
             'SIPB_1_1': 2.4116031017471104e-05, 'SIPB_0_2': -4.302508968114081e-06, 'SIPB_2_0': 5.78043080411828e-07},

    'ak06': {'SIPA_1_1': -9.220745938359626e-07, 'SIPA_0_2': -1.085016354858831e-05, 'SIPA_2_0': 1.3092279633304746e-05,
             'SIPB_1_1': 2.311880941128194e-05, 'SIPB_0_2': -1.5497196904810146e-06,
             'SIPB_2_0': -3.348831288004136e-07},

    'ak10': {'SIPA_1_1': 8.681404689671723e-07, 'SIPA_0_2': -1.0897815933209444e-05, 'SIPA_2_0': 1.2154989981759681e-05,
             'SIPB_1_1': 2.246038112054045e-05, 'SIPB_0_2': 2.698442176640065e-06, 'SIPB_2_0': 9.98348136299767e-07},


    'kb38': {'SIPA_1_1': -5.706068042156633e-05, 'SIPA_0_2': 5.239102464660156e-06, 'SIPA_2_0': -3.0336923628038486e-06,
             'SIPB_1_1': -5.785784217724779e-06, 'SIPB_0_2': -3.141027096614327e-05, 'SIPB_2_0': 2.8867430966128886e-05},


    'kb42': {'SIPA_1_1': -6.223288587920849e-05, 'SIPA_0_2': 4.694665020791541e-06, 'SIPA_2_0': -6.408923602681285e-06,
             'SIPB_1_1': -8.713441942086379e-06, 'SIPB_0_2': -3.306473868683677e-05, 'SIPB_2_0': 3.0916991418377915e-05},

    'ak14': {'SIPA_1_1': 2.579599092541112e-06, 'SIPA_0_2': -1.15833428738132e-05, 'SIPA_2_0': 1.1950471404987229e-05,
              'SIPB_1_1': 2.37986331949568e-05, 'SIPB_0_2': 3.8672274815249445e-06, 'SIPB_2_0': 5.071715952241676e-07},


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


def transformList(x, y, sip):
    """

    :param sip:
    :param x: array of x pixel coordinates
    :param y: array of y pixel coordiantes
    :return: (u,v) transformed pixels, but with CRPIX reapplied.
    """

    log.debug("undistorting source catalog with sip %s" % sip.crpix)
    uv = sip.pix2foc(np.asarray([x, y]).T, 1)
    u = uv[:, 0] + sip.crpix[0]
    v = uv[:, 1] + sip.crpix[1]
    return u, v


