import numpy as np
from astropy.wcs import Sip

akwcslookup = {

    'ak01': {'SIPA_1_1': 2.8875384573560257e-06, 'SIPA_0_2': -1.2776642259520679e-05, 'SIPA_2_0': 6.873210426347869e-06, 'SIPB_1_1': 1.8322056773537455e-05, 'SIPB_0_2': 4.49844648740455e-06, 'SIPB_2_0': 2.076956178814459e-06},
    'ak10': {'SIPA_1_1': 8.421872327419632e-07, 'SIPA_0_2': -1.090631828349261e-05, 'SIPA_2_0': 1.2125230440157089e-05,
             'SIPB_1_1': 2.253149965786546e-05, 'SIPB_0_2': 2.716433310572754e-06, 'SIPB_2_0': 1.0398387533028745e-06},

}


def getWCSForcamera(cameraname, crpix1, crpix2):
    m = 3
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
