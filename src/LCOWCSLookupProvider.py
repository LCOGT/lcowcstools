import numpy as np
from astropy.wcs import Sip

akwcslookup = {
    'ak01': {'SIPA_1_1': 2.981363486451533e-06, 'SIPA_0_2': -1.187772754353179e-05, 'SIPA_2_0': 6.9448752825068934e-06,
             'SIPB_1_1': 1.8171226242992988e-05, 'SIPB_0_2': 3.711915336052473e-06, 'SIPB_2_0': 2.4502180419716353e-06}

}


class akcameraWCSProvider:

    def getWCSForcamera(cameraname, crpix1, crpix2):

        m = 3
        sip_a = np.zeros((m + 1, m + 1), np.double)
        sip_b = np.zeros((m + 1, m + 1), np.double)

        if cameraname in akwcslookup:
            sip_a[1][1] = akwcslookup[cameraname]['SIPA_11']
            sip_a[2][0] = akwcslookup[cameraname]['SIPA_11']
            sip_a[0][2] = akwcslookup[cameraname]['SIPA_11']
            sip_b[1][1] = akwcslookup[cameraname]['SIPA_11']
            sip_b[2][0] = akwcslookup[cameraname]['SIPA_11']
            sip_b[0][2] = akwcslookup[cameraname]['SIPA_11']




        sip = Sip(sip_a, sip_b, None, None, [crpix1,crpix2])
        return sip

