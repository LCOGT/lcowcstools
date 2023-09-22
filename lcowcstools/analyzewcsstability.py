import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import json
import matplotlib.dates as mdates
from wcsfitsdatabase import wcsfitdatabase
__author__ = 'drharbeck@gmail.com'
log = logging.getLogger(__name__)

def dateformat (starttime=None,endtime=None):
    """ Utility to prettify a plot with dates.
    """

    #plt.xlim([starttime, endtime])
    plt.gcf().autofmt_xdate()
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator(bymonth=[2,3,4,5,6,7,8,9,10,11,12])  # every month
    yearsFmt = mdates.DateFormatter('%Y %b')
    monthformat = mdates.DateFormatter('%b')
    plt.gca().xaxis.set_major_locator(years)
    plt.gca().xaxis.set_major_formatter(yearsFmt)
    plt.gca().xaxis.set_minor_locator(months)
    plt.gca().xaxis.set_minor_formatter(monthformat)
    plt.setp(plt.gca().xaxis.get_minorticklabels(), rotation=45)
    plt.setp(plt.gca().xaxis.get_majorticklabels(), rotation=45)
    plt.gca().grid(which='minor')

def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='LCO WCS Tool')


    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    parser.add_argument('--database', default="wcsfits.sqlite")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


def getmeanValue(data, m=2):
    d = np.abs(data - np.median(data))
    mdev = np.median(d)
    s = d/mdev if mdev else 0.
    return np.mean (data[s<m])

def diagnosticByCamera (cameraname, args):

    wcsdb = wcsfitdatabase (args.database)
    allsolutions = wcsdb.readmeasurements(camera=cameraname)
    if not allsolutions:
        log.warning ("No results from db query. skippping")
        return
    cd11 = []
    cd12 = []
    cd21 = []
    cd22 = []

    sip_a20 = []
    sip_a11 = []
    sip_a02 = []
    sip_b20 = []
    sip_b11 = []
    sip_b02 = []

    filterset = set (allsolutions['filter'])
    dateobs = allsolutions['dateobs']
    log.info ("Set of individual filters: {}".format (filterset))

    # sort out the json content
    for wcs in allsolutions['wcs']:
        wcs = json.loads(wcs)
        cd11.append (wcs['CD1_1'])
        cd22.append (wcs['CD2_2'])
        cd12.append (wcs['CD1_2'])
        cd21.append (wcs['CD2_1'])
        sip_a20.append (wcs['sipa20'])
        sip_a11.append (wcs['sipa11'])
        sip_a02.append (wcs['sipa02'])
        sip_b20.append (wcs['sipb20'])
        sip_b11.append (wcs['sipb11'])
        sip_b02.append (wcs['sipb02'])
    cd11 = np.asarray (cd11)
    cd12 = np.asarray (cd12)
    cd21 = np.asarray (cd21)
    cd22 = np.asarray (cd22)
    sip_a20 = np.asarray (sip_a20)
    sip_a11 = np.asarray (sip_a11)
    sip_a02 = np.asarray (sip_a02)
    sip_b20 = np.asarray (sip_b20)
    sip_b11 = np.asarray (sip_b11)
    sip_b02 = np.asarray (sip_b02)

    plt.plot (dateobs, cd11, 'o', label="CD11_{}".format (filter))
    plt.plot (dateobs, cd22, 'o', label="CD22_{}".format (filter))
    plt.legend()
    dateformat()
    plt.savefig ("wcstrend_cddiag_{}.png".format(cameraname))
    plt.close()


    plt.plot (dateobs, cd12, 'o',label="CD12_{}".format (filter))
    plt.plot (dateobs, cd21, 'o',label="CD21_{}".format (filter))
    dateformat()
    plt.legend()
    plt.savefig ("wcstrend_cdnondiag_{}.png".format (cameraname))
    plt.close()






    fig = plt.figure(figsize=(10, 8))
    plt.subplot(211)
    meana11 = getmeanValue(sip_a11)
    meana02 = getmeanValue(sip_a02)
    meana20 = getmeanValue(sip_a20)
    plt.plot (dateobs, sip_a11, 'o',label="A_1_1 {: 6e}".format (meana11))
    plt.plot (dateobs, sip_a20, 'o',label="A_2_0 {: 6e}".format (meana20))
    plt.plot (dateobs, sip_a02, 'o',label="A_0_2 {: 6e}".format (meana02))
    plt.axhline(meana11)
    plt.axhline(meana20)
    plt.axhline(meana02)

    lgd1=plt.legend(bbox_to_anchor=(0,01.02,1,0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=3)
    dateformat()
    plt.ylabel ("SIP A coefficient")
    plt.title ("{}\n\n".format(cameraname))

    plt.subplot(212)
    meanb11 = getmeanValue(sip_b11)
    meanb02 = getmeanValue(sip_b02)
    meanb20 = getmeanValue(sip_b20)
    plt.plot (dateobs, sip_b11, 'o',label="B_1_1 {: 6e}".format (meanb11))
    plt.plot (dateobs, sip_b20, 'o',label="B_2_0 {: 6e}".format (meanb20))
    plt.plot (dateobs, sip_b02, 'o',label="B_0_2 {: 6e}".format (meanb02))

    lgd2 = plt.legend(bbox_to_anchor=(0,01.02,1,0.2), loc="lower left",
               mode="expand", borderaxespad=0, ncol=3)
    dateformat()
    plt.axhline(meanb11)
    plt.axhline(meanb20)
    plt.axhline(meanb02)
    plt.ylabel ("SIP B coefficient")

    plt.tight_layout()
    plt.savefig ("wcstrend_sipab_{}.png".format(cameraname),bbox_extra_artists=(lgd1,lgd2))


    sip =  {cameraname: { 'SIPA_1_1' : meana11,
                         'SIPA_0_2' : meana02,
                         'SIPA_2_0' : meana20,
                         'SIPB_1_1' : meanb11,
                         'SIPB_0_2' : meanb02,
                         'SIPB_2_0' : meanb20,
                         }}
    print (sip)

    plt.close()


    plt.subplot (221)
    plt.plot (allsolutions['azimuth'], cd11, '.' , label='CD1_1')
    plt.plot (allsolutions['azimuth'], cd22, '.' , label='CD2_2')
    plt.legend()
    plt.xlabel("Azimuth [deg]")

    plt.subplot (222)
    plt.plot (allsolutions['azimuth'],  cd21, '.', label='CD1_2')
    plt.plot (allsolutions['azimuth'],  cd12, '.', label='CD2_1')
    plt.legend()
    plt.xlabel('Azimuth')

    plt.subplot (223)
    plt.plot (allsolutions['altitude'], cd11, '.' , label='CD1_1')
    plt.plot (allsolutions['altitude'], cd22, '.' , label='CD2_2')
    plt.legend()
    plt.xlabel("Altitude [deg]")

    plt.subplot (224)
    plt.plot (allsolutions['altitude'],  cd21, '.', label='CD1_2')
    plt.plot (allsolutions['altitude'],  cd12, '.', label='CD2_1')
    plt.legend()
    plt.xlabel('Altitude')
    plt.tight_layout()
    plt.savefig("wcstrans_flexurecd_{}".format(cameraname))
    plt.close()




if __name__ == '__main__':
    args = parseCommandLine()
    db = wcsfitdatabase(args.database)
    cameras = db.getcameras()
    db.close()
    for camera in cameras:
        diagnosticByCamera(camera, args)




