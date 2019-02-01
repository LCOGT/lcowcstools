import argparse
import logging
import matplotlib.pyplot as plt
import numpy as np
import json
import matplotlib.dates as mdates
from wcsfitsdatabase import wcsfitdatabase

log = logging.getLogger(__name__)

def dateformat (starttime=None,endtime=None):
    """ Utility to prettify a plot with dates.
    """

    #plt.xlim([starttime, endtime])
    plt.gcf().autofmt_xdate()
    years = mdates.YearLocator()   # every year
    months = mdates.MonthLocator(bymonth=[4, 7, 10])  # every month
    yearsFmt = mdates.DateFormatter('%Y %b')
    monthformat = mdates.DateFormatter('%b')
    # plt.gca().xaxis.set_major_locator(years)
    # plt.gca().xaxis.set_major_formatter(yearsFmt)
    # plt.gca().xaxis.set_minor_locator(months)
    # plt.gca().xaxis.set_minor_formatter(monthformat)
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

    allfilters = set (allsolutions['filter'])
    dateobs = allsolutions['dateobs']
    log.info ("Set of individual filters: {}".format (allfilters))

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


    for filter in allfilters:

        plt.plot (dateobs, cd11, 'o', label="CD11_{}".format (filter))
        plt.plot (dateobs, cd22, 'o',label="CD22_{}".format (filter))
        plt.legend()
        dateformat()
        plt.savefig ("wcstrend_cddiag_{}.png".format(cameraname))
        plt.close()

    for filter in allfilters:
        plt.plot (dateobs,cd12, 'o',label="CD12_{}".format (filter))
        plt.plot (dateobs,cd21, 'o',label="CD21_{}".format (filter))
        dateformat()
        plt.legend()
        plt.savefig ("wcstrend_cdnondiag_{}.png".format (cameraname))
        plt.close()

    for filter in allfilters:
        plt.plot (dateobs,sip_a11, 'o',label="A_1_1_{}".format (filter))
        plt.plot (dateobs,sip_a20, 'o',label="A_2_0_{}".format (filter))
        plt.plot (dateobs,sip_a02, 'o',label="A_0_2_{}".format (filter))
        plt.legend()
        dateformat()
        plt.savefig ("wcstrend_sipa_{}.png".format(cameraname))
        plt.close()

    for filter in allfilters:
        plt.plot (dateobs,sip_b11, 'o',label="B_1_1_{}".format (filter))
        plt.plot (dateobs,sip_b20, 'o',label="B_2_0_{}".format (filter))
        plt.plot (dateobs,sip_b02, 'o',label="B_0_2_{}".format (filter))
        plt.legend()
        dateformat()
        plt.savefig ("wcstrend_sipb_{}.png".format(cameraname))
        plt.close()


if __name__ == '__main__':
    args = parseCommandLine()

    diagnosticByCamera("ak10", args)
    diagnosticByCamera("fa03", args)
