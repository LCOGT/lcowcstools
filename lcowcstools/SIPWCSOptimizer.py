import argparse
import copy
import logging
import os
import numpy as np
from astropy.wcs import Sip
from scipy import optimize as optimize

from lcowcstools.CatalogMatcher import CatalogMatcher
from lcowcstools.ReferenceCatalogProvider import refcat2, online_refcat2
from lcowcstools.wcsfitsdatabase import wcsfitdatabase

logging.getLogger('matplotlib.font_manager').disabled = True
logging.getLogger('PIL.PngImagePlugin').disabled = True

__author__ = 'drharbeck@gmail.com'

log = logging.getLogger(__name__)


class SIPOptimizer:
    """ Given a matched catalog, iteratively fit a WCS with SIP distortions up to second order.

    """

    def __init__(self, newMatchedCatalog, maxorder=2):

        self.matchedCatalog = newMatchedCatalog
        if self.matchedCatalog.wcs is None:
            log.error("Cannot proceed without a wcs in matched catalog. Aborting.")
            return None

        # bootstrap the initial SIP wcs
        self.maxorder = maxorder
        crval = self.matchedCatalog.wcs.wcs.crval
        cd = self.matchedCatalog.wcs.wcs.cd

        # 6 paramters
        self.wcsfitparams = [crval[0], crval[1], cd[0][0], cd[0][1], cd[1][0], cd[1][1]]

        # for the second order, we need 6 paramters, x^2, xy, and  y^2 for each the x and y axis
        # total of 12 paramters
        if maxorder > 1:
            self.wcsfitparams.extend([0, 0, 0, 0, 0, 0])

        # for the third order, we need to add even more parameters. This is work in progress.
        # total of  20 paramters
        if maxorder > 2:
            self.wcsfitparams.extend([0, 0, 0, 0, 0, 0, 0, 0])
        # total of 30
        if maxorder > 3:
            self.wcsfitparams.extend([0, 0, 0, 0,0,0,0,0,0,0])

        log.info (f"Now optimizing to {maxorder} order. I have a total of {len(self.wcsfitparams)} parameters.")
        # We calculate the inital merrit function before anything else to get a baseline.
        merrit = SIPOptimizer.merritFunction(self.wcsfitparams, self.matchedCatalog)
        # log.info("SIPOptimizer init: merrit function is %12.7f" % (merrit))

    @staticmethod
    def merritFunction(sipcoefficients, matchedCatalog):
        """ Calculate a merrit function for a matched catalog based on current understadning of WCS.
        This function is at the core of the optimizer, as it is the merrit function for the fitting function.
        """

        # update the matched Catalog's WCS
        # 1ST ORDER LINEAR FIT
        matchedCatalog.wcs.wcs.crval[0] = sipcoefficients[0]
        matchedCatalog.wcs.wcs.crval[1] = sipcoefficients[1]
        matchedCatalog.wcs.wcs.cd[0][0] = sipcoefficients[2]
        matchedCatalog.wcs.wcs.cd[0][1] = sipcoefficients[3]
        matchedCatalog.wcs.wcs.cd[1][0] = sipcoefficients[4]
        matchedCatalog.wcs.wcs.cd[1][1] = sipcoefficients[5]

        m = 4
        sip_a = np.zeros((m + 1, m + 1), np.double)
        sip_b = np.zeros((m + 1, m + 1), np.double)

        if len(sipcoefficients) > 6:
        # 2nd order fit
            sip_a[1][1] = sipcoefficients[6]
            sip_a[2][0] = sipcoefficients[7]
            sip_a[0][2] = sipcoefficients[8]
            sip_b[1][1] = sipcoefficients[9]
            sip_b[2][0] = sipcoefficients[10]
            sip_b[0][2] = sipcoefficients[11]

        if len(sipcoefficients) > 12:
            #3rd order fit
            sip_a[3][0] = sipcoefficients[12]
            sip_a[0][3] = sipcoefficients[13]
            sip_b[3][0] = sipcoefficients[14]
            sip_b[0][3] = sipcoefficients[15]

            sip_a[2][1] = sipcoefficients[16]
            sip_a[1][2] = sipcoefficients[17]
            sip_b[2][1] = sipcoefficients[18]
            sip_b[1][2] = sipcoefficients[19]



        if len (sipcoefficients) > 20:
            # 4th order fit
            sip_a[0][4] = sipcoefficients[20]
            sip_a[1][3] = sipcoefficients[21]
            sip_a[2][2] = sipcoefficients[22]
            sip_a[3][2] = sipcoefficients[23]
            sip_a[4][0] = sipcoefficients[24]

            sip_b[0][4] = sipcoefficients[25]
            sip_b[1][3] = sipcoefficients[26]
            sip_b[2][2] = sipcoefficients[27]
            sip_b[3][2] = sipcoefficients[28]
            sip_b[4][0] = sipcoefficients[29]


        sip = Sip(sip_a, sip_b,
                  None, None, matchedCatalog.wcs.wcs.crpix)
        matchedCatalog.wcs.sip = sip

        matchedCatalog.wcs.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        matchedCatalog.wcs.wcs.set()

        # Get the rms of the matched catalog.
        merrit = matchedCatalog.updateWCSandUpdateRMS(matchedCatalog.wcs)
        return merrit

    def improveSIP(self):

        bestfit = optimize.minimize(SIPOptimizer.merritFunction, self.wcsfitparams, args=(self.matchedCatalog), method = 'Nelder-Mead')
        merrit = SIPOptimizer.merritFunction(bestfit.x, self.matchedCatalog)
        log.debug("Optimizer return        {: 10.4f}".format(merrit))
        log.debug (f"Fit result:\n{bestfit}")


def iterativelyFitWCSmany(images, args, refcat=None):
    """ Wrapper to optimize the WCS for a set of images
    """

    if refcat is None:
        refcat =  online_refcat2(args.refcat2)

    if len(images) > 0:
        for image in images:
            iterativelyFitWCSsingle(image, args, searchradii=args.searchradii, refcat=refcat)


def iterativelyFitWCSsingle(image, args, searchradii, refcat=None):
    """ Logistic to start optimizing WCS for a single image.
    """

    log.info("Starting to process {}".format(image))

    if refcat is None:
        refcat = online_refcat2(args.refcat2)

    pngbasename = os.path.basename(image)

    if args.database:
        wcsdb = wcsfitdatabase(args.database)
        if wcsdb.checkifalreadyused(pngbasename) and not args.reprocess:
            log.info("File already measured. Not repeating work; skipping")
            wcsdb.close()
            return

    matchedCatalog = CatalogMatcher.createMatchedCatalogForLCO(
        image,
        refcat, searchradii[0], minobjects=args.minmatched, undistort=args.undistort)

    if matchedCatalog is None:
        log.info("returned empty catalog, not continuing.")
        return

    if args.makepng:
        matchedCatalog.diagnosticPlots('{:s}_prefit'.format(pngbasename))

    # Preserve the initial pointing
    initialPointing = copy.deepcopy(matchedCatalog.wcs.wcs.crval)
    # do a full fit
    if (matchedCatalog.matchedCatalog is None) or (len(matchedCatalog.matchedCatalog['x']) < args.minmatched):
        log.warning("Not enough stars in input catalog: %s found, %d are required to start. Giving up" % (
            'None' if matchedCatalog.matchedCatalog is None else len(matchedCatalog.matchedCatalog['x']), args.minmatched))
        return

    opt = SIPOptimizer(matchedCatalog, maxorder=args.fitorder)
    opt.improveSIP()

    for searchradius in searchradii[1:]:
        matchedCatalog.matchCatalogs(matchradius=searchradius)
        opt = SIPOptimizer(matchedCatalog, maxorder=args.fitorder)
        opt.improveSIP()

    if args.makepng:
        matchedCatalog.diagnosticPlots('{:s}_postfits'.format (pngbasename))

    if (args.database):
        wcsdb.addmeasurement(pngbasename, matchedCatalog.dateobs, matchedCatalog.camera, matchedCatalog.filter, None,
                             None, matchedCatalog.azimuth, matchedCatalog.altitude,
                             wcsdb.wcstojson(matchedCatalog.wcs))
    log.info(wcsdb.wcstojson(matchedCatalog.wcs))

    # Final report
    finalPointing = matchedCatalog.wcs.wcs.crval
    log.info('Fitting updated pointing at CRPIX by {}"'.format((finalPointing - initialPointing) * 3600.))
    if (args.database):
        wcsdb.close()


def parseCommandLine():
    parser = argparse.ArgumentParser(
        description='LCO WCS Tool')

    parser.add_argument('--inputfiles', type=str, nargs='+', help="FITS file for which to derive the WCS function.")
    parser.add_argument('--refcat2', type=str, default='http://phot-catalog.lco.gtn/',
                        help='Location of Atlas refcat2 catalog in slite forrmat')
    parser.add_argument('--minmatched', type=int, default=50,
                        help='Minimum number of matched stars to accept solution or even proceed to fit.')
    parser.add_argument('--fitorder', type=int, default=2)
    parser.add_argument('--makepng', action='store_true', help="Create a png output of wcs before and after fit.")
    parser.add_argument('--undistort', action='store_true',
                        help="Undistort input coordinate catalog before WCS fitting.")
    parser.add_argument('--reprocess', action='store_true',
                        help="Reprocess even though file may have been processed already.")
    parser.add_argument('--loglevel', dest='log_level', default='INFO', choices=['DEBUG', 'INFO', 'WARN'],
                        help='Set the debug level')
    parser.add_argument('--searchradii', type=float, nargs='+', default=[10, 10, 5, 3, 2])
    parser.add_argument ("--catradius", default = 0.2, help="Search radius in the reference catalog")
    parser.add_argument('--database', default="wcsfits.sqlite")
    args = parser.parse_args()

    logging.basicConfig(level=getattr(logging, args.log_level.upper()),
                        format='%(asctime)s.%(msecs).03d %(levelname)7s: %(module)20s: %(message)s')
    return args


def main():
    args = parseCommandLine()
    log.info("Processing %d input files" % len(args.inputfiles))
    iterativelyFitWCSmany(args.inputfiles, args)


if __name__ == '__main__':
    main()


