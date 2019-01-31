""" Tools to store distortion database.
    For LCO, store image file name, DATE-OBS, FILTER, RMS, WCS in json


"""
import sqlite3
import logging
import json

log = logging.getLogger(__name__)



class wcsfitdatabase:
    ''' Storage model for flat field based noise gain measurements'''

    createstatement = "CREATE TABLE IF NOT EXISTS wcsfit (" \
                      "name TEXT PRIMARY KEY, " \
                      "dateobs text," \
                      " camera text," \
                      " filter text," \
                      " rms real," \
                      " nstars int," \
                      " wcsjson text )"

    def __init__(self, fname):
        log.debug("Open data base file %s" % (fname))
        self.sqlite_file = fname
        self.conn = sqlite3.connect(self.sqlite_file)
        self.conn.execute(self.createstatement)
        self.conn.execute("PRAGMA journal_mode=WAL;")
        self.conn.commit()

    def addmeasurement(self, identifier, dateobs, camera, filter, rms, nstars, WCSjson, commit=True):
        with self.conn:
            log.debug("Inserting: %s\n %s %s %s %s %s %s" % (
                identifier, dateobs, camera, filter, rms, nstars, WCSjson))
            self.conn.execute("insert or replace into wcsfit values (?,?,?,?,?,?,?)",
                              (identifier, dateobs, camera, filter, rms, nstars, WCSjson,))

        if (commit):
            self.conn.commit()

    def checkifalreadyused(self, flat1):

        """ Check if a noisegain measurement based on two flat field expsoures already exists or not"""

        with self.conn:
            query = 'select name from noisegain where (name like ?)'
            cursor = self.conn.execute(query, ( "%{}%".format(flat1),))
            allmatch = cursor.fetchall()
            if len(allmatch) > 0:
                log.debug("match found for %s"  % (flat1))
                return True

        log.debug("no match found for %s" % (flat1))
        return False

    def getcameras(self):
        query = "select distinct camera from noisegain"

        cursor = self.conn.execute(query)
        retarray = []

        allcameras = (cursor.fetchall())
        if len(allcameras) == 0:
            log.warning("Zero results returned from query")

        for c in allcameras:
            retarray.append(c[0])
        log.info("Distinct cameras: %s" % retarray)
        return retarray

    def readmeasurements(self, camera=None, filters=None, levelratio=None):
        """

        :param camera: name of the camera to read
        :param filters: array of filters to use. None if no filter selection
        :param levelratio: maximum ratio how much the two flat field levels may vary
        :return: astropy.table with query results. May be None if no results are returned.
        """

        query = "select name,dateobs,camera,filter,extension,gain,readnoise,level,differencenoise,level1,level2 from noisegain " \
                "where (camera like ?) __FILTER__ ORDER BY dateobs"

        queryargs = [camera if camera is not None else '%', ]

        if filters is not None:
            filtercondition = 'AND (filter in (%s))' % ','.join ('?' * len(filters) )
            query = query.replace("__FILTER__",  filtercondition)
            queryargs.extend (filters)
        else:
            query = query.replace ("__FILTER__", "")

        log.debug("Read from database with query\n  %s\n  %s\n" % (query, queryargs))

        cursor = self.conn.execute(query, queryargs)

        allrows = np.asarray(cursor.fetchall())
        if len(allrows) == 0:
            log.warning("Zero results returned from query")
            return None

        t = Table(allrows,
                  names=['identifier', 'dateobs', 'camera', 'filter', 'extension', 'gain', 'readnoise', 'level',
                         'diffnoise', 'level1', 'level2'])
        t['dateobs'] = t['dateobs'].astype(str)
        t['dateobs'] = astt.Time(t['dateobs'], scale='utc', format=None).to_datetime()
        t['gain'] = t['gain'].astype(float)
        t['level'] = t['level'].astype(float)
        t['diffnoise'] = t['diffnoise'].astype(float)
        t['readnoise'] = t['readnoise'].astype(float)
        t['level1'] = t['level1'].astype(float)
        t['level2'] = t['level2'].astype(float)

        if levelratio is not None:
            t = t[abs((t['level1'] - t['level2']) / t['level']) < levelratio]

        return t



    def wcstojson (self, WCS):
        jsonwcs = {}
        jsonwcs['CRVAL1'] = WCS.wcs.crval[0]
        jsonwcs['CRVAL2'] = WCS.wcs.crval[1]
        jsonwcs['CD1_1'] = WCS.wcs.cd[0][0]
        jsonwcs['CD1_2'] = WCS.wcs.cd[0][1]
        jsonwcs['CD2_1'] = WCS.wcs.cd[1][0]
        jsonwcs['CD2_2'] = WCS.wcs.cd[1][1]
        if WCS.sip.a_order > 1:
            jsonwcs['sipa11'] = WCS.sip.a[1][1]
            jsonwcs['sipa20'] = WCS.sip.a[2][0]
            jsonwcs['sipa02'] = WCS.sip.a[0][2]
            jsonwcs['sipb11'] = WCS.sip.b[1][1]
            jsonwcs['sipb20'] = WCS.sip.b[2][0]
            jsonwcs['sipb02'] = WCS.sip.b[0][2]
        return json.dumps(jsonwcs)


    def close(self):
        log.debug("Closing data base file %s " % (self.sqlite_file))
        self.conn.commit()
        self.conn.close()