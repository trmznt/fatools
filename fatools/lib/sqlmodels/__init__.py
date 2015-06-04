

## handler

from fatools.lib.sqlmodels import schema
from fatools.lib.utils import cout, cerr
import os, sys

class SQLHandler(object):


    Marker = schema.Marker

    Batch = schema.Batch

    def __init__(self, dbfile, initial=False):
        print("Opening db: %s" % dbfile)
        if not initial and not os.path.isfile(dbfile):
            cerr('ERR - sqlite db file not found: %s' % dbfile)
            sys.exit(1)
        if initial and os.path.isfile(dbfile):
            cerr('ERR - sqlite db file already exists: %s' % dbfile)
            sys.exit(1)
        self.dbfile = dbfile
        self.engine, self.session = schema.engine_from_file(dbfile)


    def initdb(self):
        schema.Base.metadata.create_all(self.engine)
        from fatools.lib.sqlmodels.setup import setup
        setup( self.session )
        cout('Database at %s has been initialized.' % self.dbfile)


    def get_batch(self, batch_code = None):

        if not batch_code:
            cerr('ERR - batch code must be supplied!')
            sys.exit(1)

        return self.Batch.search(batch_code, self.session)


    def get_marker(self, marker_code = None):

        if not marker_code:
            cerr('ERR - marker code must be supplied')
            sys.exit(1)

        return self.Marker.search(marker_code, self.session)


    def Panel(self):
        p = schema.Panel()
        p._dbh_session_ = self.session
        return p

