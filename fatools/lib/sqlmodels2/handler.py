
import sys, os
from fatools.lib.utils import cout, cerr
from fatools.lib.sqlmodels.handler_interface import base_sqlhandler
from fatools.lib.sqlmodels2 import schema


class SQLHandler(base_sqlhandler):

    Panel = schema.Panel
    Marker = schema.Marker
    Batch = schema.Batch
    Sample = schema.Sample
    FSA = schema.FSA
    Channel = schema.Channel
    AlleleSet = schema.AlleleSet
    Allele = schema.Allele


    def __init__(self, dbfile, initial=False):
        cerr("Opening db: %s" % dbfile)
        if not initial and not os.path.isfile(dbfile):
            cerr('ERR - sqlite db file not found: %s' % dbfile)
            sys.exit(1)
        if initial and os.path.isfile(dbfile):
            cerr('ERR - sqlite db file already exists: %s' % dbfile)
            sys.exit(1)
        self.dbfile = dbfile
        self.engine, self.session = schema.engine_from_file(dbfile)


    def initdb(self, create_table = True):
        if create_table:
            schema.Base.metadata.create_all(self.engine)
        from fatools.lib.sqlmodels2.setup import setup
        setup( self )
        cerr('Database at %s has been initialized.' % self.dbfile)

