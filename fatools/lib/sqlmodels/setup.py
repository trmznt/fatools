
from fatools.lib.sqlmodels import schema
from fatools.lib.utils import cout, cerr


def setup(session):

    # create undefined marker
    marker = schema.Marker( code = 'undefined' )
    cout("INFO - marker 'undefined' created.")
    session.add(marker)

    # create ladder marker
    marker = schema.Marker( code =  'ladder' )
    cout("INFO - marker 'ladder' created.")
    session.add(marker)

    # create combined marker
    marker = schema.Marker( code = 'combined' )
    cout("INFO - marker 'combined' created.")
    session.add(marker)

    # create default panel
    panel = schema.Panel( code = 'undefined' )
    cout("INFO - panel 'undefined' created.")
    session.add(panel)


