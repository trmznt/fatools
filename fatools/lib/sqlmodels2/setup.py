

from fatools.lib.utils import cout, cerr
from fatools.lib.params import default_markers, default_panels

def setup(dbh):

    session = dbh.session()

    # create default markers
    for d in default_markers.values():
        marker = dbh.new_marker()
        marker.update(d)
        marker.remark = 'Test __slots__ for %s!' % marker.code
        marker.anything = 'abc'
        cerr("I: marker '%s' created." % marker.code)
        session.add(marker)

    # create default panels
    for d in default_panels.values():
        panel = dbh.new_panel()
        panel.update(d)
        cerr("I: panel '%s' created." % panel.code)
        session.add(panel)

    # create default batch (for bin holder)
    batch = dbh.Batch( code = 'default' )
    batch.fsa_provider = ''
    batch.species = 'X'
    cerr("I: batch 'default' created.")
    session.add(batch)



