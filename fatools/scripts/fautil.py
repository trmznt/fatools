
import sys, argparse, yaml, os

from fatools.lib.utils import cout, cerr, cexit


def init_argparser():

    p = argparse.ArgumentParser('fautil')

    p.add_argument('--info', default=False, action='store_true',
        help = 'get information on FA assay')

    p.add_argument('--view', default=False, action='store_true',
        help = 'view information')

    p.add_argument('--analyze', default=False, action='store_true',
        help = 'analyze single FSA file')

    p.add_argument('--file', default=False,
        help = 'input file')

    p.add_argument('--sqldb', default=False,
        help = 'Sqlite database file')

    p.add_argument('--sizestandard', default='LIZ600',
        help = 'Size standard')

    return p


cache_traces = {}


def main(args):

    do_fautil(args)



def do_fautil(args):

    if args.sqldb:
        dbh = get_dbhandler(args)
    else:
        dbh = None

    if args.info is not False:
        do_info(args, dbh)
    if args.view is not False:
        do_view(args, dbh)
    if args.analyze is not False:
        do_analyze(args)



def get_traces(args, dbh):

    traces = []

    if dbh is None:
        # get from infile
        infile = args.file
        if infile is False:
            cexit('E - Please provide a filename or Sqlite database path')

        abspath = os.path.abspath( args.file )

        if abspath in cache_traces:
            traces.append((abspath, cache_traces[abspath]))

        else:
            from fatools.lib.fautil.traceio import read_abif_stream
            with open( abspath, 'rb') as instream:
                t = read_abif_stream(instream)
                cache_traces[abspath] = t
                traces.append((abspath, t))

    else:
        pass

    return traces



def do_info(args, dbh):


    traces = get_traces(args, dbh)

    for abspath, trace in traces:
        cout('I - trace: %s' % abspath)
        cout('I - runtime: %s' % trace.get_run_start_time())



def do_view(args, dbh):

    traces = get_traces(args, dbh)

    from fatools.lib.gui.viewer import viewer

    for abspath, trace in traces:
        
        viewer( trace )


def do_analyze(args):
    """ open a tracefile, performs fragment analysis (scan & call only)
    """

    from fatools.lib.fautil.traceio import read_abif_stream
    from fatools.lib.fautil.traceutils import separate_channels
    from fatools.lib.fsmodels.models import Assay, Marker, Panel
    from fatools.lib import params

    scanning_parameter = params.Params()

    # create dummy markers
    ladder = Marker('ladder', 10, 600, 0, None)

    # create dummy panel
    dummy_panel = Panel( '-', {
        'ladder': args.sizestandard,
        'markers': {},
    })

    with open(args.file, 'rb') as in_stream:
        cerr('Reading FSA file: %s' % args.file)
        t = read_abif_stream(in_stream)

    # create a new Assay and add trace
    assay = Assay()
    assay.size_standard = args.sizestandard
    assay._trace = t

    # create all channels
    assay.create_channels()

    # assign all channels
    assay.assign_channels( panel = dummy_panel )


    # scan for peaks
    assay.scan(scanning_parameter)


    # scan all channels


