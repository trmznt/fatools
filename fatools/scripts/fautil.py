
import sys, argparse, yaml, os

from fatools.lib.utils import cout, cerr, cexit


def init_argparser():

    p = argparse.ArgumentParser('fautil')

    p.add_argument('--info', default=False, action='store_true',
        help = 'get information on FA assay')

    p.add_argument('--view', default=False, action='store_true',
        help = 'view information')

    p.add_argument('--file', default=False,
        help = 'input file')

    p.add_argument('--sqldb', default=False,
        help = 'Sqlite database file')

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
