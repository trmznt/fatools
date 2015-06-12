

import sys, argparse

from fatools.lib.analytics.query import Query, load_yaml
from fatools.lib.utils import cout, cerr, cexit, get_dbhandler
from pprint import pprint

def init_argparser(parser=None):

    if parser:
        p = parser
    else:
        p = argparse.ArgumentParser('analyze')

    p.add_argument('--sqldb', default=False,
            help = 'SQLite3 database filename')

    p.add_argument('--fsdb', default=False,
            help = 'directory for filesystem-based database')

    ## Commands

    p.add_argument('--samplesummary', default=False, action='store_true',
            help = 'report sample summary')


    ## Options

    p.add_argument('--yamlquery', default=False,
            help = 'YAML query file')

    p.add_argument('--outformat', default=False,
            help = 'format output type (html, tab)')


    return p


def main(args):

    do_analyze(args)


def do_analyze(args, dbhandler_func = get_dbhandler):

    dbh = dbhandler_func( args )


    if args.samplesummary:
        do_samplesummary(args, dbh)



def do_samplesummary(args, dbh):

    query = get_query( args, dbh )
    sample_sets = query.get_filtered_sample_sets()
    cout( make_sample_report(sample_sets) )


def get_sample_sets(args, dbh):
    
    pass


def get_analytical_sets( args, dbh ):

    
    query = load_yaml( open(args.yamlquery).read() )
    sample_sets = query['selector'].get_sample_sets(dbh)
    pprint(sample_sets)


def get_query( args, dbh ):

    query_params = load_yaml( open(args.yamlquery).read() )
    return Query( query_params, dbh )


def make_sample_report( sample_sets ):
    lines = []
    _ = lines.append
    _( 'SAMPLE SUMMARY' )
    _( '================================' )
    _( 'Group: %d' % len(sample_sets) )
    _( 'Total samples: %d' % sample_sets.total_samples )
    _( '--------------------------------------------' )
    for s in sample_sets:
        _( '  %-20s  %4d' % (s.label, s.N) )
    _( '--------------------------------------------' )

    return '\n'.join(lines)

    
