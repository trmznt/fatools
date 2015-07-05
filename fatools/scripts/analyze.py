#
# this script contains commands which do not need to update the database content
# commands that need to update/change the database content should be in facmd toolset
#

import sys, argparse, yaml

from fatools.lib.analytics.query import Query, load_yaml
from fatools.lib.utils import cout, cerr, cexit, get_dbhandler
from fatools.lib import params
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

    p.add_argument('--allelesummary', default=False, action='store_true',
            help = 'report allele summary')

    p.add_argument('--binsummary', default=False, action='store_true',
            help = 'report bin summary')

    p.add_argument('--adjustbins', default=False, action='store_true',
            help = 'iteratively adjust the binning')

    p.add_argument('--export', default=False, action='store_true',
            help = 'export allele data to file')

    ## Options

    p.add_argument('--yamlquery', default=False,
            help = 'YAML query file')

    p.add_argument('--outformat', default=False,
            help = 'format output type (html, tab, arlequin)')

    p.add_argument('--outfile', default=False,
            help = 'output filename, or - for stdout/console')

    p.add_argument('--outplot', default=False,
            help = 'output plot filename')

    p.add_argument('--iteration', default=2, type=int,
            help = 'iteration number')

    ## Override params

    p.add_argument('--sample_qual_threshold', default=-1, type=float,
            help = 'sample quality threshold')

    p.add_argument('--rel_threshold', default=-1, type=float,
            help = 'relative allele rfu threshold')

    p.add_argument('-m', '--markers', default='',
            help = 'markers')


    return p


def main(args):

    do_analyze(args)


def do_analyze(args, dbhandler_func = get_dbhandler):

    dbh = dbhandler_func( args )


    if args.samplesummary:
        do_samplesummary(args, dbh)
    elif args.allelesummary:
        do_allelesummary(args, dbh)
    elif args.binsummary:
        do_binsummary(args, dbh)
    elif args.export:
        do_export(args, dbh)



def do_samplesummary(args, dbh):

    query = get_query( args, dbh )
    sample_sets = query.get_filtered_sample_sets()
    cout( make_sample_report(sample_sets) )



def do_allelesummary(args, dbh):

    from fatools.lib.analytics.summary import summarize_alleles, plot_alleles

    query = get_query( args, dbh )
    analytical_sets = query.get_filtered_analytical_sets()
    report = summarize_alleles( analytical_sets )
    cout( make_sample_report( analytical_sets.get_sample_sets() ) )
    cout( make_allele_report(report) )

    if args.outplot:
        plot_alleles( report, args.outplot )


def do_binsummary(args, dbh):

    from fatools.lib.analytics.summary import summarize_bins

    scanning_parameter = params.Params()

    markers = None
    for i in range( args.iteration ):
        query = get_query( args, dbh )
        analytical_sets = query.get_filtered_analytical_sets()
        report = summarize_bins( analytical_sets )
        cerr('I: Bin summary iteration %d' % i)
        pprint(report)

        markers = []
        for (marker_id, updated_bins) in report.items():
            marker = dbh.get_marker_by_id(marker_id)
            marker.adjustbins( updated_bins )
            markers.append( marker )
        dbh.session().flush()

        # rebinning
        cerr('I: Rebinning samples')
        assay_list = []
        N = len(analytical_sets.sample_ids)
        count = 1
        for sample_id in analytical_sets.sample_ids:
            sample = dbh.get_sample_by_id(sample_id)
            cerr('\rI: [%d/%d] - Binning sample...' % (count, N), nl=False)
            for assay in sample.assays:
                assay.bin( scanning_parameter.nonladder, markers )
            count += 1
        cerr('')
        dbh.session().flush()

    if args.outfile:
        
        output_dict = {}
        for marker in markers:
            output_dict[marker.label] = {
                    'label': marker.label,
                    'bins': marker.bins
                }

        with open(args.outfile, 'wt') as f:
            yaml.dump(output_dict, f)
        cerr('I: writing bins to %s' % args.outfile)



def do_export(args, dbh):

    from fatools.lib.analytics.export import export

    query = get_query( args, dbh )
    analytical_sets = query.get_filtered_analytical_sets()
    if analytical_sets.total_samples <= 0:
        cexit('ERR - query does not yield any sample data')
    else:
        cerr('INFO - total sampel number: %d' % analytical_sets.total_samples)
    output = export( analytical_sets, dbh, outfile = args.outfile, format = args.outformat )
    cout('Done.')


def do_corralleles(args, dbh):

    from fatools.lib.analysis.correlation import correlate_alleles

    query = get_query( args, dbh )
    analytical_sets = query.get_filtered_analytical_sets()
    for marker_code in analytical_sets.marker_ids:

        report = correlate_alleles(analytical_sets[0], analytical_sets[1], marker=marker_code)
        cout( make_correlate_report( report ) )



def get_sample_sets(args, dbh):
    
    pass


def get_analytical_sets( args, dbh ):

    
    query = load_yaml( open(args.yamlquery).read() )
    sample_sets = query['selector'].get_sample_sets(dbh)
    pprint(sample_sets)


def get_query( args, dbh ):

    query_params = load_yaml( open(args.yamlquery).read() )
    if args.sample_qual_threshold >= 0:
        query_params['filter'].sample_qual_threshold = args.sample_qual_threshold
    if args.markers:
        query_params['filter'].markers = args.markers
    if args.rel_threshold >= 0:
        query_params['filter'].rel_threshold = args.rel_threshold
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


def make_allele_report( summaries ):

    #sample_sets = analytical_sets.get_sample_sets()

    #sample_report = make_sample_report( sample_sets )

    lines = []; _ = lines.append
    _('ALLELE SUMMARY')
    _('==================================')

    for label in summaries:
        summary = summaries[label]['summary']
        _('Sample Set: %s' % label)

        for marker_id in summary:
            _('    Marker ID: %d' % marker_id)
            _('    Unique alleles: %d' % summary[marker_id]['unique_allele'])
            _('    Total alleles: %d' % summary[marker_id]['total_allele'])

            for data in summary[marker_id]['alleles']:
                _('        %3d  %5.3f  %3d  %5.2f - %5.2f  %5.2f  %4.2f' %
                        (data[0], data[1], data[2], data[4], data[5], data[8], data[6]))
    
    return '\n'.join(lines)


    
