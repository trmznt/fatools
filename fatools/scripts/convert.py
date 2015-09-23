

import sys
import argparse
import csv
from collections import defaultdict

from fatools.lib.utils import cout, cerr, cexit, get_dbhandler
from fatools.lib.fautil.traceio import read_abif_stream
from pprint import pprint

def init_argparser( parser=None):

    if parser is None:
        p = argparse.ArgumentParser( 'convert' )
    else:
        p = parser

    ## commands

    p.add_argument('--fsa2tab', default=False, action='store_true',
        help = 'convert from FSA to tab file')
    p.add_argument('--genemapper2tab', default=False, action='store_true',
        help = 'convert genemapper CSV file to fatools assay info tab file')
    p.add_argument('--checkfsa', default=False, action='store_true',
        help = 'check FSA files')

    ## options

    p.add_argument('--sqldb', default=False,
        help = 'SQLITE3 database filename')
    p.add_argument('--fsdb', default=False,
        help = 'root directory for filesystem-based database')
    p.add_argument('--species', default=False,
        help = 'species for markers')
    p.add_argument('--fsadir', default=False,
        help = 'root directory for FSA files')

    ## mandatory options

    p.add_argument('infiles', nargs='+')

    return p


def main(args):

    do_convert(args)


def do_convert(args, dbh=None):

    if not dbh and (args.sqldb or args.fsdb):
        dbh = get_dbhandler(args)

    if args.fsa2tab:
        do_fsa2tab(args)
    elif args.genemapper2tab:
        do_genemapper2tab(args, dbh)
    elif args.checkfsa:
        do_checkfsa(args)
    else:
        cerr('Unknown command, nothing to do!')
        return False
    return True


def do_fsa2tab( args ):

    for infile in args.infiles:
        with open( infile, 'rb' ) as instream:
            t = read_abif_stream( instream )
        channels = t.get_channels()
        names = [ '"' + c + '"' for c in channels ]
        print("Dyes: %s" % ' '.join( channels ))
        with open( infile + '.raw.tab', 'wt') as out:
            out.write( '\t'.join( names ) )
            out.write( '\n' )
            for p in zip( * [ channels[c].raw for c in channels ] ):
                out.write( '\t'.join( str(x) for x in p ) )
                out.write( '\n' )
        with open( infile + '.base.tab', 'wt') as out:
            out.write( '\t'.join( names ) )
            out.write( '\n' )
            for p in zip( * [ channels[c].smooth() for c in channels ] ):
                out.write( '\t'.join( str(x) for x in p ) )
                out.write( '\n' )


def do_genemapper2tab(args, dbh):

    species = None
    if args.species: species = args.species

    for infile in args.infiles:

        sample_set = defaultdict(list)
        csv_in = csv.DictReader( open(infile) )
        assay_list = {}

        for row in csv_in:
            assay = row['Sample File']
            sample = row['Sample Name']
            run_name = row['Run Name']
            panel = row['Panel']
            marker = row['Marker']

            if assay in assay_list:
                if assay_list[assay] != run_name:
                    cexit('Inconsistence or duplicate FSA file name: %s' % assay)
            else:
                assay_list[assay] = run_name

            token = (sample, assay, panel)
            sample_set[token].append( marker )


        outfile = open(infile + '.tab', 'w')
        outfile.write('SAMPLE\tASSAY\tPANEL\tOPTIONS\n')

        for token in sorted(sample_set.keys()):
            sample, assay, panel = token
            markers = sample_set[token]

            db_panel = dbh.get_panel(panel)
            s_panel_markers = set( x.upper() for x in db_panel.get_marker_codes())
            s_assay_markers = set(
                    ('%s/%s' % (species, x) if (species and '/' not in x) else x).upper()
                    for x in markers )

            excludes = s_panel_markers - s_assay_markers
            if s_assay_markers - s_panel_markers:
                cexit('ERROR inconsistent marker(s) for sample %s assay %s: %s' % 
                    (sample, assay, str(s_assay_markers - s_panel_markers)))

            if excludes:
                excludes = 'exclude=%s' % ','.join(excludes)
            else:
                excludes = ''

            outfile.write('%s\t%s\t%s\t%s\n' % (sample, assay, panel, excludes))

        outfile.close()


def do_checkfsa(args):

    fsadir = args.fsadir or '.'

    for infile in args.infiles:
        data = csv.DictReader( open(infile), delimiter='\t')

        files = {}
        line = 2
        for row in data:
            sample = row['SAMPLE']
            if sample.startswith('#'):
                line += 1
                continue
            assay_file = row['ASSAY']
            panel = row['PANEL']
            if assay_file in files:
                cerr('WARN file: %s - duplicated assay: %s for sample %s panel %s' % 
                        (infile, assay_file, sample, panel))
            files[assay_file] = True
            try:
                with open( '%s/%s' % (fsadir , assay_file), 'rb' ) as instream:
                    t = read_abif_stream( instream )
                line += 1
            except:
                cerr('ERR file: %s line: %d  - sample: %s assay: %s' %
                        (infile, line, sample, assay_file))
                #raise
                line += 1

