# provide commands for Fragment Analysis (FA)

from fatools.lib import params
from fatools.lib.utils import cerr, cout, cverr, cexit, tokenize, detect_buffer, set_verbosity

import argparse, yaml, csv, os
from io import StringIO


def init_argparser(parser=None):

    p = parser if parser else argparse.ArgumentParser('facmd')

    p.add_argument('--sqldb', default=False,
            help = 'SQLite3 database filename')

    p.add_argument('--file', default=None,
            help = "Comma-separated FSA filenames (optional)")

    p.add_argument('--fsdb', default=None,
            help = 'Filesystem-based database')

    p.add_argument('--infile', default=None,
            help = 'Tab-delimited or CSV manifest file')


    # command in ascending order

    p.add_argument('--clear', default=False, action='store_true',
            help = 'clear (or remove) all peaks from FSA file')

    p.add_argument('--align', default=False, action='store_true',
            help = 'scan ladder channel, preannotate and align with size standards')

    p.add_argument('--call', default=False, action='store_true',
            help = 'scan non-ladder channels, preannotate peaks and determine their sizes')

    p.add_argument('--bin', default=False, action='store_true',
            help = 'bin non-ladder peaks')

    p.add_argument('--annotate', default=False, action='store_true',
            help = 'annotate non-ladder peaks')

    p.add_argument('--plot', default=False, action='store_true',
            help = 'plot normalized trace')

    p.add_argument('--dendogram', default=False, action='store_true',
            help = 'plot dendograms of ladders and alleles')

    # semi-mandatory

    p.add_argument('--panel', default="",
            help = 'comma-separated panel code(s)')

    p.add_argument('--marker', default="",
            help = 'comma-separated marker code(s)')

    p.add_argument('--panelfile', default="",
            help = 'YAML panel file')

    p.add_argument('--markerfile', default="",
            help = "YAML marker file")

    # options

    p.add_argument('--cluster', default=0, type=int,
            help = 'number of cluster for hierarchical clustering alignment')

    p.add_argument('--verbose', default=0, type=int,
            help = 'show verbosity')

    p.add_argument('--use-cache', default=False, action='store_true',
            help = 'prepare to use caches')

    p.add_argument('--no-cache', default=False, action='store_true',
            help = 'do not use caches')

    p.add_argument('--commit', default=False, action='store_true',
            help = 'commit to database')

    return p


def main(args):

    if args.verbose != 0:
        set_verbosity(args.verbose)

    dbh = None

    if args.file or args.infile:
        cverr(4, 'D: opening FSA file(s)')
        fsa_list = open_fsa(args)
    elif dbh is None:
        cverr(4, 'D: connecting to database')
        dbh = get_dbhandler(args)
        fsa_list = get_fsa_list(args, dbh)

    cerr('I: obtained %d FSA' % len(fsa_list))

    if args.commit:
        with transaction.manager:
            do_facmd(args, fsa_list, dbh)
            cerr('** COMMIT to database **')
    elif dbh:
        cerr('WARNING ** running without database COMMIT! All changes will be discarded!')
        if not ( args.test or args.y ):
            keys = input('Do you want to continue [y/n]? ')
            if not keys.lower().strip().startswith('y'):
                sys.exit(1)
        do_facmds(args, fsa_list, dbh)
    else:
        do_facmds(args, fsa_list)


def do_facmds(args, fsa_list, dbh=None):

    executed = 0

    if args.clear:
        do_clear( args, fsa_list, dbh )
        executed += 1
    if args.align:
        do_align( args, fsa_list, dbh )
        executed += 1
    if args.call:
        do_call( args, fsa_list, dbh )
        executed += 1
    if args.plot:
        do_plot( args, fsa_list, dbh )
        executed += 1
    if args.dendogram:
        do_dendogram( args, fsa_list, dbh)
        executed += 1

    if executed == 0:
        cerr('W: please provide a relevant command')
    else:
        cerr('I: executed %d command(s)' % executed)


def do_clear( args, fsa_list, dbh ):
    pass


def do_align( args, fsa_list, dbh ):

    cerr('I: Aligning size standards...')

    for (fsa, sample_code) in fsa_list:
        cverr(3, 'D: aligning FSA %s' % fsa.filename)
        fsa.align(params.Params())


def do_call( args, fsa_list, dbh ):

    cerr('I: Calling non-ladder peaks...')

    for (fsa, sample_code) in fsa_list:
        cverr(3, 'D: calling FSA %s' % fsa.filename)
        fsa.call(args.marker)


def do_plot( args, fsa_list, dbh ):

    cerr('I: Creating plot...')

    from matplotlib import pylab as plt

    for (fsa, sample_code) in fsa_list:
        for c in fsa.channels:
            plt.plot(c.data)

        plt.show()


def do_dendogram( args, fsa_list, dbh ):

    from fatools.lib.fautil import hclustalign
    from matplotlib import pyplot as plt

    for (fsa, sample_code) in fsa_list:

        c = fsa.get_ladder_channel()
        c.scan(params.Params()) # scan first if necessary

        ladder = fsa.panel.get_ladder()
        peaks = c.get_alleles()

        #initial_pair, P, L = hclustalign.hclust_align(peaks, ladder)
        P = hclustalign.generate_tree( [ (n.rtime, 0) for n in peaks ] )
        L = hclustalign.generate_tree( [ (e, 0) for e in ladder['sizes'] ] )

        clusters = hclustalign.fcluster(L.z, args.cluster or ladder['k'], criterion="maxclust")
        print(clusters)

        clusters = hclustalign.fcluster(P.z, args.cluster or ladder['k'], criterion="maxclust")
        print(clusters)

        plt.figure()
        plt.subplot(121)
        hclustalign.dendrogram(L.z, leaf_rotation=90, leaf_font_size=8,
                labels = [ x[0] for x in L.p ])
        plt.subplot(122)
        hclustalign.dendrogram(P.z, leaf_rotation=90, leaf_font_size=8,
                labels = [ x[0] for x in P.p ])
        plt.show()


def open_fsa( args ):
    """ open FSA file(s) and prepare fsa instances
        requires: args.file, args.panel, args.panelfile
    """

    from fatools.lib.fileio.models import Marker, Panel, FSA

    if not args.panel:
        cexit('ERR: using FSA file(s) requires --panel argument!')

    if not args.panelfile:
        cerr('WARN: using default built-in panels')
        Panel.upload(params.default_panels)
    else:
        with open(args.panelfile) as f:
            # open a YAML file that describe panel sets
            Panel.upload(yaml.load(f))

    if not args.markerfile:
        Marker.upload(params.default_markers)
    else:
        raise NotImplementedError()

    panel = Panel.get_panel(args.panel)
    fsa_list = []
    index = 1

    # prepare caching
    if args.use_cache:
        if not os.path.exists('.fatools_caches/channels'):
            os.makedirs('.fatools_caches/channels')

    if args.file:
        for fsa_filename in args.file.split(','):
            fsa_filename = fsa_filename.strip()
            fsa = FSA.from_file(fsa_filename, panel, cache = not args.no_cache)
            # yield (fsa, str(i))
            fsa_list.append( (fsa, str(index)) )
            index += 1

    elif args.infile:

        with open(args.infile) as f:
            buf, delim = detect_buffer( f.read() )
        inrows = csv.DictReader( StringIO(buf), delimiter=delim )
        line = 1
        index = 1

        for r in inrows:

            line += 1

            fsa_filename = r['FILENAME'].strip()
            if fsa_filename.startswith('#'):
                continue

            if r.get('OPTIONS', None):
                options = tokenize( r['OPTIONS'] )
            else:
                options = None

            panel_code = r.get('PANEL', None) or args.panel
            panel = Panel.get_panel(panel_code)

            fsa = FSA.from_file( fsa_filename, panel, options, cache = not args.no_cache )
            if 'SAMPLE' in inrows.fieldnames:

                # yield (fsa, r['SAMPLE'])
                fsa_list.append( (fsa, r['SAMPLE']) )
            else:

                # yield (fsa, str(index))
                fsa_list.append( (fsa, str(index)) )
                index += 1

    return fsa_list


def get_fsa_list( args, dbh ):
    """
    get fsa instance from database based on parameters in args
    """

    if not args.batch:
        cexit('ERR: using database requires --batch argument!', 1)

    batch = dbh.get_batch( args.batch )
    if not batch:
        cexit('ERR: batch %s not found!' % args.batch, 1)

    samples = []
    if args.sample:
        samples = args.sample.split(',')

    fsas = []
    if args.fsa:
        fsas = args.assay.split(',')

    panels = []
    if args.panel:
        panels = args.panel.split(',')

    markers = []
    if args.marker:
        markers = dbh.get_markers(args.panel.split(','))

    fsa_list = []
    for sample in batch.samples:
        if samples and sample.code not in samples: continue
        for assay in sample.assays:
            if assays and assay.filename not in assays: continue
            if panels and assay.panel.code not in panels: continue
            fsa_list.append( (assay, sample.code) )

    cerr('I: number of assays to be processed: %d' % len(assay_list))
    return fsa_list

