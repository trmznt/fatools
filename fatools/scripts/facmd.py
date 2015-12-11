

import sys, argparse, yaml, csv, transaction
from fatools.lib.utils import cout, cerr, get_dbhandler
from fatools.lib import params
from fatools.lib.const import assaystatus, peaktype
from fatools.lib.fautil import algo


def init_argparser(parser=None):

    if parser:
        p = parser
    else:
        p = argparse.ArgumentParser('facmd')

    p.add_argument('--sqldb', default=False,
            help = 'SQLite3 database filename')

    p.add_argument('--fsdb', default=False,
            help = 'directory for filesystem-based database')

    p.add_argument('--clear', default=False, action='store_true',
            help = 'clear / remove all peaks from assay')

    p.add_argument('--scan', default=False, action='store_true',
            help = 'scanning assay for peaks')

    p.add_argument('--preannotate', default=False, action='store_true',
            help = 'preannotate assay for overlapping peaks, stutter and broad peaks')

    p.add_argument('--alignladder', default=False, action='store_true',
            help = 'align ladder peaks with standard size')

    p.add_argument('--call', default=False, action='store_true',
            help = 'calling peaks (determining the sizes of peaks)')

    p.add_argument('--bin', default=False, action='store_true',
            help = 'binning peaks')

    p.add_argument('--postannotate', default=False, action='store_true',
            help = 'post annotate peaks')

    p.add_argument('--listpeaks', default=False, action='store_true',
            help = 'list all peaks')

    p.add_argument('--listassay', default=False, action='store_true',
            help = 'list assay information')

    p.add_argument('--showtrace', default=False, action='store_true',
            help = 'show trace as a plot')

    p.add_argument('--findpeaks', default=False, action='store_true',
            help = 'only find peaks')

    p.add_argument('--setallele', default=False, action='store_true',
            help = 'set allele type')


    p.add_argument('--batch', default=False,
            help = 'batch code')

    p.add_argument('--sample', default=False,
            help = 'sample code')

    p.add_argument('--assay', default=False,
            help = 'assay filename')

    p.add_argument('--marker', default=False,
            help = 'marker code')

    p.add_argument('--commit', default=False, action='store_true',
            help = 'commit to database')

    p.add_argument('--outfmt', default='text',
            help = 'output format, either text or tab')

    p.add_argument('--outfile', default='-',
            help = 'output filename')

    p.add_argument('--peakcachedb', default=False,
            help = 'peak cache db filename')

    p.add_argument('--method', default='',
            help = 'spesific method or algorithm to use')

    p.add_argument('--value', default='',
            help = 'bin value of alleles')

    p.add_argument('--totype', default='',
            help = 'new type of alleles')

    p.add_argument('--fromtype', default='',
            help = 'original type of alleles')


    p.add_argument('--excluded_peaks')

    p.add_argument('--stutter_ratio', default=0, type=float)

    p.add_argument('--stutter_range', default=0, type=float)

    p.add_argument('--force', default=False, action='store_true',
            help = 'force the method (even if need short-cutting)' )

    p.add_argument('--test', default=False, action='store_true',
            help = 'just testing, not need to commit to database')

    p.add_argument('-y', default=False, action='store_true',
            help = 'say yes to all interactive questions')

    p.add_argument('--abort', default=False, action='store_true',
            help = 'abort for any warning')


    p.add_argument('--showladderpca', default=False, action='store_true',
            help = 'show PCA plot for ladder peaks')

    p.add_argument('--showz', default=False, action='store_true',
            help = 'show Z plot for ladder peaks')

    return p


def main(args):

    if args.commit:
        with transaction.manager:
            do_facmd(args)
            cerr('** COMMIT to database **')
    else:
        cerr('WARNING ** running without database COMMIT! All changes will be discarded!')
        if not ( args.test or args.y ):
            keys = input('Do you want to continue [y/n]? ')
            if not keys.lower().strip().startswith('y'):
                sys.exit(1)
        do_facmd(args)


def do_facmd(args, dbh=None):

    if dbh is None:
        dbh = get_dbhandler(args)

    executed = 0
    if args.clear is not False:
        do_clear( args, dbh )
        executed += 1
    if args.findpeaks is not False:
        do_findpeaks( args, dbh )
        executed += 1
    if args.scan is not False:
        do_scan(args, dbh)
        executed += 1
    if args.preannotate is not False:
        do_preannotate(args, dbh)
        executed += 1
    if args.alignladder is not False:
        do_alignladder(args, dbh)
        executed += 1
    if args.call is not False:
        do_call(args, dbh)
        executed += 1
    if args.bin is not False:
        do_bin(args, dbh)
        executed += 1
    if args.postannotate is not False:
        do_postannotate(args, dbh)
        executed += 1
    if args.setallele is not False:
        do_setallele(args, dbh)
        executed += 1
    if args.showladderpca is not False:
        do_showladderpca( args, dbh )
        executed += 1
    if args.listassay is not False:
        do_listassay( args, dbh )
        executed += 1
    if args.listpeaks is not False:
        do_listpeaks( args, dbh )
        executed += 1
    if args.showtrace is not False:
        do_showtrace( args, dbh )
        executed += 1
    if args.showz is not False:
        do_showz( args, dbh )
        executed += 1
    if executed == 0:
        cerr('WARN - unknown command, nothing to do!')
    else:
        cerr('INFO - executed %d command(s)' % executed)


def do_clear( args, dbh ):

    cerr('Clearing peaks...')

    assay_list = get_assay_list( args, dbh )
    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('Clearing sample: %s assay %s [%d/%d]' %
                (sample_code, assay.filename, counter, len(assay_list)))

        assay.clear()
        counter += 1


def do_scan( args, dbh ):

    cerr('I: Scanning peaks...')

    scanning_parameter = params.Params()
    assay_list = get_assay_list( args, dbh )

    if args.peakcachedb:
        import leveldb
        peakdb = leveldb.LevelDB(args.peakcachedb, create_if_missing=False)
    else:
        peakdb = None

    if args.method:
        scanning_parameter.ladder.method = args.method
        scanning_parameter.nonladder.method = args.method

    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('I: [%d/%d] - Scanning: %s | %s' %
                (counter, len(assay_list), sample_code, assay.filename))

        assay.scan( scanning_parameter, peakdb = peakdb )
        counter += 1


def do_preannotate( args, dbh ):

    cerr('I: Preannotating peaks...')

    scanning_parameter = params.Params()
    assay_list = get_assay_list( args, dbh )

    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('I: [%d/%d] - Preannotating: %s | %s' %
                (counter, len(assay_list), sample_code, assay.filename))

        assay.preannotate( scanning_parameter )
        counter += 1


def do_alignladder( args, dbh ):

    cerr('Aligning ladders...')

    assay_list = get_assay_list( args, dbh )
    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('I: [%d/%d] - Aligning: %s | %s' %
                (counter, len(assay_list), sample_code, assay.filename))
        (dpscore, rss, no_of_peaks, no_of_ladders, qcscore, remarks, method) = assay.alignladder( args.excluded_peaks, force_mode = args.force )
        if qcscore < 0.9:
            msg = 'W! low ladder QC'
        else:
            msg = 'I:'
        cerr( '%s [%d/%d] - Score %3.2f %4.2f %5.2f %d/%d %s for %s | %s'
                % ( msg, counter, len(assay_list),
                    qcscore, dpscore, rss, no_of_peaks, no_of_ladders,
                    method, sample_code, assay.filename) )
        if remarks:
            cerr('%s - %s' % (msg, ' | '.join(remarks)))
        if qcscore != 1.0 and args.abort:
            sys.exit(1)

        counter += 1


def do_call(args, dbh):

    cerr('I: Calling peaks...')

    scanning_parameter = params.Params()

    assay_list = get_assay_list( args, dbh )
    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('I: [%d/%d] - Calling: %s | %s' %
                (counter, len(assay_list), sample_code, assay.filename))
        assay.call( scanning_parameter.nonladder )
        counter += 1


def do_bin(args, dbh):

    cerr('I: Binning peaks...')

    scanning_parameter = params.Params()

    if args.marker:
        markers = [ dbh.get_marker( code ) for code in args.marker.split(',') ]
    else:
        markers = None


    assay_list = get_assay_list( args, dbh )
    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('I: [%d/%d] - Binning: %s | %s' %
                (counter, len(assay_list), sample_code, assay.filename))
        assay.bin( scanning_parameter.nonladder, markers )
        counter += 1


def do_postannotate(args, dbh):

    cerr('I: Post-annotating peaks...')

    scanning_parameter = params.Params()

    if args.marker:
        markers = [ dbh.get_marker( code ) for code in args.marker.split(',') ]
    else:
        markers = None

    if args.stutter_ratio > 0:
        scanning_parameter.nonladder.stutter_ratio = args.stutter_ratio
    if args.stutter_range > 0:
        scanning_parameter.nonladder.stutter_range = args.stutter_range


    assay_list = get_assay_list( args, dbh )
    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('I: [%d/%d] - Post-annotating: %s | %s' %
                (counter, len(assay_list), sample_code, assay.filename))
        assay.postannotate( scanning_parameter.nonladder, markers )
        counter += 1



def do_findpeaks( args, dbh ):

    import leveldb
    from fatools.lib import params

    cerr('Finding and caching peaks...')

    if not args.peakcachedb:
        cexit('ERR - please provide cache db filename')

    # opening LevelDB database
    peakdb = leveldb.LevelDB(args.peakcachedb)

    scanning_parameter = params.Params()
    assay_list = get_assay_list( args, dbh )

    if args.method:
        scanning_parameter.ladder.method = args.method
        scanning_parameter.nonladder.method = args.method

    channel_list = []
    counter = 1
    cerr('', nl=False)
    for (assay, sample_code) in assay_list:
        cerr('\rI: [%d/%d] processing assay' % (counter, len(assay_list)), nl=False)
        for c in assay.channels:
            if c.marker.code == 'ladder':
                params = scanning_parameter.ladder
            else:
                params = scanning_parameter.nonladder
            channel_list.append( (c.tag(), c.data, params) )
        counter += 1
    cerr('')

    do_parallel_find_peaks( channel_list, peakdb )

    #peakdb.close()



def do_setallele( args, dbh ):

    marker_codes = args.marker.split(',')
    marker_ids = [ dbh.get_marker(code).id for code in marker_codes ]
    bin_values = [ int(x) for x in args.value.split(',') ]

    totype = getattr(peaktype, args.totype)

    assay_list = get_assay_list( args, dbh )
    counter = 1
    for (assay, sample_code) in assay_list:
        for c in assay.channels:
            if marker_ids and c.marker_id in marker_ids:
                for allele in c.alleles:
                    if not allele.bin in bin_values:
                        continue
                    if args.fromtype and allele.type != args.fromtype:
                        continue
                    allele.type = totype
                    cerr('I: - setting allele %d marker %s for sample %s' %
                            (allele.bin, c.marker.label, sample_code))



def do_showladderpca( args, dbh ):

    assay_list = get_assay_list( args, dbh )
    counter = 1
    for (assay, sample_code) in assay_list:
        cerr('Showing ladder PCA for  sample: %s assay %s [%d/%d]' %
                (sample_code, assay.filename, counter, len(assay_list)))
        assay.showladderpca()


def do_listassay( args, dbh ):

    assay_list = get_assay_list( args, dbh )

    if args.outfile != '-':
        out_stream = open(args.outfile, 'w')
    else:
        out_stream = sys.stdout
    for (assay, sample_code) in assay_list:
        printout_assay( assay, outfile = out_stream, fmt = args.outfmt )


def do_listpeaks( args, dbh ):

    assay_list = get_assay_list( args, dbh )
    if args.marker:
        markers = [ dbh.get_marker( code ) for code in args.marker.split(',') ]
    else:
        markers = None

    if markers:
        cerr('Markers: %s' % ','.join( m.code for m in markers ))

    for (assay, sample_code) in assay_list:
        cout('Sample: %s assay: %s' % (sample_code, assay.filename))
        for channel in assay.channels:
            if markers and channel.marker not in markers:
                continue
            cout('Marker => %s | %s [%d]' % (channel.marker.code, channel.dye,
                    len(channel.alleles)))
            for p in channel.alleles:
                cout( '  ' + str(p) )

def do_showtrace( args, dbh ):

    assay_list = get_assay_list( args, dbh )

    from matplotlib import pylab as plt

    for (assay, sample_code) in assay_list:
        peaks = []
        for c in assay.channels:
            plt.plot( c.raw_data )
            peaks +=  list(c.alleles)

        for p in peaks:
            plt.plot( p.rtime, p.height, 'r+' )

        plt.show()


def do_showz( args, dbh ):

    assay_list = get_assay_list( args, dbh )

    from matplotlib import pylab as plt
    import numpy as np

    for (assay, sample_code) in assay_list:
        ladder_peaks = list(assay.ladder.alleles)
        z = assay.z

        peak_pairs = [ (x.rtime, x.size) for x in ladder_peaks ]

        x = np.linspace( peak_pairs[0][0], peak_pairs[-1][0] + 200, 100 )
        f = np.poly1d(z)
        y = f(x)

        print(' => Z: ', z)
        for p in ladder_peaks:
            print(' => %6d -> %6.2f | %4d | %5.2f' % ( p.rtime, f(p.rtime), p.size,
                                                        abs( f(p.rtime) - p.size )))

        plt.plot(x,y)
        rtimes = [ x[0] for x in peak_pairs ]
        sizes = [ x[1] for x in peak_pairs ]
        plt.scatter(rtimes, sizes)
        plt.show()

# helpers

def get_assay_list( args, dbh ):

    if not args.batch:
        cerr('ERR - need --batch argument!')
        sys.exit(1)

    batch = dbh.get_batch( args.batch )
    if not batch:
        cerr('ERR - batch %s not found!' % args.batch)
        sys.exit(1)

    samples = []
    if args.sample:
        samples = args.sample.split(',')

    assays = []
    if args.assay:
        assays = args.assay.split(',')

    assay_list = []
    for sample in batch.samples:
        if samples and sample.code not in samples: continue
        for assay in sample.assays:
            if assays and assay.filename not in assays: continue
            assay_list.append( (assay, sample.code) )

    cerr('INFO - number of assays to be processed: %d' % len(assay_list))
    return assay_list


## PRINTOUT

def printout_assay( assay, outfile=sys.stdout, fmt='text' ):

    if fmt == 'tab':
        outfile.write('%s\t%s\t%f\t%f\t%f\t%d\t%d\t%s\n' %
                    ( assay.sample.code, assay.filename,
                    assay.score, assay.dp, assay.rss, assay.ladder_peaks,
                    len(assay.ladder.alleles), assay.method) )
        return ''


    buf = []
    _ = buf.append

    _( 'Assay: %s -- Sample: %s' % (assay.filename, assay.sample.code) )
    if assay.status in ( assaystatus.aligned, assaystatus.called,
                        assaystatus.annotated, assaystatus.binned ):
        _( ' => Score: %3.2f, DP: %5.2f, RSS: %5.2f, N-peak: %d' %
            ( assay.score, assay.dp, assay.rss, assay.ladder_peaks ))

    return '\n'.join( buf )


## parallel word

def do_parallel_find_peaks( channel_list, peakdb ):

    import concurrent.futures, pickle

    cerr('I: Processing channel(s)')
    total = len(channel_list)
    counter = 0
    with concurrent.futures.ProcessPoolExecutor() as executor:
        for (tag, peaks) in executor.map( find_peaks_p, channel_list ):
            peakdb.Put(tag.encode(), pickle.dumps(peaks))
            counter += 1
            cerr('I: [%d/%d] channel %s => %d peak(s)' % (counter, total, tag, len(peaks)))


def find_peaks_p( args ):
    tag, data, param = args

    return (tag, algo.find_raw_peaks(data, param))

