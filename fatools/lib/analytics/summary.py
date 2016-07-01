
from collections import defaultdict
import numpy as np
import sys
from pprint import pprint


def summarize_alleles( analytical_sets ):
    """ returns a dict of:
        set_label1:
            marker_id1:
                code:
                unique_alleles:
                alleles:
                ...
            marker_id2:
                code:
                ...
        set_label2:
            marker_id2:
                code:
                unique_alleles:
                ...
    """

    summaries = {}

    for analytical_set in analytical_sets:

        summary = summarize_allele_df( analytical_set.allele_df )
        summaries[analytical_set.label] = { 'summary': summary,
                                            'colour': analytical_set.colour }

    return summaries



def summarize_allele_df( allele_df ):
    """ return a dict of data containing:
        alleles: [ (allele, freq, count, mean_height, min_size, max_size, delta, items), ...]
            where items is [
    """

    allele_list = defaultdict(list)
    #marker_list = defaultdict(lambda x = None: ([], []))
    #print(allele_df)

    for (marker_id, allele), df in allele_df.grouped_df:

        allele_list[marker_id].append(
            (allele, len(df), np.mean( df['height'] ), min(df['size']), max(df['size']),
                list(df['sample_id']), np.mean( df['size']), (df['size'], df['height']))
        )

        #marker_list[marker_id][0].extend( df['size'] )
        #marker_list[marker_id][1].extend( df['height'] )


    # calculate other stuff

    results = {}

    for marker_id in allele_list:
        alleles = allele_list[marker_id]
        total_allele = sum( x[1] for x in alleles )
        allele_params = [
            (allele, count/total_allele, count, mean_height, min_size, max_size,
                max_size - min_size, sample_ids, mean_size, items )
            for (allele, count, mean_height, min_size, max_size, sample_ids, mean_size, items )
            in alleles
        ]

        delta_status = check_delta(allele_params)

        results[marker_id] = dict(
            code = marker_id,
            unique_allele = len(allele_params),
            total_allele = total_allele,
            alleles = allele_params,
            delta_status = delta_status,
            #items = marker_list[marker_id]
            )

    return results


def check_delta( alleles ):
    """ return True if allele bin is 1 bp adjacent to prev or nex allele bin """

    # check if only single allele
    if len(alleles) <= 1:
        return [ True ]

    threshold = 1

    delta_status = []
    if alleles[1][0] - alleles[0][0] <= threshold:
        delta_status.append( False )
    else:
        delta_status.append( True )
    for i in range(1, len(alleles) - 1):
        if (    alleles[i][0] - alleles[i-1][0] <= threshold or
                alleles[i+1][0] - alleles[i][0] <= threshold ):
            delta_status.append( False )
        else:
            delta_status.append( True )
    if alleles[-2][0] - alleles[-1][0] == 1:
        delta_status.append( False )
    else:
        delta_status.append( True )

    return delta_status


def summarize_bins( analytical_sets ):
    """ return bin summary for each marker """

    allele_summaries = summarize_alleles( analytical_sets )

    marker_summaries = defaultdict( lambda x=None: defaultdict( list ) )
    # marker_summaries[marker_id] = { bin1: [ ?, ..], bin2:

    # collect items (sizes) for each marker
    for label in allele_summaries:
        for (marker_id, allele_summary) in allele_summaries[label]['summary'].items():
            pprint(allele_summary)
            for allele_params in allele_summary['alleles']:
                bin_value = allele_params[0]
                size_items = allele_params[9][0]
                marker_summaries[marker_id][bin_value].extend( size_items )
                marker_summaries[marker_id][bin_value].sort()

    pprint(marker_summaries)
    bin_summaries = {}

    # process marker summary to obtain new bin paramater
    for marker_id, marker_summary in marker_summaries.items():
        empirical_bins = {}
        for (allele_bin, sizes) in marker_summary.items():
            # not found
            percentiles = np.percentile( sizes, [ 25, 75 ] )
            empirical_bins[ allele_bin ] = [ int(allele_bin),
                                                round(float(np.mean(sizes)), 3),
                                                round(float(percentiles[0]), 3),
                                                round(float(percentiles[1]), 3) ]
        bin_summaries[marker_id] = empirical_bins

    return bin_summaries


def plot_alleles( allele_reports, filename, rfu_height=True, dbh=None ):

    from matplotlib import pyplot as plt
    from matplotlib.figure import Figure
    from matplotlib.backends.backend_agg import FigureCanvasAgg as FigureCanvas
    from matplotlib.ticker import MultipleLocator


    marker_ids = set()
    for (label, allele_report) in allele_reports.items():
        marker_ids.update( allele_report['summary'].keys() )

    m = len(marker_ids)
    fig = plt.figure( figsize = (21, 4 * m), dpi=600 )

    axes = {}
    axhlines = set()
    binsets = {}

    for idx, allele_report in enumerate(allele_reports.values(), 0):
        pprint(allele_report)
        colour = allele_report['colour']
        for (marker_id, summary) in allele_report['summary'].items():
            if marker_id in axes:
                ax = axes[marker_id]
                bins = binsets[marker_id]
            else:
                ax = fig.add_subplot( m, 1, len(axes) + 1 )
                axes[marker_id] = ax
                bins = binsets[marker_id] = set()

            x = []
            y = []

            for allele_params in summary['alleles']:
                bins.add( allele_params[0] )
                data = allele_params[9]

                x += list(data[0])
                y += list(data[1])

            if rfu_height:
                max_y = max(y)
                y = [ (idx + h/max_y) for h in y ]
            else:
                y = [ (idx + 1) for h in y ]
            idxs = [ idx ] * len(x)

            ax.vlines( x, idxs, y, colors = colour )
            ax.vlines( x, 0, -0.05, colors = 'k')
            if not (ax, idx) in axhlines:
                # just to make sure we don't duplicate these lines
                ax.axhline( idx, color='#aaaaaa' )
                axhlines.add( (ax, idx) )


    for (marker_id, ax) in axes.items():

        ax.get_xaxis().set_tick_params( which='both', direction='out' )
        ax.get_yaxis().set_tick_params( which='both', direction='out' )
        ax.get_xaxis().set_minor_locator( MultipleLocator(1) )
        ax.get_xaxis().set_ticks( sorted(list(binsets[marker_id])) )
        ax.get_yaxis().set_minor_locator( MultipleLocator(0.25))
        ax.get_yaxis().set_major_locator( plt.NullLocator() )

        for label in ax.get_xticklabels():
            label.set_size( 'xx-small' )
        for label in ax.get_yticklabels():
            label.set_size( 'xx-small' )

        if dbh:
            ax.set_ylabel( dbh.get_marker_by_id(marker_id).label )
        else:
            ax.set_ylabel( marker_id )
        ax.set_ylim(-0.05)
        #ax.set_xlim(min(data[0]), max(data[0]))
        ax.set_xlim(auto = True)


    fig.tight_layout()

    fig.savefig( filename )
    plt.close()



def summarize_haplotypes(haplotype_sets):

    summaries = {}
    H = {}  # dictionary containing all haplotype_sets

    for hs in haplotype_sets:
        summaries[hs.label] = summarize_haplotype_df(hs.haplotype_df)

    # check all haplotypes, finding singletons (or dualtons, tripletons, etc)


def summarize_haplotype_df(haplotype_df):

    pass
