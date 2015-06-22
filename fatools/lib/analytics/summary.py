
from collections import defaultdict

import numpy as np


def summarize_alleles( analytical_sets ):

    summaries = {}

    for analytical_set in analytical_sets:

        summary = summarize_allele_df( analytical_set.allele_df.df )
        summaries[analytical_set.label] = summary

    return summaries



def summarize_allele_df( allele_df ):
    """ return a tuple of (dict, dict):
        1dict: alleles: [ (allele, freq, count, mean_height, min_size, max_size, delta), ...]
        2dict: marker: ( [ size, ...], [ height, ....] )
    """

    allele_list = defaultdict(list)
    marker_list = defaultdict(lambda x = None: ([], []))
    grouped = allele_df.groupby( ['marker_id', 'value'] )

    for (marker_id, allele), df in grouped:

        allele_list[marker_id].append(
            (allele, len(df), np.mean( df['height'] ), min(df['size']), max(df['size']),
                list(df['sample_id']), np.mean( df['size'] ))
        )

        marker_list[marker_id][0].extend( df['size'] )
        marker_list[marker_id][1].extend( df['height'] )
            

    # calculate other stuff

    results = {}

    for marker_id in allele_list:
        alleles = allele_list[marker_id]
        total_allele = sum( x[1] for x in alleles )
        allele_params = [
            (allele, count/total_allele, count, mean_height, min_size, max_size,
                max_size - min_size, sample_ids, mean_size )
            for (allele, count, mean_height, min_size, max_size, sample_ids, mean_size )
            in alleles
        ]

        delta_status = check_delta( allele_params)

        results[marker_id] = dict(
            code = marker_id,
            unique_allele = len(allele_params),
            total_allele = total_allele,
            alleles = allele_params,
            delta_status = delta_status,
            items = marker_list[marker_id])

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

