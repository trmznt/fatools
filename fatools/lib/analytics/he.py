

import numpy as np
from pandas import DataFrame
from scipy.stats import wilcoxon, kruskal


def summarize_he( analytical_sets ):

    results = {}
    he = {}

    for analytical_set in analytical_sets:
        he[analytical_set.label] = calculate_he(analytical_set.allele_df)

    he_df = DataFrame( he )
    labels = list(he_df.columns)
    if len(labels) == 2:
        # use Mann-Whitney / Wilcoxon test
        results['test'] = 'Wilcoxon test (paired)'
        results['stats'] = wilcoxon( he_df[labels[0], labels[1]])

    elif len(labels) > 2:
        # use Kruskal Wallis
        results['test'] = 'Kruskal-Wallis test'
        results['stats'] = kruskal( * [he_df[x] for x in labels])

    results['data'] = he_df
    results['mean'] = he_df.mean()
    results['stddev'] = he_df.std()

    return results


def calculate_he(allele_df, adjust=True):
    """ He is calculated using major allele """
    marker_he = {}

    dist = allele_df.dominant_df_distribution
    for marker_id in dist.index.levels[0]:
        total = dist[marker_id].sum()
        he = 1.0 - sum( (x/total)**2 for x in dist[marker_id] )
        if adjust and total > 1:
            he = he *total / (total-1)
        marker_he[marker_id] = he

    return marker_he
