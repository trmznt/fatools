

import numpy as np
from pandas import DataFrame
from scipy.stats import wilcoxon, kruskal

def summarize_he( analytical_sets ):

    results = {}
    he = {}

    for analytical_set in analytical_sets:
        he[analytical_set.label] = calculate_he(analytical_set.allele_df)

    he_df = DataFrame( results )
    if len(analytical_sets) == 2:
        # use Mann-Whitney / Wilcoxon test
        results['test'] = 'Wilcoxon test (paired)'
        results['stats'] = wilcoxon( he_df[labels[0], labels[1]])

    elif len(analytical_sets) > 2:
        # use Kruskal Wallis
        results['test'] = 'Kruskal-Wallis test'
        results['stats'] = kruskal( * [he_df[x] for x in labels])

    results['data'] = he_df

    return results


def calculate_he(allele_df):
    pass

