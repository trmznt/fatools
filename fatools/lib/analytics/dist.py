# genetic distance matrix
#
#

import numpy as np

def simple_distance( genotable ):
    """ return proportion of difference alleles """

    n = len(genotable)
    m = np.zeros( (n,n) )
    g = genotable.values
    v = []
    l = len(g[0])

    for i in range(0, n):
        for j in range(i, n):
            d = sum( 1 if x != y else 0 for (x,y) in zip( g[i], g[j]) )
            m[i,j] = m[j,i] = d/l
            v.append( d )

    return (m, v)


class DistanceMatrix(object):
    """ this class holds distance matrix result
    """

    def __init__(self, haplotype_sets):
        self._haplotype_sets = haplotype_sets
        self.H = None   # haplotypes as pandas DataFrame
        self.M = None   # the distance matrix
        self.S = None   # set index


    def calculate_distance(self, dfunc):

        for idx, hs in enumerate(self._haplotype_sets):
            if hs.total_samples <= 0:
                continue
            if self.H is None:
                self.L = [idx] * hs.total_samples
                self.H = hs.haplotype_df
                continue
            self.L += [idx] * hs.total_samples
            self.H = self.H.append(hs.haplotype_df)

        (m, v) = dfunc(self.H)
        self.M = m
        self.V = v


    @property
    def total_samples(self):
        return len(self.L)



def get_distance_matrix(haplotype_sets, dfunc = simple_distance):
    """ returning a distance matrix object
        dfunc is a function that returns (matrix, values)
    """

    dm = DistanceMatrix(haplotype_sets)
    dm.calculate_distance( dfunc )

    return dm

