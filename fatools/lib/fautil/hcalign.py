
# hclustalign.py
# perform hierarchical clustering

from fatools.lib.utils import cerr, cverr, cexit
from fatools.lib import const
from fatools.lib.fautil.alignutils import ( estimate_z, AlignResult, align_dp,
            pair_sized_peaks, generate_similarity )

from scipy.cluster.hierarchy import dendrogram, linkage, fcluster
from collections import defaultdict
from functools import reduce
import operator
import numpy as np
import attr

@attr.s
class T(object):
    p = attr.ib()       # original peaks
    z = attr.ib()       # z matrix, containing tree
    c = attr.ib()       # list of peaks contained in each node/branch

    def cluster(self, idx):
        if idx == None:
            return []
        if idx < len(self.p):
            return [ self.p[idx] ]
        return self.c[ idx - len(self.p) ]

    def binodes(self, idx):
        if idx < 0:
            idx = len(self.p) + len(self.z) + idx
        if idx < len(self.p):
            return (idx, None)
        n = self.z[ idx - len(self.p) ]
        return (int(n[0]), int(n[1]))

    def biclusters(self, idx):
        pass

    def bitrees(self, idx):
        """ return (n1, c1, n2, c2) """
        n1, n2 = self.binodes( idx )
        return n1, self.cluster(n1), n2, self.cluster(n2)



def generate_tree( p ):

    z = linkage(p, 'average')
    l = len(p)
    c = []
    for n in z:
        cc = []
        n1, n2 = int(n[0]), int(n[1])
        if n1 >= l:
            cc += c[n1 - l]
        else:
            cc.append( n1 )
        if n2 >= l:
            cc += c[n2 - l]
        else:
            cc.append( n2 )
        cc.sort()
        c.append( cc )


    return T(p, z, c)


def generate_cluster(T, k):

    grouping = fcluster(T.z, k, criterion='maxclust')

    groups = defaultdict(list)
    for i, e in enumerate(grouping):
        groups[e].append( T.p[i][0] )

    clusters = sorted( list(groups.values()), key = lambda x: x[0] )
    cverr(3, str(clusters))

    return clusters


def align_hc( peaks, ladder):
    """ peaks: list of rtime, in ascending order
        ladders: list of size from ladders, in ascending order

        returns: (score, msg, result, method)
    """

    #import pprint; pprint.pprint(peaks)

    # generate P for ladder
    if 'C' not in ladder:
        if 'T' not in ladder:
            ladder['T'] = generate_tree( [ (n,0) for n in ladder['sizes'] ] )
        ladder['C'] = generate_cluster(ladder['T'], ladder['k'])
    ladder_clusters = ladder['C']
    ladder_sizes = ladder['sizes']

    P = generate_tree( [ (n.rtime, 0) for n in peaks ] )
    peak_clusters = generate_cluster( P, ladder['k'] )

    # generate cluster should use so-called balance tree

    print(peak_clusters)

    if len(peak_clusters[-1]) == 1:
        if len( reduce(operator.add, peak_clusters ) ) > len(ladder_sizes):
            del peak_clusters[-1]
            #del peaks[-1]
    if len(peak_clusters[0]) == 1:
        if len( reduce(operator.add, peak_clusters ) ) > len(ladder_sizes):
            del peak_clusters[0]
            #del peaks[0]

    if len(peak_clusters) < ladder['k']:
        P = generate_tree( [ (n, 0) for n in reduce(operator.add, peak_clusters) ] )
        peak_clusters = generate_cluster(P, ladder['k'])

    # short cut, in case we have good high quality peaks
    if sum( len(c) for c in peak_clusters ) == len(ladder_sizes):
        hq_peaks = sum(peak_clusters, [])
        #hq_pairs = zip(hq_peaks, ladder_sizes)
        zres = estimate_z(hq_peaks, ladder_sizes)
        dp_result = align_dp( hq_peaks, ladder_sizes, [1.0] * len(hq_peaks),
                                    zres.z, zres.rss )
        dp_result.sized_peaks = pair_sized_peaks(peaks, dp_result.sized_peaks)
        score, msg = ladder['qcfunc']( dp_result, method = 'relax')
        if score > 0.9:
            return AlignResult(score, msg, dp_result, const.alignmethod.hcm_strict)

    #print(">>> clusters:\n", peak_clusters)
    cluster_pairings, expected_missing = align_clusters( peak_clusters,
            ladder_clusters )

    #print(">>> cluster pairs:\n", cluster_pairings)
    # check each cluster pairing

    initial_pairs = []
    for pairs in cluster_pairings:
        if is_good_pairing(pairs):
            initial_pairs.extend( pairs )
        else:
            cverr(3, '>> this pairings is not included:\n%s' % pairs)

    cverr(3, '>> initial pairs:\n%s' % initial_pairs)

    if not initial_pairs:
        return AlignResult(-1, 'E: no initial pairs defined!', None, None)

    # try to dp align the initial pairs as a shortcut for good sample or peaks

    rtimes, sizes = zip( *initial_pairs )
    zres = estimate_z(rtimes, sizes)

    dp_result = align_dp( [p.rtime for p in peaks], ladder_sizes,
                            generate_similarity(peaks), zres.z, zres.rss )
    dp_result.sized_peaks = pair_sized_peaks(peaks, dp_result.sized_peaks)
    score, msg = ladder['qcfunc']( dp_result, method = 'strict')
    if score > 0.9:
        return AlignResult(score, msg, dp_result, const.alignmethod.hcm_strict)

    return AlignResult(-1, 'ERR: alignment needs minimization', None, None,
                initial_pairs=initial_pairs)


def align_clusters( peak_clusters, ladder_clusters, anchor=None ):

    pairs = []
    expected_missing = 0
    last_i = -1
    for i in range(len(ladder_clusters)):
        if len(peak_clusters) < i:
            expected_missing += len( ladder_clusters[i] )
            break
        c_peak = peak_clusters[i]
        c_ladder = ladder_clusters[i]
        if len(c_peak) == len(c_ladder):
            if last_i < 0 or i - last_i == 1:
                pairs.append ( list(zip( c_peak, c_ladder)) )
                last_i = i
        elif len(c_peak) < len(c_ladder):
            expected_missing += len(c_ladder) - len(c_peak)

    return pairs, expected_missing


def is_good_pairing( pairs ):

    rtimes, bpsizes = zip( *pairs )
    zres = estimate_z(rtimes, bpsizes, 1)
    #f = np.poly1d(zres.z)
    #for (rtime, bpsize) in pairs:
    #    print(':', rtime, bpsize, (f(rtime)-bpsize)**2)
    #print(rss)
    # check if total rss is less than the threshold
    return (zres.rss < (0.5 * len(bpsizes)))


def align_clusters2( peak_clusters, ladder_clusters, anchor):
    pairs = []
    pairs.append( (peak_clusters[anchor-1][-1], ladder_clusters[anchor-1][-1]) )
    pairs.extend( zip( peak_clusters[anchor], ladder_clusters[anchor]) )
    return pairs, 0

