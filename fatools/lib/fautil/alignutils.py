
from fatools.lib.utils import cerr, cout
from fatools.lib.fautil.dpalign import dp

import numpy as np
import attr
import math


@attr.s
class AlignResult(object):
    score = attr.ib()
    msg = attr.ib()
    dpresult = attr.ib()
    method = attr.ib()
    initial_pairs = attr.ib(default=None)


@attr.s
class DPResult(object):
    dpscore = attr.ib()
    rss = attr.ib()
    z = attr.ib()
    sized_peaks = attr.ib()

    @property
    def ztranspose(self):
        if self.sized_peaks:
            sizes = []
            rtimes = []
            for size, allele in self.sized_peaks:
                sizes.append(size)
                rtimes.append(allele.rtime)
            zres = estimate_z(sizes, rtimes, len(self.z) - 1)
            return zres.z
        return []


@attr.s
class ZResult(object):
    z = attr.ib()
    rss = attr.ib()
    f = attr.ib()


class PeakPairs(object):

    def __init__(self, peak_pairs):
        self.pairs = peak_pairs
        self.r2s = {}
        self.s2r = {}
        for (rtime, size) in self.pairs:
            self.r2s[rtime] = size
            self.s2r[size] = rtime


def estimate_z( x, y, degree = 3 ):
    """ estimate z and rss
            x: peak rtime
            y: standard sizes
        return (z, rss, func)

        y ~ f(x) where f = poly1d(z)
        rss ~ SUM( (f(x) - y)**2 ) for all (x,y)
    """
    z = np.polyfit( x, y, degree )
    p = np.poly1d( z )
    y_p = p(x)
    rss = ( (y_p - y) ** 2 ).sum()

    return ZResult(z, rss, p)


def generate_similarity( peaks ):

    rfus = [ p.rfu for p in peaks ]

    # use the 2nd highest peaks since the 1st or 2nd may be noises
    highest_rfu = list(sorted(rfus, reverse=True))[2]

    N = len(rfus)
    similarity = list( [ (np.log10( rfu/ highest_rfu ) + N) / N if rfu < highest_rfu
                            else 1.0
                            for rfu in rfus ] )

    for (p, score) in zip(peaks, similarity):
        p.qscore = score
    return similarity


def pair_sized_peaks( peaks, peak_pairs ):
    """ translate (rtime, size) to (size, peak)

        return: [ (size, peak), ...]
    """

    d_peaks = {}
    for p in peaks:
        d_peaks[p.rtime] = p

    sized_peaks = []
    for (bpsize, rtime) in peak_pairs:
        sized_peaks.append( (bpsize, d_peaks[rtime]) )

    sized_peaks.sort( key = lambda x: x[0] )
    return sized_peaks


def pair_f(f, rtimes, std_sizes, similarity, deviation=False):
    """ match rtimes to std_sizes

        return: [ (rtime, size), ... ] or
                [ (rtime, size, f(rtime), dev), ... ]

    """

    rtimes = list(reversed(rtimes))
    std_sizes = list(reversed(std_sizes))
    similarity = list(reversed(similarity))

    S = generate_scores( std_sizes, rtimes, similarity, f )

    result = dp(S, -5e-3)

    matches = result['matches']
    aligned_peaks = [ (rtimes[j], std_sizes[i]) for i, j in matches ]

    if not deviation:
        return aligned_peaks

    peak_pairs = []
    for (rtime, size) in aligned_peaks:
        rtime_size = f(rtime)
        peak_pairs.append( (rtime, size, rtime_size, (size - rtime_size)**2) )

    return peak_pairs


def generate_scores_xxx(sizes, rtimes, func, tolerance = 4):
    """ return a numpy matrix for scoring peak similarity
            func -> polynomial fit funcs
            size[bp] = func(rtime[sec])

        Score matrix is

                  peak1   peak2   peak3
        ladder1
        ladder2
        ladder3

        S[ladder][peak] = 1 if ladder & peak are similar

    """
    M = np.zeros( (len(sizes), len(rtimes)), dtype='d' )

    _TOL = 0.001
    cutoff = tolerance * math.sqrt( -1.0 * math.log(_TOL))
    ladder_N = len(sizes)

    row = 0
    col = 0

    for peak in rtimes:
        size = func(peak)
        for ladder in sizes:
            M[row][col] = math.exp( - ((size - ladder)/(tolerance))**2  / 2 )
            row += 1
        col += 1
        row = 0

    return M


def generate_scores(sizes, rtimes, similarity, func, tolerance = 4):
    """ return a numpy matrix for scoring peak similarity
            func -> polynomial fit funcs
            size[bp] = func(rtime[sec])

        Score matrix is

                  peak1   peak2   peak3
        ladder1
        ladder2
        ladder3

        S[ladder][peak] = 1 if ladder & peak are similar

    """
    M = np.zeros( (len(sizes), len(rtimes)), dtype='d' )

    _TOL = 0.001
    cutoff = tolerance * math.sqrt( -1.0 * math.log(_TOL))
    ladder_N = len(sizes)

    for c in range(0, len(rtimes)):
        size = func( rtimes[c] )
        for r in range(0, len(sizes)):
            M[r][c] = similarity[c] * math.exp( - ((size - sizes[r])/(tolerance))**2  / 2 )

    return M


def plot(rtimes, sizes, z, peak_pairs):
    """ plot rtimes, sizes, z and peak pairs
    """

    from matplotlib import pyplot as plt

    fig = plt.figure()
    max_x = max(rtimes) + 100
    min_x = min(rtimes) - 100
    range_x = max_x - min_x
    x = np.linspace(min_x, max_x, range_x+1)
    f = np.poly1d(z)
    y = f(x)

    plt.plot(x, y)
    plt.hlines(sizes, min_x, max_x)
    plt.vlines(rtimes, 0, max(sizes) + 50)

    plt.scatter( [x[0] for x in peak_pairs], [x[1] for x in peak_pairs])
    plt.ylim([0, max(sizes) + 50])
    plt.xlim([min_x, max_x])
    plt.show()


def align_dp( rtimes, sizes, similarity, z, rss, order = 3):
    """ align ladders with peaks using dynamic programming (global alignment)
        return (dpscore, RSS, Z, ladder_aligned_peaks)
    """

    sizes = list( sorted(sizes, reverse=True) )
    rtimes = list( sorted(rtimes, reverse=True) )

    dpscore = -1

    while True:

        S = generate_scores( sizes, rtimes, similarity, np.poly1d(z))

        result = dp(S, -5) #-5e-3)

        cur_dpscore = result['D'][-1][-1]
        matches = result['matches']

        aligned_peaks = [ (sizes[i], rtimes[j]) for i, j in matches ]

        # realign

        std_size, peak_sizes = zip(*aligned_peaks)
        cur_zres = estimate_z( peak_sizes, std_size, order )

        if cur_dpscore < dpscore:
            cerr('W: dynamic programming did not converge!!')
            break

        if cur_dpscore == dpscore:
            break

        z = cur_zres.z
        rss = cur_zres.rss
        dpscore = cur_dpscore
        sized_peaks = aligned_peaks

    return DPResult(dpscore, rss, z, sized_peaks)
