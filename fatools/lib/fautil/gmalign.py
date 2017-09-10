
import numpy as np
from scipy.optimize import minimize, differential_evolution

from fatools.lib.utils import cerr
from fatools.lib import const
from fatools.lib.fautil.alignutils import (estimate_z, pair_f, align_dp,
            pair_sized_peaks, DPResult, AlignResult, generate_similarity, plot)

class ZFunc(object):

    def __init__(self, peaks, sizes, anchor_pairs):
        """
        """

class ZFunc(object):

    def __init__(self, peaks, sizes, anchor_pairs, estimate=False):
        """
        peaks, sizes and anchor_pairs must be in ascending order
        """
        self.set_sizes(sizes)
        self.peaks = list(sorted(peaks))
        self.rtimes = [ p.rtime for p in self.peaks]
        self.anchor_rtimes = [ a[0] for a in anchor_pairs]
        self.anchor_sizes = [ a[1] for a in anchor_pairs]
        self.anchor_pairs = anchor_pairs
        self.penalty = 4 if estimate else 1
        #if estimate:
        self.similarity = generate_similarity( self.peaks )
        #else:
        #self.similarity = [ 1.0 ] * len(self.peaks)

    def set_sizes(self, sizes):
        self.sizes = list(sorted(sizes))


    def get_initial_z(self):
        """
        return (initial_z, initial_rss) based on anchor pairs
        """

        # check which order we want to use for initial Z
        if (    (self.anchor_sizes[-1] - self.anchor_sizes[0]) >
                0.2 * (self.sizes[-1] - self.sizes[0])
                and len(self.anchor_pairs) >= 5
            ):
                orders = [1, 2]
        else:
                orders = [1]

        zresults = []
        for order in orders:
            zresult = estimate_z( self.anchor_rtimes, self.anchor_sizes, order )

            zres = align_dp(
                    self.rtimes, self.sizes, zresult.z, zresult.rss, order)

            zresults.append( zres )

        zresults.sort( key = lambda x: x.rss)
        return zresults[0]


    def get_pairs(self, z):

        f = np.poly1d(z)
        pairings = pair_f(f, self.rtimes, self.sizes, self.similarity, deviation=True)

        pairs = []
        rss = 0
        for (rtime, bpsize, rsize, err) in pairings:
            pairs.append( (rtime, bpsize) )
            rss += err

        return pairs, rss


    def get_sized_peaks(self, peak_pairs):

        return pair_sized_peaks(self.peaks, peak_pairs)


    def __call__(self, z):
        """
        z is array for polynomial curve
        return rss
        """

        # prepare z function
        f = np.poly1d(z)
        pairs = pair_f(f, self.rtimes, self.sizes, self.similarity, deviation=True)


        # calculate anchors
        rss = 0.0
        for (rtime, bpsize) in self.anchor_pairs:
            rss += ( bpsize - f(rtime) ) ** 2

        # calculate the rest of peaks
        for (rtime, bpsize, rsize, err) in pairs:
            if rtime in self.anchor_rtimes:
                continue
            rss += err

        missing_peaks = len(self.sizes) - len(pairs) + 1

        if missing_peaks / len(self.sizes) > 0.5:
            # increase the penalty for missing ladders
            score = 1e3 *  missing_peaks ** 4
        else:
            score = rss * missing_peaks**self.penalty

        return score


def align_gm( peaks, ladder, anchor_pairs, z=None):

    cerr('I: generalized minimization method is running!')

    sizes = ladder['sizes']

    f = ZFunc( peaks, sizes, anchor_pairs)
    if z is None:
        zresult = f.get_initial_z()
        z = zresult.z

    # try pair f
    #result = pair_f( np.poly1d(z), f.rtimes, sizes)
    #import pprint; pprint.pprint(result)

    rss = -1
    prev_rss = 0
    #print('>>> Initial rss: ', rss)
    #plot(f.rtimes, f.sizes, z, result)

    #minimizer_kwargs = {'method': 'BFGS'}
    #bounds = [ (-1e-10,1e-10), (-1e-5,1e-5), (0,1e-3), (-1e2,1e2) ]
    bounds = [ (-1e-10, 1e-10), (-1e-5, 1e-5), (0.05, 0.18), (-175, 10) ]

    niter = 1
    results = []
    while abs(rss - prev_rss) > 1e-3:

        prev_rss = rss

        #res = minimize(f, z, method='Powell', tol=1e-6)
        #res = minimize(f, z, method='SLSQP', tol = 1e-6, bounds=bounds)
        res = minimize(f, z, method='Nelder-Mead', tol=1e-6)

        pairs, final_rss = f.get_pairs(res.x)

        rtimes, bpsizes = zip( *pairs)
        zresult = estimate_z(rtimes, bpsizes, niter if niter < 3 else 3)
        rss = zresult.rss
        z = zresult.z
        cerr('I: GM iter: %2d  - pairs: %2d  - Cur RSS: %6.2f' % (niter, len(pairs), rss))
        niter += 1
        results.append( zresult )

        #plot(f.rtimes, f.sizes, z, pairs)

        if rss < len(pairs) * 1.0:
            break

    # get the best result
    results.sort( key = lambda x: x.rss )
    zresult = results[0]

    # last dp
    dp_result = align_dp(f.rtimes, f.sizes, zresult.z, zresult.rss)
    #import pprint; pprint.pprint(dp_result.sized_peaks)
    #plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])

    dp_result.sized_peaks = f.get_sized_peaks(dp_result.sized_peaks)

    score, msg = ladder['qcfunc'](dp_result, method='strict')
    if score > 0.9:
        return AlignResult(score, msg, dp_result, const.alignmethod.gm_strict)

    score, msg = ladder['qcfunc'](dp_result, method='relax')
    return AlignResult(score, msg, dp_result, const.alignmethod.gm_relax)


def align_sh( peaks, ladder ):
    """ stochastic-heuristic method
        this is a combination of differential evolution method on a subset of data
        and then minimization for further processing
    """

    cerr('I: semi-heuristic method is running!')

    subset_peaks = [ p for p in peaks if 1500 < p.rtime < 5000 ]

    initial_pairs, initial_z = estimate_de( subset_peaks, ladder['signature'])

    result = align_gm( peaks, ladder, initial_pairs, initial_z )
    result.method = ( const.alignmethod.sh_strict
                        if result.method == const.alignmethod.gm_strict
                        else const.alignmethod.sh_relax )
    return result

def align_de( peaks, ladder, initial_pair=[] ):
    """ differential evolution method
    """

    cerr('I: differential evolution method is running!')

    sizes = ladder['sizes']

    f = ZFunc( peaks, sizes, initial_pair )
    bounds = [ (-1e-10, 1e-10), (-1e-5, 1e-5), (0.01, 0.1), (-75, 75) ]

    niter = 0
    results = []

    while niter < 3:

        #prev_rss = rss

        res = differential_evolution(f, bounds, tol=1e-5, mutation=(0.4, 1.5),
                popsize=30, recombination=0.8)

        pairs, final_rss = f.get_pairs(res.x)
        rtimes, bpsizes = zip( *pairs)
        zres = estimate_z(rtimes, bpsizes, niter if niter < 3 else 3)

        niter += 1
        cerr('I: DE iter: %2d  - pairs: %2d  - Cur RSS: %6.2f' % (niter, len(pairs), zres.rss))
        results.append( zres )

        if zres.rss < len(bpsizes) * 1.0:
            break

    results.sort( key = lambda x: x.rss )
    zres = results[0]

    # last dp
    dp_result = align_dp(f.rtimes, f.sizes, zres.z, zres.rss)
    #plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])
    #import pprint; pprint.pprint(dp_result.sized_peaks)

    dp_result.sized_peaks = f.get_sized_peaks(dp_result.sized_peaks)

    score, msg = ladder['qcfunc'](dp_result, method='relax')
    return AlignResult(score, msg, dp_result, const.alignmethod.de_relax)


def estimate_de( peaks, sizes ):

    f = ZFunc( peaks, sizes, [], estimate=True )
    bounds = [ (0.01, 0.5), (-275, 75) ]


    niter = 0
    results = []

    while niter < 3:

        #prev_rss = rss

        res = differential_evolution(f, bounds, tol=1e-5, mutation=(0.3, 1.7),
                popsize=45, recombination=0.5, strategy='rand1bin')

        pairs, final_rss = f.get_pairs(res.x)
        pairs.sort()
        rtimes, bpsizes = zip( *pairs)
        zres = estimate_z(rtimes, bpsizes, 1)

        niter += 1
        cerr('I: DE iter: %2d  - pairs: %2d  - Cur RSS: %6.2f' % (niter, len(pairs), zres.rss))
        results.append( (zres, pairs ) )

        if zres.rss < len(bpsizes) * 1.0:
            break

    results.sort( key = lambda x: x[0].rss )
    zres, pairs = results[0]

    plot(f.rtimes, f.sizes, zres.z, pairs)

    return pairs, zres.z
