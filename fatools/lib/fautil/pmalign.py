"""
pair minimization algorithm
"""


import numpy as np
import itertools
from scipy.optimize import minimize

from fatools.lib.utils import cout, cerr, cverr, is_verbosity
from fatools.lib.fautil.alignutils import (estimate_z, pair_f, align_dp,
        pair_sized_peaks, DPResult, AlignResult, plot)
from fatools.lib.fautil.gmalign import ZFunc, align_gm
from fatools.lib import const

ANCHOR_RTIME_LOWER_BOUND = 1400
ANCHOR_RTIME_UPPER_BOUND = 5000
PEAK_RTIME_UPPER_BOUND = 11000


def align_pm(peaks, ladder, anchor_pairs=None):

    if not anchor_pairs:
        longest_rtime_peak = max([p.rtime for p in peaks])
        if longest_rtime_peak > PEAK_RTIME_UPPER_BOUND:
            bound_adjust_ratio = longest_rtime_peak / PEAK_RTIME_UPPER_BOUND
            anchor_start = ANCHOR_RTIME_LOWER_BOUND * bound_adjust_ratio
            anchor_end = ANCHOR_RTIME_UPPER_BOUND * bound_adjust_ratio
        else:
            anchor_start = ANCHOR_RTIME_LOWER_BOUND
            anchor_end = ANCHOR_RTIME_UPPER_BOUND
        anchor_peaks = [ p for p in peaks if anchor_start < p.rtime < anchor_end ]
        anchor_pairs, initial_z = estimate_pm( anchor_peaks, ladder['signature'] )

    else:
        rtimes, bpsizes = zip( *anchor_pairs )
        initial_z = estimate_z(rtimes, bpsizes, 1)

    anchor_pairs.sort()
    pairs, z, rss, f = align_upper_pm(peaks, ladder, anchor_pairs, initial_z)
    #print(pairs)
    pairs, z, rss, f = align_lower_pm(peaks, ladder, pairs, initial_z)

    #print(rss)
    #plot(f.rtimes, f.sizes, z, pairs)
    # last dp
    dp_result = align_dp(f.rtimes, f.sizes, f.similarity, z, rss)
    if is_verbosity(1):
        import pprint; pprint.pprint(dp_result.sized_peaks)
    if is_verbosity(4):
        plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])

    dp_result.sized_peaks = f.get_sized_peaks(dp_result.sized_peaks)

    score, msg = ladder['qcfunc'](dp_result, method='strict')
    if score > 0.9:
        return AlignResult(score, msg, dp_result, const.alignmethod.pm_strict)

    score, msg = ladder['qcfunc'](dp_result, method='relax')
    return AlignResult(score, msg, dp_result, const.alignmethod.pm_relax)


    f = ZFunc(peaks, ladder['sizes'], anchor_pairs)

    z = initial_z
    score = last_score = 0
    last_z = None

    for order in [1, 2, 3]:

        last_rss = -1
        rss = 0

        niter = 0
        while abs(rss - last_rss) > 1e-3:

            niter += 1
            cverr(5, 'Iter: %d' % niter)

            cverr(5, z)
            score = f(z)
            if last_score and last_score < score:
                # score does not converge; just exit
                cverr(5, 'does not converge!')
                break

            pairs, cur_rss = f.get_pairs(z)
            rtimes, bpsizes = zip( *pairs )
            zres = estimate_z(rtimes, bpsizes, order)

            last_z = z
            z = zres.z
            last_rss = rss
            rss = zres.rss
            cverr(5, rss)

    dp_result = align_dp(f.rtimes, f.sizes, last_z, last_rss)

    return align_gm2(peaks, ladder, anchor_pairs, dp_result.z)



    new_anchor_pairs = []
    zf = np.poly1d(dp_result.z)
    for p in dp_result.sized_peaks:
        if (p[0] - zf(p[1]))**2 < 2:
            new_anchor_pairs.append( (p[1], p[0]) )
    #import pprint; pprint.pprint(dp_result.sized_peaks)
    plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])

    return align_gm(peaks, ladder, anchor_pairs, dp_result.z)


def align_lower_pm(peaks, ladder, anchor_pairs, anchor_z):

    # anchor pairs must be in asceding order
    print(anchor_pairs)
    raise RuntimeError

    last_rtime = anchor_pairs[-1][0]
    last_size = anchor_pairs[-1][1]
    anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
    anchor_rtimes = list(anchor_rtimes)
    anchor_bpsizes = list(anchor_bpsizes)

    # preparing all "lower peaks", including anchor peaks
    lower_peaks= [ p for p in peaks if p.rtime <= last_rtime ]
    lower_sizes = [ s for s in ladder['sizes'] if s <= last_size]

    # we try to pair-minimize lower_peaks and lower_sizes

    f = ZFunc(peaks, ladder['sizes'], anchor_pairs, estimate=True)

    # check the first
    print(anchor_z)
    est_first_bpsize = anchor_z[2] * lower_peaks[0].rtime + anchor_z[3]
    print('est_first_bpsize:', est_first_bpsize)
    first_bpsize = [ s for s in lower_sizes if s >= est_first_bpsize ][0]

    scores = []
    for first_peak in lower_peaks[:-2]:
        if first_peak.rtime >= anchor_pairs[0][0]:
            break

        for first_bpsize in ladder['sizes'][:2]:
            zres = estimate_z( [ first_peak.rtime ] + anchor_rtimes, [ first_bpsize ] + anchor_bpsizes, 3 )
            #print('rss:', zres.rss)
            #plot(f.rtimes, f.sizes, zres.z, [ (first_peak.rtime, first_bpsize), ] )
            score, z = minimize_score(f, zres.z, 3)

            scores.append( (score, z) )
            #plot(f.rtimes, f.sizes, z, [ (first_peak.rtime, first_bpsize), ] )

    scores.sort( key = lambda x: x[0] )
    #import pprint; pprint.pprint( scores[:10] )

    z = scores[0][1]
    pairs, rss = f.get_pairs(z)

    print(rss)
    #plot(f.rtimes, f.sizes, z, pairs )

    return pairs, z, rss, f


    raise RuntimeError


def align_lower_pm(peaks, ladder, anchor_pairs, anchor_z):

    # this is another attempt to perform ladder - size standard alignment one peak by one


    while True:

        anchor_pairs = sorted(anchor_pairs)
        anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
        anchor_rtimes = list(anchor_rtimes)
        anchor_bpsizes = list(anchor_bpsizes)
        remaining_sizes = [x for x in ladder['sizes'] if x < anchor_bpsizes[0]]
        if not remaining_sizes:
            return pairs, z, rss, f

        current_sizes = [ remaining_sizes[-1] ] + anchor_bpsizes
        print('current_sizes:', current_sizes)
        f = ZFunc(peaks, current_sizes, anchor_pairs, estimate=True)
        zres = estimate_z(anchor_rtimes, anchor_bpsizes, 3)
        score, z = minimize_score(f, zres.z, 3)
        pairs, rss = f.get_pairs(z)
        plot(f.rtimes, f.sizes, z, pairs )
        anchor_pairs = pairs


def align_lower_pm(peaks, ladder, anchor_pairs, anchor_z):

    # this is another attempt to perform ladder - size standard alignment one peak by one

    anchor_pairs = sorted(anchor_pairs)
    anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
    anchor_rtimes = list(anchor_rtimes)
    anchor_bpsizes = list(anchor_bpsizes)
    remaining_sizes = [x for x in ladder['sizes'] if x < anchor_bpsizes[0]]
    current_sizes = anchor_bpsizes
    zscore = estimate_z(anchor_rtimes, anchor_bpsizes, 3)
    z = zscore.z
    rss = zscore.rss
    f = ZFunc(peaks, current_sizes, anchor_pairs)

    while True:

        if not remaining_sizes:
            return pairs, z, rss, f

        current_sizes.insert(0, remaining_sizes.pop(-1))
        f.set_sizes(current_sizes)
        score, next_z = minimize_score(f, z, 3)
        next_pairs, next_rss = f.get_pairs(next_z)

        # if delta rss (current rss - prev rss) is above certain threshold,
        # then assume the latest peak standar is not appropriate, and
        # use previous z and rss
        if (next_rss - rss) > 20:
            current_sizes.pop(0)
        else:
            z = next_z
            rss = next_rss
            pairs = next_pairs

        if is_verbosity(5):
            plot(f.rtimes, f.sizes, z, pairs )


def align_upper_pm(peaks, ladder, anchor_pairs, anchor_z):

    # anchor pairs must be in asceding order


    anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
    anchor_rtimes = list(anchor_rtimes)
    anchor_bpsizes = list(anchor_bpsizes)

    # we try to pair-minimize higher peaks and sizes

    first_rtime = anchor_rtimes[0]
    first_bpsize = anchor_bpsizes[0]
    peaks = [ p for p in peaks if p.rtime >= first_rtime]
    sizes = [ s for s in ladder['sizes'] if s >= first_bpsize]
    remaining_sizes = [ s for s in ladder['sizes'] if s > anchor_bpsizes[-1] ]

    #sizes = ladder['sizes']
    f = ZFunc(peaks, sizes, anchor_pairs, estimate=True)

    # check the first
    #print(peaks[-1])
    est_last_bpsize = np.poly1d(anchor_z)(peaks[-1].rtime)
    #est_last_bpsize = anchor_z[1] * peaks[-1].rtime**2 + anchor_z[2] * peaks[-1].rtime + anchor_z[3]
    #print(est_last_bpsize)
    last_bpsize = max( remaining_sizes[1], [ s for s in sizes if s < est_last_bpsize ][-3] )
    #print('last_bpsize:', last_bpsize)
    #plot(f.rtimes, f.sizes, anchor_z, [])

    scores = []
    #print(peaks)
    for last_peak in reversed(peaks[-14:]):
        if last_peak.rtime <= anchor_pairs[-1][0]:
            break

        zres = estimate_z(anchor_rtimes + [last_peak.rtime], anchor_bpsizes + [last_bpsize], 2)
        #plot(f.rtimes, f.sizes, zres.z, [ (last_peak.rtime, last_bpsize)] )
        score, z = minimize_score(f, zres.z, 2)
        #print(score)
        #plot(f.rtimes, f.sizes, z, [] )

        scores.append( (score, z) )

    scores.sort( key = lambda x: x[0] )
    #import pprint; pprint.pprint( scores[:10] )

    z = scores[0][1]
    pairs, rss = f.get_pairs(z)

    #print(rss)
    #plot(f.rtimes, f.sizes, z, pairs )

    return pairs, z, rss, f


def align_upper_pm(peaks, ladder, anchor_pairs, anchor_z):

    # this is another attempt to perform ladder - size standard alignment one peak by one

    anchor_pairs = sorted(anchor_pairs)
    anchor_rtimes, anchor_bpsizes = zip( *anchor_pairs )
    anchor_rtimes = list(anchor_rtimes)
    anchor_bpsizes = list(anchor_bpsizes)
    remaining_sizes = [x for x in ladder['sizes'] if x > anchor_bpsizes[-1]]
    current_sizes = anchor_bpsizes
    order = ladder['order']
    zres = estimate_z(anchor_rtimes, anchor_bpsizes, order)
    z,rss = zres.z, zres.rss
    f = ZFunc(peaks, current_sizes, anchor_pairs)

    while remaining_sizes:

        current_sizes.append( remaining_sizes.pop(0) )
        if ( remaining_sizes and
             (remaining_sizes[-1] - current_sizes[-1]) < 100 and
             (remaining_sizes[0] - current_sizes[-1]) < 11 ):
            current_sizes.append( remaining_sizes.pop(0) )

        f.set_sizes(current_sizes)
        score, next_z = minimize_score(f, z, order)
        next_pairs, next_rss = f.get_pairs(z)

        if (next_rss - rss) < 70:
            z = next_z
            rss = next_rss
            pairs = next_pairs

        if is_verbosity(5):
            plot(f.rtimes, f.sizes, z, pairs )

    # finalize the alignment with stringent criteria
    dp_result = align_dp(f.rtimes, f.sizes, f.similarity, z, rss)
    if dp_result.rss - rss > 50:
        return pairs, z, rss, f
    dp_pairs = [(x[1], x[0]) for x in dp_result.sized_peaks]
    if is_verbosity(5):
        plot(f.rtimes, f.sizes, dp_result.z, dp_pairs)

    return dp_pairs, dp_result.z, dp_result.rss, f


def minimize_score( f, z, order ):

    last_score = score = 0

    niter = 1
    while niter  < 50:

        score = f(z)
        #print(score)

        if last_score and abs(last_score - score) < 1e-6:
            break

        pairs, rss = f.get_pairs(z)
        rtimes, bpsizes = zip( *pairs )
        zres = estimate_z(rtimes, bpsizes, order)

        z = zres.z
        last_score = score
        niter += 1

    return last_score, z



def estimate_pm(peaks, bpsizes):

    rtimes = [ p.rtime for p in peaks ]

    rtime_points = prepare_rtimes( rtimes )
    bpsize_pair = [ bpsizes[1], bpsizes[-2]]

    f = ZFunc(peaks, bpsizes, [], estimate = True)

    scores = []
    for rtime_pair in rtime_points:
        if rtime_pair[0] >= rtime_pair[1]:
            continue

        # y = ax + b
        # y1 = ax1 + b
        # y2 = ax2 + b
        # ------------ -
        # y1 - y2 = a(x1 - x2)
        # a = (y1 - y2)/(x1 - x2)
        # b = y1 - ax1

        #slope = (bpsize_pair[1] - bpsize_pair[0]) / (rtime_pair[1] - rtime_pair[0])
        #intercept = bpsize_pair[0] - slope * rtime_pair[0]
        #z = [ slope intercept ]
        zres = estimate_z(rtime_pair, bpsize_pair, 1)
        score = f(zres.z)
        scores.append( (score, zres) )
        if is_verbosity(5):
            plot(f.rtimes, f.sizes, zres.z, [] )

    scores.sort( key = lambda x: x[0] )
    #import pprint; pprint.pprint(scores[:5])
    zresult = scores[0][1]

    dp_result = align_dp(f.rtimes, f.sizes, f.similarity, zresult.z, zresult.rss)
    #import pprint; pprint.pprint(dp_result.sized_peaks)
    if is_verbosity(5):
        plot(f.rtimes, f.sizes, dp_result.z, [(x[1], x[0]) for x in dp_result.sized_peaks])

    return ( [(x[1], x[0]) for x in dp_result.sized_peaks], dp_result.z )



def prepare_rtimes(rtimes):
    # prepare combination of begin and end rtimes

    mid_size = round(len(rtimes)/2)
    return list( itertools.product(rtimes[:mid_size], rtimes[mid_size-2:]) )
