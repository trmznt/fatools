
from fatools.lib.utils import cout, cerr
from fatools.lib.const import alignmethod
from fatools.lib.fautil import dpalign as dp

import numpy as np
import pprint


cdbg = cerr

def fast_align( data, ladders, peaks , qcfunc ):

    peaks = sorted( peaks, key = lambda x: x.rtime )

    hq_peaks = [ p for p in peaks if p.qscore >= 0.99 ]

    #cerr( ' -> HQ peak: %d / %d' % (len(hq_peaks), len(peaks)) )

    hq_score = mq_score = -1

    if len(hq_peaks) > 2 and ( len(hq_peaks) >= len(ladders) - 1
                                or len(peaks) == len(ladders) ):
        cerr(' -> fast_align() HQ peak: %d => ' % len(hq_peaks), nl=False)
        hq_result = do_fast_align( hq_peaks, ladders, peaks )
        (hq_score, hq_msg) = qcfunc( hq_result, method = 'strict' )
        cerr('%d | %4.2f | %6.2f | %3.2f' %
                (len(hq_result[3]), hq_result[0], hq_result[1], hq_score))
        if hq_score > 0.95:
            return (hq_score, hq_msg, hq_result, alignmethod.fast_hq)

    mq_peaks = [ p for p in peaks if p.qscore >= 0.75 ]

    if len(mq_peaks) > len(hq_peaks) and len(mq_peaks) >= len(ladders) - 1:
        cerr(' -> fast_align() MQ peak: %d => ' % len(mq_peaks), nl=False)
        mq_result = do_fast_align( mq_peaks, ladders, peaks )
        (mq_score, mq_msg) = qcfunc( mq_result, method = 'strict' )
        cerr('%d | %4.2f | %6.2f | %3.2f' %
                (len(mq_result[3]), mq_result[0], mq_result[1], mq_score))
        if mq_score > 0.95:
            return (mq_score, mq_msg, mq_result, alignmethod.fast_mq)

    hq_relax_score = mq_relax_score = -1
    if hq_score >= 0:
        hq_relax_score, hq_relax_msg = qcfunc(hq_result, method='relax')
    if mq_score >= 0:
        mq_relax_score, mq_relax_msg = qcfunc(mq_result, method='relax')

    if hq_relax_score < 0 and mq_relax_score < 0:
        return (0, None, None, None)

    if hq_relax_score >= mq_relax_score:
        return (hq_relax_score, hq_relax_msg, hq_result, alignmethod.fast_hqr)
    return (mq_relax_score, mq_relax_msg, mq_result, alignmethod.fast_mqr)



def shift_align( data, ladders, peaks, qcfunc):
    """ perform fast_align but with shifted peaks aka removing last peaks """

    for i in [1, 2, 3, 4, 5]:
        score = fast_align( data, ladders[:-i], peaks, qcfunc)
        if score[0] > 0.5:
            return score

    return (0, None, None, None)


def greedy_align( data, ladders, peaks, qcfunc ):

    peaks = sorted( peaks, key = lambda x: x.rtime )

    # just use the first and last peaks to determine Z

    delta = len(peaks) - len(ladders)
    if delta > 0:
        # high likely of having noises / overlapping peaks
        I, J = 4 + delta, 5 + delta
    elif delta < 0:
        # losing some ladder peaks, however we opt to always obtain max_no_peaks!
        I, J = 4, 5
        ladders = ladders[:delta]
    else:
        I, J = 4, 5

    initial_dp_results = []

    for i in range(0, I):
        for j in range(1, J):

            # use the first & last peak to estimate Z
            initial_peak_pairs = [  (peaks[i].rtime, ladders[1]),
                                    (peaks[-j].rtime, ladders[-1]) ]

            # use this for z_align, and obtain rss
            peak_pairs = z_align( peaks, ladders, initial_peak_pairs )
            z, rss = dp.estimate_z( * zip( *peak_pairs ) )

            initial_dp_results.append( (rss, peak_pairs, z, initial_peak_pairs) )

    # sort based on the smallest rss
    seed_peek_pairs = sorted(initial_dp_results, key = lambda x: (x[0], -len(x[1])))

    reports = []

    for (rss, peak_pairs, z, initial_peak_pairs) in seed_peek_pairs:

        # perform iterative alignment
        cerr('=> Iterative alignment')
        pprint.pprint( initial_peak_pairs )
        result = iterative_align( peaks, ladders, peak_pairs, peaks )
        (score, msg) = qcfunc(result[:4], 'relax')
        if score > 0.99:
            return (score, msg, result[:4], alignmethod.greedy_filtered)

        dp_score, dp_rss, dp_z, dp_peaks = result[:4]

        cdbg('Current best score: %3.2f dp: %3.2f rss: %4.2f' % (score, dp_score, dp_rss))
        cdbg('Start aligning by shifting pairs...')

        counter = 0
        while True:
            peak_pairs = shift_peak_pairs( dp_peaks, dp_z )
            result_shift = iterative_align( peaks, ladders, peak_pairs, peaks )
            score_shift, msg_shift = qcfunc( result_shift[:4], 'relax' )
            if score_shift > 0.99:
                return (score_shift, msg_shift, result_shift[:4], alignmethod.greedy_shifted)
            if result_shift[0] <= dp_score and result_shift[1] >= dp_rss and score_shift <= score:
                break
            if result_shift[0] <= dp_score and counter > 50:
                break
            dp_score, dp_rss, dp_z, dp_peaks = result_shift[:4]
            score = score_shift
            msg = msg_shift
            counter += 1

        reports.append( (score, msg, (dp_score, dp_rss, dp_z, dp_peaks),
                                            alignmethod.greedy_scored) )

    reports.sort( reverse = True, key = lambda x: x[0] )
    return reports[0]


def z_align( peaks, ladders, peak_pairs ):
    """ generated Z (3rd polynomial curve function) based on peak pairs,
        then realign ladders & peaks """

    # peak_pairs is [ (rtime, size), ... ]

    last_rss = -1
    step = 1

    while True:

        # optimize locally just by using simple peak assignment

        z, _rss_ = dp.estimate_z( [ x[1] for x in peak_pairs ],
                                [ x[0] for x in peak_pairs ],
                                3 if len(peak_pairs) > 2 else 1)

        peak_pairs, rss = estimate_peak_pairs( peaks, ladders, z )
        #cdbg(' z_align() step %d ==> %5.2f' % (step, rss))
        #pprint.pprint(peak_pairs)

        if last_rss < 0:
            last_rss = rss
        elif last_rss - rss <= 1:
            break

        last_rss = rss
        step += 1

    return peak_pairs


def do_fast_align( initial_peaks, ladders, all_peaks ):

    # simple fast aligner

    #if len(initial_peaks) < len(ladders):
    #    raise RuntimeError("shouldn't be here")
    peak_pairs = list(zip( reversed( [p.rtime for p in initial_peaks] ), reversed( ladders ) ))
    peak_pairs.sort()
    #print('fast_align with %d peaks.' % len(peak_pairs))
    #pprint.pprint(peak_pairs)
    if len(initial_peaks) < len(ladders):
        initial_peaks = all_peaks
    result = progressive_align(initial_peaks, ladders, peak_pairs)
    #print(' * dp_score =', dp_score )
    #pprint.pprint(dp_peaks)
    return result[:4]


def progressive_align( peaks, ladders, peak_pairs ):
    """ perform one-step dynamic programming alignment """

    peak_pairs = z_align( peaks, ladders, peak_pairs )
    z, rss = dp.estimate_z( * zip( *peak_pairs ) )
    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp.align_peaks( ladders, peaks, z, rss )
    cdbg(' progressive_align() ==> %4.2f | %5.2f | %d' % (dp_score, dp_rss, len(dp_peaks)))
    return (dp_score, dp_rss, dp_z, dp_peaks, S, D)


def iterative_align( peaks, ladders, peak_pairs, all_peaks ):
    """ peaks is initial peaks """


    dp_results = []

    last_dp = last_rss = -1
    last_z = last_peaks = last_S = last_D = None
    while True:

        dp_score, dp_rss, dp_z, dp_peaks, S, D = progressive_align(peaks, ladders, peak_pairs)
        #print(' * dp_score =', dp_score )

        if abs(dp_score - last_dp) < 0.1:
            #print(' ==> latest dp_score: ', dp_score )
            break
        elif dp_score < last_dp:
            dp_score = last_dp
            dp_rss = last_rss
            dp_z = last_z
            dp_peaks = last_peaks
            S = last_S
            D = last_D
            break

        #cerr(' ** ===> reiterate peak alignment' )
        last_dp = dp_score
        last_rss = dp_rss
        last_z = dp_z
        last_peaks = dp_peaks
        last_S = S
        last_D = D

        # recreate peak_pairs based on db_peaks
        #peak_pairs = [ (x[1].rtime, x[0]) for x in dp_peaks ]
        #if dp_rss > min_rss:
        peaks = all_peaks
        peak_pairs = reassign_peaks( peaks, ladders, dp_z )
        #pprint.pprint(peak_pairs)

    return( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )


def reassign_peaks( peaks, ladders, z ):
    """ reassign peaks & ladders, with (z)(rtime) -> bp size """

    #print(' => Z:', z)
    f = np.poly1d(z)
    peak_rtime = [ ( f(p.rtime), p.rtime ) for p in peaks ]
    peak_pairs = []
    for l in ladders:
        distances = [ ( abs(x[0] - l), x[1] ) for x in peak_rtime ]
        distances.sort()
        peak_pairs.append( ( distances[0][1], l ) )
    peak_pairs.sort()
    return peak_pairs


def shift_peak_pairs( peak_pairs, z ):
    """ shift peak pairs based on the highest deviation """

    f = np.poly1d(z)
    max_i = 0
    max_value = 0

    for (idx, peak_pair) in enumerate(peak_pairs):
        deviation = abs( f( peak_pair[1].rtime ) - peak_pair[0] )
        #cerr('  %2d  %3d  %6d  %5.2f' % (idx, peak_pair[0], peak_pair[1].rtime, deviation))
        if deviation > max_value:
            max_value = deviation
            max_i = idx

    new_peak_pairs = []
    if max_i > len(peak_pairs)/2:
        # on the right hand
        for i in range(0, max_i):
            #new_peak_pairs.append( peak_pairs[i] )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i][0]) )
        for i in range(max_i+1, len(peak_pairs)):
            #new_peak_pairs.append( (peak_pairs[i-1][0], peak_pairs[i][1]) )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i-1][0]) )

    else:
        for i in range(0, max_i):
            #new_peak_pairs.append( (peak_pairs[i+1][0], peak_pairs[i][1]) )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i+1][0]) )

        for i in range(max_i+1, len(peak_pairs)):
            #new_peak_pairs.append( peak_pairs[i] )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i][0]) )

    #cerr('after shifting...')
    #for (idx, peak_pair) in enumerate(new_peak_pairs):
    #    deviation = abs( f( peak_pair[0] ) - peak_pair[1] )
    #    #cerr('  %2d  %3d  %6d  %5.2f' % (idx, peak_pair[1], peak_pair[0], deviation))

    return new_peak_pairs


def estimate_peak_pairs( peaks, ladders, z ):
    """ estimate peak pairs from peaks & ladders with (z)(bp size) -> rtime
        return peak pairs [ (rtime, bpsize), ... ]
    """

    # generate synthetic standard peaks
    standard_peaks = []
    f = np.poly1d(z)
    for ladder in ladders:
        standard_peaks.append( ( f(ladder), ladder ) )

    #pprint.pprint( standard_peaks )

    # align peaks, must be done twice

    peak_pairs = []
    rss = 0

    for (rtime, ladder) in standard_peaks:
        min_dev = abs( rtime - peaks[0].rtime )
        aln_peak = peaks[0]
        for p in peaks[1:]:
            c = abs( rtime - p.rtime )
            if c < min_dev:
                min_dev = c
                aln_peak = p
        peak_pairs.append( (aln_peak.rtime, ladder) )
        rss += min_dev

    peak_pairs.sort()
    return (peak_pairs, rss)

