
import numpy as np
import math

from scipy.signal import find_peaks_cwt
from scipy.optimize import leastsq, curve_fit
from scipy.interpolate import UnivariateSpline

from fatools.lib.const import peaktype, binningmethod, allelemethod
from fatools.lib.fautil.dpalign import estimate_z, align_peaks, plot_z
from fatools.lib.utils import cout, cerr, cverr

import fatools.lib.fautil.peakalign as pa

from sortedcontainers import SortedListWithKey

import pprint, sys, pickle
from matplotlib import pylab as plt
from bisect import bisect_left


def find_raw_peaks( raw_data, params ):

    max_height = max(raw_data)
    width_ratio = max(1, round(math.log( max_height/params.width_ratio )))
    widths = params.widths

    if params.method == 'cwt':
        from scipy.signal import find_peaks_cwt
        indices = find_peaks_cwt( raw_data, widths,
                                    min_snr = params.min_snr )
        #cerr('find_peaks_cwt() found %d peaks' % len(indices))
        #pprint.pprint(indices)

    elif params.method == 'relmax':
        indice_set = []
        from scipy.signal import argrelmax
        for i in (params.widths * width_ratio):
            indice_set.append( argrelmax( raw_data, order=i+5 )[0] )
        # get consensus
        indices = filter_by_snr( get_consensus_indices( indice_set ),
                                raw_data, params.min_snr * 3.5 )
        #print('indices => %d' % len(indices))
        #pprint.pprint( indices )

    elif params.method == 'mlpy':
        indice_set = []
        from mlpy import findpeaks_win
        for i in params.widths:
            indice_set.append( findpeaks_win( raw_data, span=i ) )
        # get consensus
        indices = filter_by_snr( get_consensus_indices( indice_set ),
            raw_data, params.min_snr )

    elif params.method == 'pd':
        from peakutils import indexes
        indices = indexes( raw_data, 1e-5, 10 )
        #pprint.pprint(indices)
        cverr(3, 'indice size: %d' % len(indices))
        cverr(3, 'indices => %s' % repr(indices))

    else:
        raise RuntimeError('unknown peak finding method: %s' % params.method)

    if indices is None or len(indices) == 0:
        return []


    # filter for absolute heights within proximity

    # special cases for pd (peak detect) method:

    if params.method == 'pd':
        return [ ( int(i), int(raw_data[i]) )
                    for i in indices if raw_data[i] > params.min_height and
                        params.min_rtime < i < params.max_rtime ]

    raw_peaks = []
    max_len = len(raw_data)
    for idx in indices:

        if not (params.min_rtime < idx < params.max_rtime):
            continue

        height, index = max( [ (raw_data[i], i)
                                for i in range(max(0, idx-3), min(max_len,idx+3) ) ] )

        if height < params.min_height: continue
        if (index, height) in raw_peaks: continue
        raw_peaks.append( (index, height) )

    return raw_peaks


def find_peaks( raw_data,  params, raw_peaks = None ):
    """
    find all peaks based on the criteria defined in params, and assign as peak-scanned
    raw_data is baseline-normalized & smoothed trace

    parameters used are:
    method: 'cwt' or 'mlpy'
    widths: window size for peak scanning
    cwt_min_snr:
    min_height:
    min_relative_ratio:
    max_relative_ratio:
    min_height_ratio:
    max_peak_number:

    """

    if raw_peaks is None:
        raw_peaks = find_raw_peaks( raw_data, params )

    # check for any peaks
    if not raw_peaks:
        return raw_peaks

    # only retain 2 * max_peak_number and discard the rest
    raw_peaks = sorted( raw_peaks, key = lambda x: x[1],
            reverse = True )[:params.max_peak_number * 2]

    if params.min_relative_ratio > 0 or params.max_relative_ratio > 0:
        med = np.median( list(p[1] for p in raw_peaks) )
        if params.min_relative_ratio > 0:
            median_min = med * params.min_relative_ratio
            raw_peaks = [ p for p in raw_peaks if p[1] > median_min ]
        if params.max_relative_ratio > 0:
            median_max = med * params.max_relative_ratio
            raw_peaks = [ p for p in raw_peaks if p[1] < median_max ]

    if not raw_peaks:
        return raw_peaks

    # filter for minimum height ratio

    if params.min_height_ratio > 0:
        min_height = max( list( p[1] for p in raw_peaks) ) * params.min_height_ratio
        raw_peaks = [ p for p in raw_peaks if p[1] > min_height ]

    # calculate area

    (q50, q75) = np.percentile( raw_data, [ 50, 75 ] )
    peaks = []
    for (peak, height) in raw_peaks:
        area, brtime, ertime, srtime, ls, rs = calculate_area( raw_data, peak, 5e-2, q50 )
        wrtime = ertime - brtime
        if wrtime < 3:
            continue
        beta = area / height
        theta = height / wrtime
        if height >= 25 and beta * theta < 6: #10:
            continue
        if height < 25 and beta * theta < 3: #6:
            continue
        peaks.append( (peak, height, area, brtime, ertime, srtime, beta, theta) )

    peaks.sort()
    cverr(3, 'peaks stage 1 size: %d' % len(peaks))
    cverr(3, 'peaks stage 1: %s' % repr(peaks))

    non_artifact_peaks = []

    for idx in range(len(peaks)):
        peak = peaks[idx]

        if idx > 0:
            prev_p = peaks[idx-1]
            if peak[3] - prev_p[4] < 5 and peak[1] < params.artifact_ratio * prev_p[1]:
                # we are artifact, just skip
                continue
        if idx < len(peaks)-1:
            next_p = peaks[idx+1]
            if next_p[3] - peak[4] < 5 and peak[1] < params.artifact_ratio * next_p[1]:
                # we are another artifact, just skip
                continue

        non_artifact_peaks.append( peak )

    cverr(3, 'max_peak_number: %d' % params.max_peak_number)

    sorted_peaks = sorted( non_artifact_peaks, key = lambda x: (x[1], x[6] * x[7]),
                        reverse=True )[:params.max_peak_number]
    peaks = sorted( sorted_peaks )
    cverr(3, 'peaks stage 3 size: %d' % len(peaks))
    cverr(3, 'peaks stage 3: %s' % repr(peaks))

    return peaks


##
## the methods below are operated with either channel data type or list of channels
##
## channel.data -> smoothed, normalized trace
## channel.new_allele() -> function to create & register new allele
## channel.get_alleles() -> function to get all alleles
##
## allele datatype


def scan_peaks( channel, params, peakdb ):
    """
    scan for peaks based on the criteria defined in params, set as peak-scanned,
    and prepare the channel data structure
    """

    if peakdb:
        raw_peaks = pickle.loads(peakdb.Get(channel.tag().encode()))
    else:
        raw_peaks = None

    initial_peaks = find_peaks( channel.data, params, raw_peaks)
    # peaks = ( rtime, height, area, brtime, ertime )
    #cerr('DEBUG - initial peaks: %d' % len(initial_peaks))

    cverr(3, 'initial peaks: %d' % len(initial_peaks))

    # perform futher cleaning for ladder channels
    if params.expected_peak_number:
        epn = params.expected_peak_number
        peak_qualities = sorted([ (p[6] * p[7], p) for p in initial_peaks ], reverse=True)
        low_scores = [ q[0] for q in peak_qualities[round(epn/3):round(epn * 1.5)] ]
        avg_low_score = sum(low_scores) / len(low_scores)
        ratio_low_score = (avg_low_score - low_scores[-1]) / low_scores[-1]
        if avg_low_score < 75:
            # questionable quality, please use more peaks
            score_threshold = 4 #avg_low_score * 0.1
            height_threshold = 6
        else:
            if avg_low_score - low_scores[-1] > low_scores[-1]:
            # peaks' height are likely not to evenly distributed
                score_threshold = max(low_scores[-1] * 0.90, 4)
            else:
                score_threshold = avg_low_score * 0.25
            height_threshold = 10
            cverr(3, 'using score threshold: %f' % score_threshold)
            cverr(3, 'using height_threshold: %d' % height_threshold)
        peaks = [ q for q in peak_qualities
                            if q[0] > score_threshold and q[1][1] > height_threshold ]
        cverr(3, 'after peak quality filtering: %d' % len(peaks))
        if len(peaks) > 1.5 * params.expected_peak_number:
            # try to remove peaks further
            saved_peaks = peaks
            while len(peaks) - len(saved_peaks) < 0.30 * len(peaks) and height_threshold < 20:
                height_threshold += 1
                saved_peaks = [ q for q in saved_peaks if q[0] > height_threshold ]
            peaks = saved_peaks
            cverr(3, 'after reducing peaks number by height: %d' % len(peaks))
        peaks = sorted( [ q[1] for q in peaks ] )

    else:
        peaks = initial_peaks


    # create alleles based on these peaks
    alleles = []
    for peak in peaks:
        ( rtime, height, area, brtime, ertime, srtime, beta, theta ) = peak
        wrtime = ertime - brtime
        height = round(height)
        allele = channel.new_allele(    rtime = rtime,
                                        height = height,
                                        area = area,
                                        brtime = brtime,
                                        ertime = ertime,
                                        wrtime = wrtime,
                                        srtime = srtime,
                                        beta = beta,
                                        theta = theta )
        allele.type = peaktype.scanned
        allele.method = binningmethod.notavailable
        allele.marker = channel.marker
        alleles.append( allele )

    return alleles



def preannotate_channels( channels, params ):
    """
    pre-annotate peaks as peak-scanned / peak-broad / peak-stutter/ peak-overlap
    based on criteria defined in params
    """

    # peak_1 is overlap of peak_2 if
    #   brtime2 < ertime1 and ertime1 < ertime2
    #   and height1 at rtime1 is a_fraction of height2 at rtime1 and
    #   height1 at rtime2 is a_fraction of height2 at rtime2.

    # peak is broad if beta > beta_broad_threshold

    channel_peaks = [ (list(channel.alleles), np.median(channel.data)) for channel in channels ]

    # reset all peak type, score the peaks and set the peak type to peak-noise,
    # peak-broad

    # collect beta * theta first, and used beta * theta as descriptor for noise
    # also if height at either brtime or ertime is higher than 50% at rtime, it is
    # likely a noise

    for (peaks, med_baseline) in channel_peaks:

        if len(peaks) == 0: continue

        beta_theta = sorted([ p.beta * p.theta for p in peaks ])
        sampled_beta_theta = beta_theta[2:len(beta_theta)-2]
        if len(sampled_beta_theta) == 0: sampled_beta_theta = beta_theta
        avg_beta_theta = sum(sampled_beta_theta) / len(sampled_beta_theta)


        for p in peaks:
            p.type = peaktype.scanned
            p.size = -1
            p.bin = -1

            peak_beta_theta = p.beta * p.theta
            score = 1.0

            # extreme noise

            if p.height < 2 * med_baseline:
                p.qscore = 0.25
                p.type = peaktype.noise
                continue

            if p.wrtime < 6 or (p.wrtime < 10 and peak_beta_theta < 0.275 * avg_beta_theta):
                p.qscore = 0.25
                p.type = peaktype.noise
                continue

            # moderately noise

            if peak_beta_theta < 0.33 * avg_beta_theta:
                if (    p.channel.data[p.brtime] > 0.5 * p.height or
                        p.channel.data[p.ertime] > 0.5 * p.height ):
                    p.qscore = 0.25
                    p.type = peaktype.noise
                    continue
                score -= 0.15

            score = 1.0
            if p.beta > params.max_beta:
                p.type = peaktype.broad
                score -= 0.20
            elif p.beta < 5:
                score -= 0.20

            # check theta
            if p.theta < 4:
                # perhaps an artifact
                score -= 0.20

            # penalty by height
            if p.height < 75:
                # decrease the score
                score -= 0.1
            if p.height < 50:
                # decrease even further
                score -= 0.1

            # penalty by symmetrics
            if not ( -1.32 < p.srtime < 1.32 ):
                score -= 0.1

            p.qscore = score
            if p.qscore < 0.5 and p.type == peaktype.scanned:
                p.type = peaktype.noise

            if p.qscore < 0:
                p.qscore = 0.0  # reset to zero

    # checking overlaps against channel !

    for channel in channels:
        for channel_r in channels:
            if channel == channel_r:
                continue

            for p in channel.alleles:
                if p.type == peaktype.noise:
                    continue

                if p.ertime - p.brtime < 3:
                    brtime = p.brtime
                    ertime = p.ertime
                elif p.ertime - p.brtime < 6:
                    brtime = p.brtime + 1
                    ertime = p.ertime - 1
                else:
                    brtime = p.brtime + 3
                    ertime = p.ertime - 3

                if brtime > p.rtime: brtime = p.rtime
                if ertime < p.rtime: ertime = p.rtime

                brtime = max(0, brtime)
                ertime = min(len(channel.data), len(channel_r.data), ertime)

                #cerr('checking %d | %s with channel %s' % (p.rtime, channel.dye,
                #            channel_r.dye))

                if (    channel.data[brtime] < channel_r.data[brtime] and
                        channel.data[ertime] < channel_r.data[ertime] and
                        p.height < channel_r.data[p.rtime] ):

                    # check how much is the relative height of this particular peak
                    rel_height = p.height / channel_r.data[p.rtime]
                    if rel_height > 1.0:
                        continue

                    (o_state, o_ratio, o_sym) = calc_overlap_ratio( channel.data,
                                                        channel_r.data, p.rtime,
                                                        brtime, ertime )

                    # if not really overlap, just continue reiteration
                    if not o_state:
                        continue

                    print('peak: %d | %s | %s <> %f | %f | %f' % (p.rtime, p.channel.dye, p.type, rel_height, o_ratio, o_sym))
                    if rel_height < 0.15:
                        if p.type != peaktype.noise:
                            p.type = peaktype.overlap
                            print('peak: %d | %s -> overlap' % (p.rtime, p.channel.dye))
                        p.qscore -= 0.10
                        continue

                    if ((rel_height < params.overlap_height_threshold and -0.5 < o_sym < 0.5) or
                        (o_ratio < 0.25 and -1.5 < o_sym < 1.5 ) or
                        (o_ratio < 0.75 and -0.5 < o_sym < 0.5 )):
                        if p.type != peaktype.noise:
                            p.type = peaktype.overlap
                            print('peak: %d | %s -> overlap' % (p.rtime, p.channel.dye))
                        p.qscore -= 0.10
                        continue

    # checking for stutter peaks based on minimum rtime & rfu

    for (peaks, med_baseline) in channel_peaks:
        alleles = sorted( [ p for p in peaks ],
                        key = lambda x: x.rtime )

        for idx in range( len(alleles) ):
            allele = alleles[idx]
            if idx > 0:
                allele_0 = alleles[idx-1]
                if allele.rtime - allele_0.rtime < params.stutter_rtime_threshold:
                    if allele_0.height * params.stutter_height_threshold > allele.height:
                        allele.type = peaktype.stutter
                        allele.qscore -= 0.2
            if idx < len(alleles) - 1:
                allele_1 = alleles[idx+1]
                if allele_1.rtime - allele.rtime < params.stutter_rtime_threshold:
                    if allele_1.height * params.stutter_height_threshold > allele.height:
                        allele.type = peaktype.stutter
                        allele.qscore -= 0.2



def size_peaks( channel, params, ladders, qcfunc = None ):

    data = channel.data
    scores = []

    # perform fast_align with both clean, high quality peaks and good peaks
    (score_0, msg_0, result_0, method_0) = pa.fast_align( data, ladders,
                                                channel.alleles, qcfunc )


    if score_0 > 0.99:
        return (score_0, msg_0, result_0, method_0)
    cerr('fast_align(): %4.2f' % score_0)
    scores.append( (score_0, msg_0, result_0, method_0) )

    # perform shift_align with both clean, high quality peaks and good peaks
    (score_0, msg_0, result_0, method_0) = pa.shift_align( data, ladders,
                                                channel.alleles, qcfunc )


    if score_0 > 0.99:
        return (score_0, msg_0, result_0, method_0)
    cerr('shift_align(): %4.2f' % score_0)
    scores.append( (score_0, msg_0, result_0, method_0) )

    # perform greedy alignment
    (score_1, msg_1, result_1, method_1) = pa.greedy_align( data, ladders,
                                                channel.alleles, qcfunc )

    scores.append( (score_1, msg_1, result_1, method_1) )

    scores.sort( key = lambda x: x[0] )
    return scores[-1]



def call_peaks( channel, params, func, min_rtime, max_rtime ):
    """
    call (determine size) each of peaks with type peak-scanned, and annotate as either
    peak-called or peak-unassigned
    """

    for allele in channel.alleles:
        if not min_rtime < allele.rtime < max_rtime:
            if allele.type == peaktype.scanned:
                allele.type = peaktype.unassigned
            continue
        size, deviation, qcall, method = func(allele.rtime)
        allele.size = size
        allele.bin = round(size)
        allele.deviation = deviation
        allele.qcall = qcall
        if allele.type == peaktype.scanned:
            allele.type = peaktype.called
        allele.method = binningmethod.notavailable



def bin_peaks(channel, params, marker):

    #sortedbins = marker.sortedbins
    sortedbins = marker.get_sortedbins(channel.batch)
    threshold = float(marker.repeats)/2 * 1.5

    for peak in channel.alleles:

        if peak.size < 0: continue

        if not marker.min_size < peak.size < marker.max_size:
            peak.type = peaktype.unassigned
            continue

        size = peak.size
        idx = sortedbins.bisect_key_right( size )

        if idx==0:
            curr_bin = sortedbins[0]
        elif idx == len(sortedbins):
            curr_bin = sortedbins[-1]
        else:
            left_bin = sortedbins[idx-1]
            right_bin = sortedbins[idx]

            if size - left_bin[3] < right_bin[2] - size:
                # belongs tp left_bin
                curr_bin = left_bin
            else:
                curr_bin = right_bin

        peak.bin = curr_bin[0]

        # only assigned peak as bin if it unassigned or called
        if peak.type in [ peaktype.unassigned, peaktype.called ]:
            peak.type = peaktype.bin


def postannotate_peaks( channel, params ):
    """
    post annotate binned peaks with peak stutter and broad signals
    """

    # peak is stutter if the range < params.stutter_range and
    # ratio < params.stutter_ratio

    alleles = sorted(list(channel.alleles), key = lambda x: x.height, reverse=True)
    prev_alleles = []

    for allele in alleles:
        if allele.type != peaktype.bin: continue
        if allele.beta > params.max_beta:
            allele.type = peaktype.broad
            continue
        for prev_allele in prev_alleles:
            if (    abs(prev_allele.size - allele.size) < params.stutter_range and
                    allele.height/prev_allele.height < params.stutter_ratio ):
                allele.type = peaktype.stutter
            elif ( abs(prev_allele.size - allele.size) < (params.stutter_range/2+1)):
                allele.type = peaktype.stutter
            elif prev_allele.bin == allele.bin:
                allele.type = peaktype.stutter

        # added this allele to prev_alleles
        prev_alleles.append( allele )


    # XXX: Note on stutter recognition
    # real stutters can be defined as a peak having at least 2 minor peak in consecutive
    # or around it, minor peak < stutter_ratio


# helper functions

def get_consensus_indices( indices_set ):

    n = len(indices_set)
    indices = {}
    for indices_item in indices_set:
        for i in indices_item:
            try:
                indices[i] += 1
            except KeyError:
                indices[i] = 1
    threshold = n/2
    real_indices = []
    for (k,v) in indices.items():
        if v > threshold:
            real_indices.append( k )

    return sorted( real_indices )


def filter_by_snr( indices, raw_data, snr ):

    size = len(raw_data)
    background_height = np.mean( np.sort(raw_data)[ size/4 : size - size/4 ] )
    min_height = background_height * snr
    return [ i for i in indices if raw_data[i] > min_height ]


def calculate_area(y, t, threshold, baseline):
    """ return (area, brtime, ertime, srtime)
        area: area
        brtime: begin rtime
        ertime: end rtime
    """

    # right area
    data = y[t:]
    r_area, ertime, r_shared = half_area(data, threshold, baseline)

    # left area
    data = y[:t+1][::-1]
    l_area, brtime, l_shared = half_area(data, threshold, baseline)


    return ( l_area + r_area - y[t], t - brtime, ertime + t, math.log2(r_area / l_area),
                l_shared, r_shared )


def half_area(y, threshold, baseline):
    """ return (area, ertime, shared_status)
    """

    winsize = 3
    threshold = threshold/2
    shared = False
    area = y[0]
    edge = float(np.sum(y[0:winsize]))/winsize
    old_edge = 2 * edge

    index = 1
    limit = len(y)

    while ( edge > area * threshold and edge < old_edge and
            index < limit and y[index] >= baseline ):
        old_edge = edge
        area += y[index]
        edge = float(np.sum(y[index:index+winsize]))/winsize
        index += 1
    if edge >= old_edge:
        shared = True
    index -= 1

    return area, index, shared


def is_overlap(peak_1, peak_2):

    if peak_1.brtime > peak_2.brtime:
        peak_1, peak_2 = peak_2, peak_1

    if peak_2.brtime < peak_1.ertime:
        return True

    return False


def is_definitive_overlap(peak_1, data_1, peak_2, data_2):
    """ is peak_1 is overlapped of peak_2 """
    for x in range(peak_1.brtime + 3, peak_1.ertime - 3):
        print('x: %d 1[x]: %d 2[x]: %d' % (x, data_1[x], data_2[x]))
        if data_1[x] > data_2[x]:
            return False
    return True


def calc_overlap_ratio(data, data_r, rtime, brtime, ertime):
    """ calculate the difference ratio of overlapping area to the left and right area
        of rtime
        return ( boolean, total_ratio, log(right_ratio/left_ratio) )
    """

    lr = rr = 0.0
    lc = rc = 0
    for x in range(brtime, rtime+1):
        if data[x] > data_r[x]:
            return (False, 0, 0)
        lr += data[x]/data_r[x]
        lc += 1
    for x in range(rtime, ertime+1):
        if data[x] > data_r[x]:
            return (False, 0, 0)
        rr += data[x]/data_r[x]
        rc += 1
    lrc = lr / lc
    rrc = rr / rc

    return (True, (lrc + rrc)/2, math.log2(rrc/lrc))


def generate_scoring_function( strict_params, relax_params ):

    def _scoring_func( alignment_result, method ):
        # alignment_result is (dp_score, dp_rss, dp_z, dp_peaks)
        dp_score, dp_rss, dp_z, dp_peaks = alignment_result

        if method == 'strict':
            if ( dp_score >= strict_params['min_dpscore'] and
                    dp_rss <= strict_params['max_rss'] and
                    len(dp_peaks) >= strict_params['min_sizes'] ):
                return (1, None)
            return (0, None)
        elif method == 'relax':
            msg = []
            # scoring based on parts of results

            # score based on DP score compared to minimum DP score
            delta_score = relax_params['min_dpscore'] - dp_score
            if delta_score <= 0:
                dp_score_part = 1
            else:
                dp_score_part = 1e-2 ** (1e-2 * delta_score)

            # score based on RSS compared to the maximum allowed RSS
            delta_rss = dp_rss - relax_params['max_rss']
            if delta_rss <= 0:
                dp_rss_part = 1
            else:
                dp_rss_part = 1e-2 ** ( 1e-3 * delta_rss )
                msg.append( 'RSS > %d' % ( relax_params['max_rss'] ) )

            # score based on how many peaks we might miss compared to minimum number of peaks
            delta_peaks = relax_params['min_sizes'] - len(dp_peaks)
            if delta_peaks <= 0:
                dp_peaks_part = 1
            else:
                dp_peaks_part = max( 0, - delta_peaks / 0.5 * relax_params['min_sizes'] - 1)
                msg.append( 'Missing peaks = %d' % delta_peaks )

            # total overall score
            score = 0.3 * dp_score_part + 0.5 * dp_rss_part + 0.2 * dp_peaks_part
            return (score, msg)

        raise RuntimeError("Shouldn't be here!")


    return _scoring_func


def least_square( ladder_alleles, z ):
    """ 3rd order polynomial resolver
    """
    ladder_allele_sorted = SortedListWithKey( ladder_alleles, key = lambda k: k.rtime )
    f = np.poly1d(z)

    def _f( rtime ):
        size = f(rtime)
        # get the left-closest and right-closest ladder

        #left_idx = ladder_allele_sorted.bisect_key_left( rtime )
        right_idx = ladder_allele_sorted.bisect_key_right( rtime )
        left_idx = right_idx - 1
        left_ladder = ladder_allele_sorted[left_idx]
        right_ladder = ladder_allele_sorted[right_idx]
        #cerr(' ==> rtime: %d/%4.2f  [ %d/%4.2f | %d/%4.2f ]' % ( rtime, size,
        #            left_ladder.rtime, left_ladder.size,
        #            right_ladder.rtime, right_ladder.size))

        return (size, (left_ladder.deviation + right_ladder.deviation) / 2,
                        min( left_ladder.qscore, right_ladder.qscore ),
                        allelemethod.leastsquare)

    return _f


def cubic_spline( ladder_alleles ):
    """ cubic spline interpolation
        x is peaks, y is standard size
    """

    ladder_allele_sorted = SortedListWithKey( ladder_alleles, key = lambda k: k.rtime )

    ladder_peaks = []
    ladder_sizes = []
    for ladder_allele in ladder_allele_sorted:
        ladder_peaks.append( ladder_allele.rtime )
        ladder_sizes.append( ladder_allele.size )
    f = UnivariateSpline(ladder_peaks, ladder_sizes, k=3, s=0)

    def _f( rtime ):
        size = f(rtime)

        right_idx = ladder_allele_sorted.bisect_key_right( rtime )
        left_idx = right_idx - 1
        left_ladder = ladder_allele_sorted[left_idx]
        right_ladder = ladder_allele_sorted[right_idx]

        return (size, (left_ladder.deviation + right_ladder.deviation) / 2,
                        min( left_ladder.qscore, right_ladder.qscore),
                        allelemethod.cubicspline)

    return _f



def local_southern( ladder_alleles ):
    """ southern local interpolation """

    ladder_allele_sorted = SortedListWithKey( ladder_alleles, key = lambda k: k.rtime )
    x = [ p.rtime for p in ladder_allele_sorted ]
    y = [ p.size for p in ladder_allele_sorted ]

    def _f( rtime ):
        """ return (size, deviation)
            deviation is calculated as delta square between curve1 and curve2
        """

        idx = ladder_allele_sorted.bisect_key_right( rtime )

        # left curve
        z1 = np.polyfit( x[idx-2:idx+1], y[idx-2:idx+1], 2)
        size1 = np.poly1d( z1 )(rtime)
        min_score1 = min( x.qscore for x in ladder_allele_sorted[idx-2:idx+1] )

        # right curve
        z2 = np.polyfit( x[idx-1:idx+2], y[idx-1:idx+2], 2)
        size2 = np.poly1d( z2 )(rtime)
        min_score2 = min( x.qscore for x in ladder_allele_sorted[idx-1:idx+2] )

        return ( (size1 + size2)/2, (size1 - size2) ** 2, (min_score1 + min_score2)/2,
                allelemethod.localsouthern)

    return _f

