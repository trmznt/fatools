
# re-imagining the peakutils

import numpy as np

from scipy.signal import find_peaks_cwt
from scipy.optimize import leastsq, curve_fit
from scipy.interpolate import UnivariateSpline

from matplotlib import pyplot as plt

from bisect import bisect_left
from operator import itemgetter
from pprint import pprint

from .dpalign import align_peaks, estimate_z
from .models import *


def set_debug(flag):
    global __printdebug__
    __printdebug__ = flag


def D( *msgs ):
    if __printdebug__:
        print( *msgs )



def filter_peak_number( initial_peaks, max_peak_number, ladder_number ):
    """ filter peaks by height to get number of peaks
        return [ (peaks, height), ... ]
    """

    peak_sets = []

    filtered_peaks = []
    used_min_height = 0
    peak_number = -1

    # gather all possible peak lists

    min_height = 0
    while min_height < 5000:
        min_height += 1
        peaks = [ x for x in initial_peaks if x.height >= min_height ]

        #print('min_height = %d, peaks = %d' % (min_height, len(peaks)))

        if len(peaks) > max_peak_number:
            # too much noises
            continue

        if peak_number < 0:
            peak_number = len(peaks)
            filtered_peaks = peaks
            used_min_height = min_height
            continue

        if peak_number == len(peaks):
            used_min_height = min_height
            continue

        if peak_number != len(peaks):
            peak_sets.append( (used_min_height, filtered_peaks) )
            if len(peaks) < ladder_number:
                peak_sets.append( (min_height, peaks) )
                break

        peak_number = len(peaks)
        filtered_peaks = peaks
        used_min_height = min_height

    return peak_sets


def filter_retention_time( ladders, peaks, parameter ):
    """ return proper filtered peaks in reversed order
    """

    if parameter.init_separation_time < 0:
        separation_slope = 1/estimate_slope( ladders, peaks )
    else:
        separation_slope = parameter.init_separation_time
    print("peak count: %d separation_slope: %2.3f" % (len(peaks), separation_slope))
    min_spacing, max_spacing = spacing( ladders )
    min_time_separation = min_spacing * separation_slope * 0.90
    #print('min time separation: %6.3f' % min_time_separation)
    peaks.reverse()

    ladder_alleles = [ peaks[0] ]
    idx = 1
        
    while idx < len(peaks):
        prev_peak = ladder_alleles.pop()
        #print('comparing peaks: %d & %d' % (prev_peak.peak, peaks[idx].peak))
        if prev_peak.rtime - peaks[idx].rtime > min_time_separation: # minimum time separation
            ladder_alleles.append( prev_peak )
            ladder_alleles.append( peaks[idx] )
        else:
            #print('=> range < min_time_separation')
            if prev_peak.height > peaks[idx].height:
                ladder_alleles.append( prev_peak )
                #print('=> use peak %d' % prev_peak.peak)
            else:
                ladder_alleles.append( peaks[idx] )
                #print('=> use peak %d' % peaks[idx].peak)
        idx += 1

    return ladder_alleles


def peak_fit_and_align( simple_fits, ladders ):

    dp_fits = []
    for (initial_z, initial_rss, initial_peaks, height, paired_peaks) in simple_fits[:6]:

        print('=> initial RSS: %3.3f' % initial_rss)
        dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders,
                                                        initial_peaks,
                                                        initial_z = initial_z )
        dp_fits.append( (dp_score, dp_rss, dp_z, dp_peaks, initial_peaks, height, S, D, paired_peaks) )

    dp_fits.sort( key = itemgetter(1) )
    dp_fits.sort( key = itemgetter(0), reverse=True )

    return dp_fits[0]


def xxx_scan_ladder_peaks( channel, ladders, parameter = None ):

    if parameter is None:
        parameter = LadderScanningParameter()

    parameter.min_height = min( parameter.ladder_height_win )
    initial_peaks = scan_peaks( channel, parameter = parameter )

    peak_sets = filter_peak_number( initial_peaks,
                                    len(ladders) * 1.5,
                                    len(ladders) )

    # sanity check
    if len(peak_sets) == 0 or len(peak_sets[0][1]) < len(ladders) * parameter.min_peak_number:
        return (None, -1, 0, 
            ('Not enough initial signals', 'Too much noises or too few signals'), [])

    # use all peak_sets
    print('peak_sets: %d' % len(peak_sets))

    separated_peaksets = []
    for (height, peaks) in peak_sets:
        proper_peaks = filter_retention_time( ladders, peaks, parameter )
        if not proper_peaks or ( len(proper_peaks) < len(ladders) * parameter.min_peak_number or
                                len(proper_peaks) > len(ladders) * 1.5 ):
            continue
        separated_peaksets.append( (height, proper_peaks) )


    if not separated_peaksets:
        return (None, -1, 0, ('Not enough filtered peaks', 'Bad peak separation time'), [])

    print('separated peaks contains %d sets' % len(separated_peaksets))

    
    simple_fits = []
    for (height, peaks) in separated_peaksets:
        print('=> height: %d, peaks: %d' % (height, len(peaks)))
        fit_results = simple_fit( ladders, peaks )
        for (z, rss, paired_peaks) in fit_results:
            simple_fits.append( (z, rss, peaks, height, paired_peaks) )
    simple_fits.sort( key = itemgetter(1) )

    dp_fits = []
    i = 0
    optimal_dpscore = 0
    optimal_rss = 1000


    while simple_fits[i:i+3] and (optimal_rss > len(ladders)*1.5 and optimal_dpscore < len(ladders)):
        dp_fit = peak_fit_and_align( simple_fits[i:i+3], ladders )
        optimal_dpscore, optimal_rss = dp_fit[0], dp_fit[1]
        dp_fits.append( dp_fit )
        i += 3

    if len(dp_fits) > 1:
        dp_fits.sort( key = itemgetter(1) )
        dp_fits.sort( key = itemgetter(0), reverse=True )

    optimal_dpscore, optimal_rss, optimal_z, optimal_peaks, initial_peaks, height, S, D, paired_peaks = dp_fits[0]

    score, reports = score_ladder( optimal_rss, len(optimal_peaks), len(ladders) )
    reports.append('Height = %d' % height)

    print('=> initial pairing:')
    pprint(paired_peaks)
    #np.savetxt('S.csv', S, delimiter=',')
    #np.savetxt('D.csv', D, delimiter=',')

    # set peaks to ladder-peak
    optimal_peaks.sort()
    for (s, p) in optimal_peaks:
        p.size = p.value = s
        p.type = 'peak-ladder'
        p.method = 'allele-autoladder'
        channel.alleles.append( p )

    print('=> initial peaks')
    pprint(initial_peaks)
    print('=> aligned ladder peaks')
    pprint(optimal_peaks)
    print(' => aligned peaks: %d of %d ladders' % (len(optimal_peaks), len(ladders)))
    return (    optimal_z, optimal_rss, score, optimal_dpscore, reports,
                [ p for (s,p) in optimal_peaks ] )


def scan_ladder_peaks( channel, ladders, parameter = None ):

    if parameter is None:
        parameter = LadderScanningParameter()

    parameter.min_height = min( parameter.ladder_height_win )
    initial_peaks = scan_peaks( channel, parameter = parameter )

    peak_sets = []  # (min_height, peaks)

    filtered_peaks = []
    used_min_height = 0
    peak_number = -1

    # gather all possible peak lists

    if parameter.height > 0:
        # just use peaks with this minimal height
        peak_sets = [ (parameter.height,
                        list( p for p in initial_peaks if p.height > parameter.height ) ) ]

    else:
        peak_sets = filter_peak_number( initial_peaks,
                                    len(ladders) * 2,
                                    len(ladders) )

    # sanity check
    if len(peak_sets) == 0:
        return (None, -1, 0, -1, 
            ('Not enough initial signals', 'Too much noises or too few signals'), [])

    pprint( peak_sets )


    separated_peaks = []
    for (height, peaks) in peak_sets:
        proper_peaks = filter_retention_time( ladders, peaks, parameter )
        if not proper_peaks or ( len(proper_peaks) < len(ladders) * parameter.min_peak_number or
                                len(proper_peaks) > len(ladders) * 1.5 ):
            continue
        separated_peaks.append( (len(proper_peaks), proper_peaks, height) )


    if not separated_peaks:
        return (None, -1, 0, -1,
            ('Not enough filtered peaks', 'Bad peak separation time'), [])

    print('separated peaks contains %d sets' % len(separated_peaks))

    best_fits = []
    for (peak_len, proper_peaks, min_height) in separated_peaks:

        dpscore, rss, z, aligned_peaks = adaptive_peak_alignment( ladders, proper_peaks )
        #rss, z, aligned_peaks = align_peaks( ladders, proper_peaks )
        best_fits.append( (dpscore, rss, z, aligned_peaks, min_height, proper_peaks) )

    #for i in range(len(best_fits)):
    #    rss, z, aligned_peaks, height, proper_peaks = best_fits[i]
    #    dp_rss, dp_z, dp_peaks = align_peaks( ladders, proper_peaks, aligned_peaks )
    #    if dp_rss < rss:
    #        print(' => DP alignment replaced RSS: %3.2f -> %3.2f' % ( rss, dp_rss ))
    #        best_fits[i] = (dp_rss, dp_z, dp_peaks, height, proper_peaks)



    best_fits.sort( key = itemgetter(1) )
    best_fits.sort( key = itemgetter(0), reverse=True )
    print(' => best fits:')
    for (dpscore, rss, z, peaks, height, proper_peaks) in best_fits:
        print('   dpscore: %3.3f rss: %3.2f  peaks: %2d  total peaks: %2d' %
                    (dpscore, rss, len(peaks), len(proper_peaks)) )
        pprint(peaks)
    optimal_dpscore, optimal_rss, optimal_z, optimal_peaks, height, proper_peaks = best_fits[0]

    score, reports = score_ladder( optimal_rss, len(optimal_peaks), len(ladders) )
    reports.append('Height = %d' % height)

    # set peaks to ladder-peak
    optimal_peaks.sort()
    for (s, p) in optimal_peaks:
        p.size = p.value = s
        p.type = 'peak-ladder'

    for p in proper_peaks:
        p.method = 'allele-autoladder'
        if p.value is None:
            p.size = p.value = -1
            p.type = 'peak-unassigned'
        channel.alleles.append( p )

    pprint(proper_peaks)

    print(' => aligned peaks: %d of %d ladders' % (len(optimal_peaks), len(ladders)))
    return (optimal_z, optimal_rss, score, optimal_dpscore, reports, [ p for (s,p) in optimal_peaks ])

    ## NOT USED FROM HERE
    
    # use peaks & min_height later on

    rss, z, aligned_peaks = align_peaks( ladders, filtered_peaks )

    score, reports = score_ladder( rss, len(aligned_peaks), len(ladders) )
    reports.append('Height > %d' % used_min_height)

    # set peaks to ladder-peak
    aligned_peaks.sort()
    for (s, p) in aligned_peaks:
        p.size = p.value = s
        p.type = 'peak-ladder'
        p.method = 'allele-autoladder'
        channel.alleles.append( p )

    print('=> aligned peaks: %d of %d ladders' % (len(aligned_peaks), len(ladders)))
    return (z, rss, score, reports, [ p for (s,p) in aligned_peaks ])


def scan_ladder_peaks_xxx( channel, ladders, parameter = None):
    """ scan and assign ladder peaks
        ladders: [ 100, 150, 200, ... ]
        return: (z, rss, score, report, assigned_peaks)
    """

    if parameter is None:
        parameter = LadderScanningParameter()

    parameter.min_height = min( parameter.ladder_height_win )
    initial_peaks = scan_peaks( channel, parameter = parameter )

    filtered_peaks = []

    prev_peak_len = -1

    for min_height in parameter.ladder_height_win:

        peaks = [ x for x in initial_peaks if x.height >= min_height ]

        if len(peaks) == prev_peak_len:
            continue

        if len(peaks) > len(ladders) * 2:
            # too much noises
            continue

        if len(peaks) < len(ladders) * parameter.min_peak_number:
            # not much signals
            break

        filtered_peaks.append( (min_height, peaks) )
        prev_peak_len = len(peaks)

    if not filtered_peaks:
        # no available data
        return (None, -1, 0, 
            ('Not enough initial signals', 'Too much noises or too few signals'), [])

    # debugging
    if True:
        for fp in filtered_peaks:
            print('height = %d, peaks = %d' % (fp[0], len(fp[1])))
            for peak in fp[1]:
                print(peak)


    separated_peaks = []
    for (min_height, peaks) in filtered_peaks:
        proper_peaks = filter_retention_time( ladders, peaks, parameter )
        if not proper_peaks or len(proper_peaks) < len(ladders) * parameter.min_peak_number:
            continue
        separated_peaks.append( (len(proper_peaks), proper_peaks, min_height) )


    if not separated_peaks:
        return (None, -1, 0, ('Not enough filtered peaks', 'Bad peak separation time'), [])

    ladders = list(reversed(ladders))
    best_fits = []

    for (peak_len, proper_peaks, min_height) in separated_peaks:

        rss, z, aligned_peaks = adaptive_peak_alignment( ladders, proper_peaks )
        #rss, z, aligned_peaks = align_peaks( ladders, proper_peaks )
        best_fits.append( (rss, z, aligned_peaks, min_height, proper_peaks) )

    best_fits.sort( key = itemgetter(0) )
    optimal_rss, optimal_z, optimal_peaks, min_height, proper_peaks = best_fits[0]

    dp_rss, dp_z, dp_peaks = align_peaks( ladders, proper_peaks, optimal_peaks )

    if dp_rss < optimal_rss:
        optimal_rss, optimal_z, optimal_peaks = dp_rss, dp_z, dp_peaks
        print(' => DP alignment replaced optimal RSS')

    score, reports = score_ladder( optimal_rss, len(optimal_peaks), len(ladders) )
    reports.append('Height > %d' % min_height)

    # set peaks to ladder-peak
    optimal_peaks.sort()
    for (s, p) in optimal_peaks:
        p.size = p.value = s
        p.type = 'peak-ladder'
        p.method = 'allele-autoladder'
        channel.alleles.append( p )

    print('=> aligned peaks: %d of %d ladders' % (len(optimal_peaks), len(ladders)))
    pprint(optimal_peaks)
    return (optimal_z, optimal_rss, score, reports, [ p for (s,p) in optimal_peaks ])


def filter_separation_time( ladders, peaks, parameter ):
    """ return proper filtered peaks in reversed order
    """

    separation_slope = estimate_slope( ladders, peaks )
    min_spacing, max_spacing = spacing( ladders )
    min_time_separation = min_spacing * separation_slope * 0.90
    #print('min time separation: %6.3f' % min_time_separation)
    peaks.reverse()

    ladder_alleles = [ peaks[0] ]
    idx = 1
        
    while idx < len(peaks):
        prev_peak = ladder_alleles.pop()
        #print('comparing peaks: %d & %d' % (prev_peak.peak, peaks[idx].peak))
        if prev_peak.peak - peaks[idx].peak > min_time_separation: # minimum time separation
            ladder_alleles.append( prev_peak )
            ladder_alleles.append( peaks[idx] )
        else:
            #print('=> range < min_time_separation')
            if prev_peak.height > peaks[idx].height:
                ladder_alleles.append( prev_peak )
                #print('=> use peak %d' % prev_peak.peak)
            else:
                ladder_alleles.append( peaks[idx] )
                #print('=> use peak %d' % peaks[idx].peak)
        idx += 1

    return ladder_alleles


def adaptive_peak_alignment( ladders, peaks ):
    """ return (rss, z, assigned_peaks)
    """

    # sort in reverse order
    ladders = sorted( ladders, reverse = True )
    peaks = sorted( peaks, key = lambda x: x.rtime, reverse = True )

    adaptive_results = []
    for i in range(-3, 8):
        if i >= 0:
            result_peaks = list(zip(ladders[i:], peaks))
        else:
            result_peaks = list(zip(ladders, peaks[ abs(i): ]))
        z, rss = estimate_z( result_peaks )

        adaptive_results.append( (rss, z, result_peaks) )
        #print('i=%d, rss=%4.3f' % (i, rss))

    # we use the smallest rss
    adaptive_results.sort()
    dp_results = []

    # only use the best 3 RSS result

    for (rss, z, aligned_peaks) in adaptive_results[:3]:

        #if rss > 5000:
            # too much deviation, just ignore
        #    continue

        # perform dynamic programming alignment

        dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders, peaks, aligned_peaks )
        if dp_rss < rss:
            print(' => DP alignment replaced RSS: %3.2f -> %3.2f' % ( rss, dp_rss ))
        dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks) )
        #dp_results.append( (rss, z, aligned_peaks) )

    if len(dp_results) <= 0:
        return adaptive_results[0]

    dp_results.sort( key = itemgetter(1) )
    dp_results.sort( key = itemgetter(0), reverse=True  )

    return dp_results[0]



def scan_ladder_peaks_xxx( channel, ladders, parameter=None):
    """ return a list of ladder alleles
        ladders: [ 100, 150, 200, ... ]
        return: (z, rss, score, report, ladder_alleles)
    """

    if parameter is None:
        parameter = LadderScanningParameter()

    parameter.min_height = min( parameter.ladder_height_win )
    minheight_peaks = scan_peaks( channel, method=None, parameter=parameter )
    prev_peaks = []
    last_peak_no = -1
    min_spacing, max_spacing = spacing( ladders )


    # find the proper height to get the optimum number of peaks
    for min_height in parameter.ladder_height_win:
        
        print('processing for min_height = %d' % min_height)

        peaks = [ peak for peak in minheight_peaks if peak.height >= min_height ]
        if len(peaks) > len(ladders) + 25:
            # too much noise
            print('too much noises')
            continue

        if len(peaks) < 5:
            # too few signal
            return (None, -1, 0, ("Too few signals",), None)

        if last_peak_no == len(peaks):
            # if current peak number is similar to above, just continue increase
            # min_height
            continue

        last_peak_no = len(peaks)

        print('initial peaks: %d' % len(peaks))
        pprint( [ (x.peak, x.height) for x in peaks ] )

        # based on the longest peak, estimate the minimum time separation range
        separation_slope = estimate_slope( ladders, peaks )
        min_time_separation = min_spacing * separation_slope * 0.90
        print('min time separation: %6.3f' % min_time_separation)
        peaks.reverse()

        ladder_alleles = [ peaks[0] ]
        idx = 1
        
        while idx < len(peaks):
            prev_peak = ladder_alleles.pop()
            #print('comparing peaks: %d & %d' % (prev_peak.peak, peaks[idx].peak))
            if prev_peak.peak - peaks[idx].peak > min_time_separation: # minimum time separation
                ladder_alleles.append( prev_peak )
                ladder_alleles.append( peaks[idx] )
            else:
                #print('=> range < min_time_separation')
                if prev_peak.height > peaks[idx].height:
                    ladder_alleles.append( prev_peak )
                    #print('=> use peak %d' % prev_peak.peak)
                else:
                    ladder_alleles.append( peaks[idx] )
                    #print('=> use peak %d' % peaks[idx].peak)
            idx += 1

        peaks = ladder_alleles
        print('after range filtering: %d' % len(peaks))
        print( [x.peak for x in peaks] )


        print('min height: %d -- peaks: %d' % (min_height, len(peaks)))

        if len(peaks) <= len(ladders):
            if not prev_peaks:
                if len(ladders) - len(peaks) >= 10:
                    return ( None, -1, 0, [ 'Missing peaks > 10' ], None )
                break
            # check which number is smaller, additional peaks or missing peaks
            if abs(len(peaks) - len(ladders)) > abs(len(prev_peaks) - len(ladders)):
                print('using previous peaks!!')
                peaks = prev_peaks
            break

        prev_peaks = peaks

    #peaks = list(reversed(peaks))
    print('ASSIGNING PEAKS')
    print('ladder peaks: %d of %d' % (len(peaks), len(ladders)))



    # perform sizing
    # assumption:
    # - no drop peaks
    # - the maximum number of missing peaks at the later retention time is 5 (shifted out)

    #peaks = align_ladder_peaks( ladders, ladder_alleles )
    #return

    ladders = list( reversed( ladders ) )
    peaks = list( sorted( peaks, key = lambda x: x.peak, reverse=True) )
    adaptive_results = []
    for i in range(-3, 8):
        if i >= 0:
            result_peaks = list(zip(ladders[i:], peaks))
        else:
            result_peaks = list(zip(ladders, peaks[ abs(i): ]))
        pprint( result_peaks )
        z, rss = estimate_z( result_peaks )
        adaptive_results.append( (result_peaks, z, rss) )
        print('i=%d, rss=%4.3f' % (i, rss))

    # we use the smallest rss
    adaptive_results.sort( key = itemgetter(2) )
    (optimal_peaks, optimal_z, optimal_rss) = adaptive_results[0]
    optimal_peaks = list(optimal_peaks)

    # scoring
    quality_score = 1.0
    reports = []
    if optimal_rss <= 0:
        quality_score -= 0.25
        reports.append( 'QC RSS < 0' )
    if optimal_rss > 50:
        quality_score -= 0.05
        reports.append( 'QC RSS > 50' )
    if optimal_rss > 255:
        quality_score -= 0.10
        reports.append( 'QC RSS > 255')
    if optimal_rss > 1024:
        quality_score -= 0.10
        reports.append( 'QC RSS > 1024' )
    missing_peaks = len(ladders) - len(optimal_peaks)
    if missing_peaks > 0:
        quality_score -= 0.025 * (missing_peaks - 1)
        reports.append( 'QC missing %d peaks' % missing_peaks )
    reports.reverse()

    # assign peaks
    optimal_peaks.sort()
    for (s, p) in optimal_peaks:
        p.size = p.value = s
        p.type = 'peak-ladder'
        p.method = 'allele-autoladder'
        channel.alleles.append( p )

    return (optimal_z, optimal_rss, quality_score, reports, [ p for (s,p) in optimal_peaks ])


def estimate_z_xxx( peak_pairs ):
    """ estimate z and rss based on [y, x]
        x: peaks
        y: ladder sizes
        return (z, rss)
    """

    x, y = [], []
    for (s, p) in peak_pairs:
        x.append( p.rtime )
        y.append( s )

    z = np.polyfit( x, y, 3 )
    p = np.poly1d( z )
    y_p = p(x)
    rss = ( (y_p - y) ** 2 ).sum()

    return z, rss


def score_ladder( optimal_rss, peak_len, ladder_len):
    """ return ladder score
    """
    # scoring
    quality_score = 1.0
    reports = []
    if optimal_rss <= 0:
        quality_score -= 0.25
        reports.append( 'QC RSS < 0' )
    if optimal_rss > 50:
        quality_score -= 0.05
        reports.append( 'QC RSS > 50' )
    if optimal_rss > 255:
        quality_score -= 0.10
        reports.append( 'QC RSS > 255')
    if optimal_rss > 1024:
        quality_score -= 0.10
        reports.append( 'QC RSS > 1024' )
    missing_peaks = ladder_len - peak_len
    if missing_peaks > 0:
        quality_score -= 0.025 * (missing_peaks - 1)
        reports.append( 'QC missing peaks = %d' % missing_peaks )
    reports.reverse()

    return quality_score, reports


def scan_peaks( channel, parameter=None ):
    """ return a list of alleles after scanned
    """

    if parameter is None:
        parameter = ScanningParameter()

    peaks = find_peaks( channel.data, parameter, channel.data )
    # peaks is [ (x, height, area), ... ]

    if parameter.max_peak_number > 0:
        min_height = 0
        while len(peaks) > parameter.max_peak_number and min_height < 100:
            min_height += 1
            peaks = [ p for p in peaks if p[1] > min_height ]

    # create alleles based on these peaks
    alleles = []
    for peak in peaks:
        right_skew = peak[4] - peak[0]
        if right_skew == 0:
            right_skew = 1
        left_skew = peak[0] - peak[3]
        if left_skew == 0:
            left_skew = 1
        allele = channel.get_allele_class()(    rtime = int(peak[0]),
                                                height = int(peak[1]),
                                                area = peak[2],
                                                brtime = peak[3],
                                                ertime = peak[4],
                                                wrtime = peak[4]-peak[3],
                                                srtime = right_skew/left_skew,
                                                beta = peak[2]/int(peak[1])
                                            )
        allele.type = 'peak-scanned'
        allele.method = 'binning-unavailable'
        alleles.append( allele )

    return alleles


def call_peaks_XXX( channel, method, parameter=None, peaks=None):
    pass

    if method is None:
        # presumably ladder channel
        # just return the alleles
        return alleles

    else:

        # estimate size of the alleles
        for allele in alleles:
            estimate_allele_size( allele, method )

        # remove small fragment size
        alleles = [ allele for allele in alleles if allele.size > parameter.min_size ]

        # classify alleles to peak-stutter or peak-called
        check_stutter_peaks( alleles, threshold = parameter.stutter_threshold )

        return alleles


def call_peaks( alleles, method, ladders, parameter=None):
    """ return a list of called (sized) alleles
        
    """
    
    called_alleles = []
    for allele in alleles:

        # only call peaks within the ladders, otherwise just assigned -1 to the value
        if not ( ladders[0].rtime < allele.rtime < ladders[-1].rtime ):
            allele.size = -1
            allele.type = 'peak-unassigned'
            continue

        size = float(method(allele.rtime))
        if np.isnan(size):
            allele.size = -1
        else:
            allele.size = size
        allele.value = round(allele.size)
        allele.method = 'binning-unavailable'

        # remove small fragment size
        if allele.size < parameter.min_size or allele.size > parameter.max_size:
            allele.type = 'peak-unassigned'
            continue
        allele.type = 'peak-called'
        called_alleles.append( allele )

    # classify alleles to peak-stutter or peak-called for the list of alleles that
    # were not peak-unassigned
    check_stutter_peaks( called_alleles, threshold = parameter.stutter_threshold )

    return called_alleles


def bin_peaks( channel, marker ):
    """ bin peaks that have been called (have peak size > 0),
        return the list of alleles, after allele has been set as peak-bin or peak-artifact
    """

    alleles = list(channel.alleles)

    for allele in alleles:

        if allele.size > 0:
            if allele.type != 'peak-overlap':
                if allele.type != 'peak-called':
                    allele.type = 'peak-called'
        else:
            continue

        #if allele.type == 'peak-unassigned':
        #    continue

        if not marker.min_size < allele.size < marker.max_size:
            allele.type == 'peak-unassigned'
            continue

        if allele.type not in ['peak-called', 'peak-overlap']:
            continue

        size = allele.size
        binlist = marker.bins
        binpos = list( x[0] for x in binlist )
        threshold = float(marker.repeats) / 2 * 1.5

        pos = bisect_left( binpos, size )

        if pos == 0:
            value = binlist[0]
        elif pos == len(binlist):
            value = binlist[-1]
        else:
            before = binlist[pos - 1]
            after = binlist[pos]
            if after[0] - size < size - before[0]:
                value = after
            else:
                value = before
        if abs(value[0] - size) > threshold:
            print('WARN: binned peak with size: %3.2f for value: %3.2f is above range threshold: %2.1f'
                    % (size, value[0], threshold) )

        allele.value = value[1]
        allele.type = 'peak-bin'


        

def check_stutter_peaks( alleles, threshold ):
    """ assign allele.type to peak-stutter or peak-called
    """

    # sort the alleles to ensure it is sorted !
    alleles.sort( key = lambda x: x.size )

    for idx in range( len(alleles) ):

        allele = alleles[idx]

        if idx > 0:
            allele_0 = alleles[idx-1]
            if allele.size - allele_0.size < threshold:
                if allele_0.height > allele.height:
                    allele.type = 'peak-stutter'
        if idx < len(alleles) - 1:
            allele_1 = alleles[idx+1]
            if allele_1.size - allele.size < threshold:
                if allele_1.height > allele.height:
                    allele.type = 'peak-stutter'


def is_overlap(peak_1, peak_2):
    
    if peak_1.brtime > peak_2.brtime:
        peak_1, peak_2 = peak_2, peak_1

    if peak_2.brtime < peak_1.ertime:
        return True

    return False


def check_overlap_peaks( channels, threshold ):
    """ assign allele.type to peak-overlap
    """

    channel_peaks = [ list(channel.alleles) for channel in channels ]

    for i, peaks in enumerate( channel_peaks ):
    
        for j in range( len(peaks) ):
            peak = peaks[j]

            if peak.type not in ['peak-called', 'peak-bin']:
                continue

            for k in range( len(channel_peaks) ):

                if k == i:
                    continue

                for peak_r in channel_peaks[k]:
                    if peak_r.type not in ['peak-called', 'peak-bin']:
                        continue

                    if is_overlap( peak, peak_r ):
                        # find whether the height of any peak is inside any peak
                        channel = peak.alleleset.channel
                        channel_r = peak_r.alleleset.channel
                        if (peak.height < channel_r.data[ peak.rtime ] or
                            peak_r.height < channel.data[ peak_r.rtime ]):
                                if peak.height < peak_r.height:
                                    peak.type = 'peak-overlap'

                    #if abs(peak.size - peak_r.size) < threshold:
                    #    if peak.height < peak_r.height:
                    #        peak.type = 'peak-overlap'


def find_peaks( signal, params, ssignal = None ):
    """ find peaks from signal
        returns [ (peak, height, area), ... ]
        TODO: this method can be perform remotely!
    """
    #print('Finding peaks with min_height = %d' % params.min_height)

    # find all peaks by cwt-based algorithm with Signal-Noise Ratio = 1
    indices = find_peaks_cwt( signal, params.cwt_widths, min_snr = params.cwt_min_snr )
    print('find_peaks_cwt() found peaks: %d' % len(indices))

    if not indices:
        return []

    # filter for absolute heights

    raw_peaks = []
    for idx in indices:
        for i in range(3, -1, -1):
            try:
                height, index = max( [ (signal[i], i) for i in range(idx-3, idx+3) ] )
            except IndexError:
                continue
            break
        if height < params.min_height:
            continue
        if index < 0:
            continue
        raw_peaks.append( (index, height) )

    #pprint(raw_peaks)

    # check for any peaks
    if not raw_peaks:
        return raw_peaks

    # filter for relative heights with median of peak height

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


    if ssignal is None:
        ssignal = smooth_signal( signal )

    # calculate area

    threshold = np.percentile( ssignal, 75 )
    if threshold < 0:
        threshold = -1
    peaks = []
    #pprint(raw_peaks)
    for (peak, height) in raw_peaks:
        area, brtime, ertime = calculate_area( ssignal, peak, 5e-2 )
        peaks.append( (peak, height, area, brtime, ertime) )

    #print('raw peaks after height filtering:')
    #pprint( [ (x[0], x[1], x[2]) for x in peaks ] )


    # filter for peak area


    return peaks


def calculate_area_xxx(y, p, threshold):

        area = 0
        x = p + 1
        if x < len(y):
            h = y[x]
            while h > threshold:
                area += h
                x = x + 1
                if x >= len(y):
                    break
                h = y[x]

        x = p - 1
        if x >= 0:
            h = y[x]
            while h > threshold:
                area += h
                x = x - 1
                if x < 0:
                    break
                h = y[x]

        area += y[ p ]
        
        return int(area)


def calculate_area(y, t, threshold):
    """ return (area, brtime and ertime)
    """

    # right area
    data = y[t:]
    r_area, ertime, r_shared = half_area(data, threshold)

    # left area
    data = y[:t+1][::-1]
    l_area, brtime, l_shared= half_area(data, threshold)

    return ( l_area + r_area - y[t], t - brtime, ertime + t )



def half_area(y, threshold):
    """ return (area, ertime)
    """

    winsize = 3
    threshold = threshold/2
    shared = False
    area = y[0]
    edge = float(np.sum(y[0:winsize]))/winsize
    old_edge = 2 * edge

    index = 1
    limit = len(y)

    while edge > area * threshold and edge < old_edge and index < limit:
        old_edge = edge
        area += y[index]
        edge = float(np.sum(y[index:index+winsize]))/winsize
        index += 1
    if edge >= old_edge:
        shared = True
    index -= 1

    return area, index, shared


def spaces( a ):
    """ return spaces between members
    """
    return [ j-i for i,j in zip( a[:-1], a[1:] ) ]


def spacing( a ):
    """ return the min and max length between consecutive elements """
    spaces = [ j-i for i,j in zip( a[:-1], a[1:] ) ]
    return min( spaces ), max( spaces )


def least_square( z ):
    """ 3rd order polynomial resolver
    """
    return np.poly1d( z )


def cubic_spline( ladder_alleles ):
    """ cubic spline interpolation
        x is peaks, y is standard size
    """
    ladder_alleles = list(sorted(ladder_alleles, key = lambda x: x.size))

    ladder_peaks = []
    ladder_sizes = []
    for ladder_allele in ladder_alleles:
        ladder_peaks.append( ladder_allele.rtime )
        ladder_sizes.append( ladder_allele.size )
    return UnivariateSpline(ladder_peaks, ladder_sizes, s=3)


def estimate_allele_size( allele, method ):
    """ estimate allele size based on the method
    """

    allele.size = float(method(allele.rtime))
    allele.value = int(allele.size)
    allele.type = 'peak-called'
    allele.method = 'binning-unavailable'


def align_ladder_peaks_xxx( ladders, peaks, indels_penalty = (0,0) ):
    """ align ladders & peaks
        return: [ (ladder_size, peak), ... ]
    """

    # sort ladders & peaks in reverse orders
    ladders = reversed( sorted(ladders) )
    peaks = reversed( sorted( peaks, key = lambda x: x.peak ) )

    # create a normalized data structure
    norm_ladders = [ (x/ladders[0], x) for x in ladders ]
    norm_peaks = [ (x.peak/peaks[0].peak, x) for x in peaks ]

    # matrix n+1 * m+1, n ~ length of ladders, m ~ length of peaks
    # matrix scoring
    n = len(norm_ladders) + 1
    m = len(norm_peaks) + 1
    S = np.zeros((n,m))     # score
    T = np.zeros((n,m))     # traceback


def estimate_slope(ladders, peaks):
    """ estimate slope based on 50% sample data """

    # sort ladders & peaks in reverse orders
    ladders = list(reversed( sorted(ladders) ))
    peaks = list(reversed( sorted( peaks, key = lambda x: x.rtime ) ))

    # get 25% - 75% sample size
    N_ladders = len(ladders)
    N_peaks = len(peaks)

    range_ladders = sorted( spaces( ladders ) )
    range_peaks = sorted( spaces( [ x.rtime for x in peaks ] ) )

    sampled_range_ladders = range_ladders[int(N_ladders * 0.33) : int(N_peaks * 0.66)]
    sampled_range_peaks = range_peaks[int(N_peaks * 0.33) : int(N_peaks * 0.66) ]

    slope = np.mean( sampled_range_peaks ) / np.mean( sampled_range_ladders )
    return 1.0/slope


def align_ladder_peaks( ladders, peaks ):
    """ align ladders & peaks
        return [ (ladder_size, peak), ... ]
    """

    # sort ladders & peaks in reverse orders
    ladders = list(reversed( sorted(ladders) ))
    peaks = list(reversed( sorted( peaks, key = lambda x: x.peak ) ))

    # get 25% - 75% sample size
    N_ladders = len(ladders)
    N_peaks = len(peaks)

    range_ladders = sorted( spaces( ladders ) )
    range_peaks = sorted( spaces( [ x.peak for x in peaks ] ) )

    sampled_range_ladders = range_ladders[int(N_ladders * 0.25) : int(N_peaks * 0.75)]
    sampled_range_peaks = range_peaks[int(N_peaks * 0.25) : int(N_peaks * 0.75) ]

    slope = np.mean( sampled_range_peaks ) / np.mean( sampled_range_ladders )

    # linear curve fit with the approximate slope
    # y -> peaks, x -> size

    initial_y = make_linear_func( slope )

    # iteration for simulating missing last peak(s)
    for i in range(-3, 3):

        z, rss, peak_assignment = fit_ladder_peaks( ladders, peaks, initial_y, i )

        print('RESULTS:')
        print(z, rss)
        pprint(peak_assignment)


def fit_ladder_peaks( ladders, peaks, initial_y, index):

    # initial fit firs (linear curve fit)
    param, covar, rss = linear_fit_ladder_peaks( ladders, peaks, initial_y, index )

    # reassign ladders with peaks for consecutive peakfit
    peak_assignment = assign_ladder_peak( ladders, peaks, lambda x: initial_y(x, param) )

    # use polynomial_fit
    z, rss = estimate_z( peak_assignment )

    return (z, rss, peak_assignment)


def assign_ladder_peak( ladders, peaks, func_y ):

    assigned_peaks = {}

    for ladder_size in ladders:
        y_ladder = func_y( ladder_size )
        peak, delta = find_closest_peak( peaks, y_ladder )
        if peak in assigned_peaks:
            (other_ladder, other_delta) = assigned_peaks[peak]
            if other_delta > delta:
                assigned_peaks[peak] = (ladder_size, delta)
        else:
            assigned_peaks[peak] = (ladder_size, delta)

    peak_assignment = [ (v[0], k) for (k, v) in assigned_peaks.items() ]
    peak_assignment.sort()

    return peak_assignment


def find_closest_peak( peaks, y ):

    min_d = abs(peaks[0].peak - y)
    min_peak = peaks[0]

    for peak in peaks[1:]:
        delta = abs(peak.peak - y)
        if delta < min_d:
            min_d = delta
            min_peak = peak

    return min_peak, min_d


def linear_fit_ladder_peaks( ladders, peaks, initial_y, idx):

    if idx >= 0:
        ladders = ladders[idx:]
    else:
        peaks = peaks[ abs(idx): ]

    data_length = min( len(ladders), len(peaks) )
    p0 = 1
    x = np.array( ladders[:data_length] )
    y = np.array( [ p.peak for p in peaks ][:data_length] )
    #linsq = leastsq( initial_residuals, p0, args = ( x, y)  )
    param, covar = curve_fit( initial_y, x, y, p0 )

    # calculate RSS
    y_p = initial_y(x, param)
    print('y_p:', y_p)
    print('x:', x)
    rss = ( (y_p - y) ** 2 ).sum()

    return (param, covar, rss)

    # plot the fitt
    print('idx = %d, param = %8.3f, covar = %8.3f, rss = %8.3f' % (idx, param, covar, rss))
    x_lines = np.linspace( x[0], x[-1], 1000 )
    y_lines = initial_y(x_lines, param)
    plt.figure()
    plt.plot( x, y, 'ro', x_lines, y_lines )
    plt.savefig('plot#%d.png' % idx)
    # based on this linsq, reassigned all peaks

    #raise RuntimeError(linsq)
    return (param, covar, rss)


def make_linear_func( slope = 1 ):
    def _y(x, B):
        return (x * slope) + B
    return _y


def make_residual_func( linear_func ):
    def _residuals(B, y, x):
        return (y - linear_func(x, B))
    return _residuals


def linear_fit(ladders, peaks):
    """ make a linear fit on all possible alignment
        return [ (z, rss), ... ]
    """

    # assume that ladders and peaks are already sorted in ascending orders

    # calculate the linear slope
    A = estimate_slope( ladders, peaks )
    linear_func = make_linear_func( A )
    print('=> estimated linear slope: %3.3f with peaks: %d' % (A, len(peaks)))

    # fit the slope to data

    linear_fits = []    # store [ (z, rss) ]
    for i in range(-10, 10):

        param, rss = align_and_fit( ladders, peaks, linear_func, i )
        linear_fits.append( ( np.array([A, param[0]]), rss ) )

    linear_fits.sort( key = itemgetter(1) )

    return linear_fits
        
        
def align_and_fit( ladders, peaks, func, idx ):
    """ align and fit the ladders to peaks
        func -> size[bp] = func( rtime[sec] )
            x : rtime
            y : size or ladder
    """

    if idx >= 0:
        ladders = ladders[idx:]
    else:
        peaks = peaks[ abs(idx): ]

    data_length = min( len(ladders), len(peaks) )
    p0 = 1
    y = np.array( ladders[:data_length] )
    x = np.array( [ p.rtime for p in peaks ][:data_length] )
    #linsq = leastsq( initial_residuals, p0, args = ( x, y)  )
    param, covar = curve_fit( func, x, y, p0 )

    # calculate RSS
    y_p = func(x, param)
    #print('y_p:', y_p)
    #print('x:', x)
    rss = ( (y_p - y) ** 2 ).sum()

    return (param, rss)


def simple_fit(ladders, peaks):
    """ return [ (z, rss), ... ] """

    fits = []
    for idx in range(-5, +5):
        
        if idx >= 0:
            pair_peaks = list(zip( ladders[idx:], peaks ))
        else:
            pair_peaks = list(zip( ladders, peaks[ abs(idx): ] ))

        z, rss = estimate_z( pair_peaks )
        fits.append( (z, rss, pair_peaks) )

    #fits.sort( key = itemgetter(1) )
    
    return fits


# ---------------------------------------------------------------------------------------
# new strategy for peak alignment
#

def find_ladder_peak(channel, N, parameter=None):
    """ return peakset with max_N number
        return values: [ (used_peak_height, peaks), ... ]
    """

    if parameter is None:
        parameter = LadderScanningParameter()

    parameter.min_height = min( parameter.ladder_height_win )
    initial_peaks = scan_peaks( channel, parameter = parameter )

    # gather all possible peak lists

    if parameter.height > 0:
        # just use peaks with this minimal height
        return [ (parameter.height,
                        list( p for p in initial_peaks if p.height > parameter.height ) ) ]

    else:
        return = filter_peak_number( initial_peaks,
                                    N + parameter.additional_peaks,
                                    N )
    


def scan_ladder_peak_2(channel, ladder, parameter=None, find_peaks = True):

    if parameter = None:
        parameter = LadderScanningParameter()

    if find_peaks:
        peak_sets = find_ladder_peak( channel, len(ladder), parameter )
    else:
        # use the old peaks, set all into peak-unassigned, and redo the peak alignment
        pass


    results = []

    for peak_set in peak_sets:

        results.append( align_ladder_peaks( peak_set, ladder, parameter ) )
        
        # find similarity using PCA-based similarity

        # assign the first peak set for the first iteration of DP

        # do the DP alignment, get the result
        # use eigenvalues of cosine product to evaluate peak similarity


    # for each results, find the least RSS and the highest DP score and quality score

    # set the peaks for peak-ladders


def align_ladder_peaks( peak_set, ladder, parameter ):

    pass
