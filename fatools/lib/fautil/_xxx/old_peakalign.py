
from fatools.lib import dpalign as dp
from fatools.lib.const import *
from fatools.lib.utils import cout, cerr
import numpy as np
import pprint, sys


def fast_align( trace, initial_peaks, avg_height, ladders, all_peaks ):

    # simple fast aligner

    #if len(initial_peaks) < len(ladders):
    #    raise RuntimeError("shouldn't be here")
    peak_pairs = list(zip( reversed( [p.rtime for p in initial_peaks] ), reversed( ladders ) ))
    peak_pairs.sort()
    #print('fast_align with %d peaks.' % len(peak_pairs))
    #pprint.pprint(peak_pairs)
    if len(initial_peaks) < len(ladders):
        initial_peaks = all_peaks
    dp_score, dp_rss, dp_z, dp_peaks, S, D = align(initial_peaks, ladders, peak_pairs)
    #print(' * dp_score =', dp_score )
    #pprint.pprint(dp_peaks)
    return (dp_score, dp_rss, dp_z, dp_peaks)


def greedy_align( trace, initial_peaks, avg_height, ladders, all_peaks ):

    # greedy aligner

    return greedy_1c( trace, initial_peaks, avg_height, ladders, all_peaks, 50)




def greedy_1d( trace, initial_peaks, avg_height, ladders, all_peaks, min_rss ):

    data = trace

    delta = len(initial_peaks) - len(ladders)
    if delta > 0:
        # there is possibility of having noise peaks
        I, J = 2 + delta, 3 + delta
    elif delta < 0:
        # there is possibility of losing ladder peaks
        I, J = 2, 3
        ladders = ladders[:delta]
    else:
        I, J = 2, 3

    dp_results = []

    for i in range(0, I):
        for j in range(1, J):

            # prepare initial peak-size pairing
            peaks = initial_peaks
            peak_pairs = generate_initial_peak_pairs( peaks, ladders, i, j )
            pprint.pprint(peak_pairs)

            last_dp = last_rss = -1
            last_z = last_peaks = last_S = last_D = None
            while True:

                dp_score, dp_rss, dp_z, dp_peaks, S, D = align(peaks, ladders, peak_pairs)
                print(' * dp_score =', dp_score )

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

                cerr(' ** ===> reiterate peak alignment' )
                last_dp = dp_score
                last_rss = dp_rss
                last_z = dp_z
                last_peaks = dp_peaks
                last_S = S
                last_D = D

                # recreate peak_pairs based on db_peaks
                peak_pairs = [ (x[1].rtime, x[0]) for x in dp_peaks ]
                #if dp_rss > min_rss:
                #peaks = all_peaks
                #peak_pairs = reassign_peaks( peaks, ladders, dp_z )
                #pprint.pprint(peak_pairs)

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]

    return dp_score, dp_rss, dp_z, dp_peaks



def greedy_1c( trace, initial_peaks, avg_height, ladders, all_peaks, min_rss ):
    """ the best peak-size alignment yet, but slow """

    #pprint.pprint( peaks )

    # prepare data
    data = trace
    delta = len(initial_peaks) - len(ladders)
    if delta > 0:
        # there is possibility of having noise peaks
        I, J = 2 + delta, 3 + delta
    elif delta < 0:
        # there is possibility of losing ladder peaks
        I, J = 2, 3
        ladders = ladders[:delta]
    else:
        I, J = 2, 3

    dp_results = []

    for i in range(0, I):
        for j in range(1, J):

            # prepare initial peak-size pairing
            peaks = initial_peaks
            peak_pairs = generate_initial_peak_pairs( peaks, ladders, i, j )
            #pprint.pprint(peak_pairs)

            last_dp = last_rss = -1
            last_z = last_peaks = last_S = last_D = None
            while True:

                dp_score, dp_rss, dp_z, dp_peaks, S, D = align(peaks, ladders, peak_pairs)
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

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]

    return dp_score, dp_rss, dp_z, dp_peaks


def reassign_peaks( peaks, ladders, z ):

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


def reassign_peaks_xxx( peaks, ladders, z ):
    """ using z, reassign peak-size pairing [ (rtime, size), ... ] """

    print(' => Z:', z)
    ladder_dict = {}
    f = np.poly1d(z)
    for p in peaks:
        s = f(p.rtime)
        # find the closest size
        d = [ ( abs(x - s), x ) for x in ladders ]
        d.sort()
        d0, s0 = d[0]
        print(' **-> %3d <> %4d' % (s0, p.rtime))
        if s0 in ladder_dict:
            if ladder_dict[s0][0] > d0:
                ladder_dict[s0] = ( d0, p.rtime )
        else:
            ladder_dict[s0] = ( d0, p.rtime )

    peak_pairs = []
    for (s, data) in ladder_dict.items():
        peak_pairs.append( (data[1], s) )
    peak_pairs.sort()

    return peak_pairs


def greedy_1b( trace, peaks, avg_height, ladders, all_peaks, min_rss ):

    data = trace

    # make pairing of peaks & ladders

    delta = len(peaks) - len(ladders)
    if delta > 0:
        # there is possibility of having noise peaks
        I, J = 2 + delta, 3 + delta
    elif delta < 0:
        # there is possibility of losing ladder peaks
        I, J = 2, 3
        ladders = ladders[:delta]
    else:
        I, J = 2, 3

    dp_results = []

    for i in range(0, I):
        for j in range(1, J):


            # repairing peaks
            used_peaks = peaks[i:len(peaks) - j + 1]
            print('Peaks: %d => %d' % (len(peaks), len(used_peaks)))
            peak_pairs = []
            min_size = round( min( len(used_peaks), len(ladders) ) / 2 )
            for idx in range(min_size):
                peak_pairs.append( (used_peaks[idx].rtime, ladders[idx]) )
                peak_pairs.append( (used_peaks[-1-idx].rtime, ladders[-1-idx]) )
            peak_pairs.sort()
            #pprint.pprint( peak_pairs )

            dp_score, dp_rss, dp_z, dp_peaks, S, D = align( peaks, ladders, peak_pairs )

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    #for dp_result in dp_results:
    #    print(' =>', dp_result[0], dp_result[1], len(dp_result[3]) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]

    if dp_rss > 70:
        
        # try again
        peak_pairs = [ (peak.rtime, size) for size, peak in dp_peaks ]
        dp_score2, dp_rss2, dp_z2, dp_peaks2, S2, D2 = align( peaks, ladders, peak_pairs)
        if dp_score2 > dp_score:
            dp_score = dp_score2
            dp_rss = dp_rss2
            dp_z = dp_z2
            db_peaks = dp_peaks2

    return dp_score, dp_rss, dp_z, dp_peaks


def simple_align(peaks, ladders, peak_pairs):

    peak_pairs = z_align( peaks, ladders, peak_pairs )
    z, rss = dp.estimate_z( * zip( *peak_pairs ) )
    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp.align_peaks( ladders, peaks, z, rss )
    #print(' ==>', dp_score, dp_rss, len(dp_peaks))
    return (dp_score, dp_rss, dp_z, dp_peaks, S, D)

align = simple_align


def iter_align(peaks, ladders, peak_pairs, all_peaks):

    dp_results = []

    last_dp = last_rss = -1
    last_z = last_peaks = last_S = last_D = None
    while True:

        dp_score, dp_rss, dp_z, dp_peaks, S, D = simple_align(peaks, ladders, peak_pairs)
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



def greedy_1a( trace, peaks, avg_height, ladders, all_peaks, min_rss ):

    data = trace

    # make pairing of peaks & ladders

    delta = len(peaks) - len(ladders)
    if delta > 0:
        # there is possibility of having noise peaks
        I, J = 2 + delta, 3 + delta
    elif delta < 0:
        # there is possibility of losing ladder peaks
        I, J = 2, 3
        ladders = ladders[:delta]
    else:
        I, J = 2, 3

    dp_results = []

    for i in range(0, I):
        for j in range(1, J):


            # repairing peaks
            used_peaks = peaks[i:len(peaks) - j + 1]
            print('Peaks: %d => %d' % (len(peaks), len(used_peaks)))
            peak_pairs = []
            min_size = round( min( len(used_peaks), len(ladders) ) / 2 )
            for idx in range(min_size):
                peak_pairs.append( (used_peaks[idx].rtime, ladders[idx]) )
                peak_pairs.append( (used_peaks[-1-idx].rtime, ladders[-1-idx]) )
            peak_pairs.sort()
            #pprint.pprint( peak_pairs )

            z, rss = dp.estimate_z( [ x[1] for x in peak_pairs ],
                                    # => should adapt to N of peaks
                                    [ x[0] for x in peak_pairs ],
                                    3 )

            #peak_pairs = estimate_peak_pairs( peaks, ladders, z )
            #print('Initial peak pairs ==>')
            #pprint.pprint(peak_pairs)

            last_rss = -1
            stage = 1

            while True:

                # optimize locally just by using simple peak assignment
            
                z, rss = dp.estimate_z(    [ x[1] for x in peak_pairs ],
                                        [ x[0] for x in peak_pairs ] )

                peak_pairs = estimate_peak_pairs( peaks, ladders, z )
                print('Stage %d ==>' % stage, rss)
                #pprint.pprint(peak_pairs)

                if last_rss < 0:
                    last_rss = rss
                elif last_rss - rss <= 1:
                    break

                last_rss = rss
                stage += 1

            peak_pairs.sort()
            z, rss = dp.estimate_z( * zip( *peak_pairs ) )
            #pprint.pprint(peak_pairs)


            # use Dynamic Programming to get optimal RSS and DP score
            dp_score, dp_rss, dp_z, dp_peaks, S, D = dp.align_peaks( ladders, peaks, z, rss )
            print('==>', dp_score, dp_rss, len(dp_peaks))
            #pprint.pprint(dp_peaks)

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    #for dp_result in dp_results:
    #    print(' =>', dp_result[0], dp_result[1], len(dp_result[3]) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]

    return dp_score, dp_rss, dp_z, dp_peaks


def z_align( peaks, ladders, peak_pairs ):
    """ based on Z generated by peak pairs, realign peak & ladder """

    # peak_pairs is [ (rtime, size), ... ]

    last_rss = -1
    stage = 1

    while True:

        # optimize locally just by using simple peak assignment
            
        z, rss = dp.estimate_z(    [ x[1] for x in peak_pairs ],
                                        [ x[0] for x in peak_pairs ] )

        peak_pairs = estimate_peak_pairs( peaks, ladders, z )
        #print('Stage %d ==>' % stage, rss)
        #pprint.pprint(peak_pairs)

        if last_rss < 0:
            last_rss = rss
        elif last_rss - rss <= 1:
            break

        last_rss = rss
        stage += 1

    return peak_pairs


def shift_peak_pairs( peak_pairs, z ):
    """ shift a peak pairs based on the biggest deviation, return [ (rtime, size), ... ] """

    f = np.poly1d(z)
    max_i = 0
    max_value = 0

    pprint.pprint(peak_pairs)

    for (idx, peak_pair) in enumerate(peak_pairs):
        deviation = abs( f( peak_pair[1].rtime ) - peak_pair[0] )
        if deviation > max_value:
            max_value = deviation
            max_i = idx

    new_peak_pairs = []
    if max_i > len(peak_pairs)/2:
        # on the right hand
        for i in range(0, max_i):
            #new_peak_pairs.append( peak_pairs[i] )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i][0]) )
        for i in range(max_i, len(peak_pairs)):
            #new_peak_pairs.append( (peak_pairs[i-1][0], peak_pairs[i][1]) )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i-1][0]) )

    else:
        for i in range(0, max_i):
            #new_peak_pairs.append( (peak_pairs[i+1][0], peak_pairs[i][1]) )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i+1][0]) )

        for i in range(max_i, len(peak_pairs)):
            #new_peak_pairs.append( peak_pairs[i] )
            new_peak_pairs.append( (peak_pairs[i][1].rtime, peak_pairs[i][0]) )

    pprint.pprint(new_peak_pairs)
    return new_peak_pairs


def shift_peak_pairs( peak_pairs, z ):

    f = np.poly1d(z)
    max_i = 0
    max_value = 0

    for (idx, peak_pair) in enumerate(peak_pairs):
        deviation = abs( f( peak_pair[1].rtime ) - peak_pair[0] )
        cerr('  %2d  %3d  %6d  %5.2f' % (idx, peak_pair[0], peak_pair[1].rtime, deviation))
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

    cerr('after shifting...')
    for (idx, peak_pair) in enumerate(new_peak_pairs):
        deviation = abs( f( peak_pair[0] ) - peak_pair[1] )
        cerr('  %2d  %3d  %6d  %5.2f' % (idx, peak_pair[1], peak_pair[0], deviation))

    return new_peak_pairs
        

def generate_initial_peak_pairs( peaks, ladders, start_pos, end_pos ):

            used_peaks = peaks[start_pos:len(peaks) - end_pos + 1]
            print('Peaks: %d => %d' % (len(peaks), len(used_peaks)))
            peak_pairs = []
            min_size = round( min( len(used_peaks), len(ladders) ) / 2 )
            for idx in range(min_size):
                peak_pairs.append( (used_peaks[idx].rtime, ladders[idx]) )
                peak_pairs.append( (used_peaks[-1-idx].rtime, ladders[-1-idx]) )
            peak_pairs.sort()

            return peak_pairs



def generate_peak_pairs( ladders, peaks, start_pos, end_pos ):
    """ return [ (ladder, peak), (ladder, peak), ... ] """

    used_peaks = peaks[start_pos:len(peaks) - end_pos + 1]
    #print('Used peaks:')
    #pprint.pprint(used_peaks)
    min_size = round(min( len(ladders), len(used_peaks) ) / 2)

    print('min_size = %d; start = %d; end = %d;' % (min_size, start_pos, end_pos))
    peak_pairs = []
    for i in range(min_size):
        peak_pairs.append( (ladders[i], used_peaks[i]) )
        peak_pairs.append( (ladders[-1-i], used_peaks[-1-i]) )

    peak_pairs.sort()
    #pprint.pprint(peak_pairs)
    return peak_pairs
        


def greedy_2( trace, peaks, avg_height, ladders, all_peaks, min_rss ):

    # XXX: when peak number is less than ladders, only use ladders - N

    print('Peaks: %d with avg_height: %4.2f' % (len(peaks), avg_height))
    #pprint.pprint(peaks)

    data = trace

    # when peak number == ladder number, just assume that every peak matches with each ladder
    delta = len(peaks) - len(ladders)
    if delta == 0:
        i_max, j_max = 2, 3
    elif delta > 0:
        i_max, j_max = 2 + delta, 3 + delta
    else:
        i_max, j_max = 2, 3
        ladders = ladders[:delta]

    dp_results = []
    for i in range(0, i_max):
        for j in range(1, j_max):
            
            # estimate z for degree = 3

            peak_pairs = generate_peak_pairs( ladders, peaks, i, -j )
            #pprint.pprint(peak_pairs)


            z, rss = dp.estimate_z(    [ x[0] for x in peak_pairs ],
                                    # => should adapt to N of peaks
                                    [ x[1].rtime for x in peak_pairs ],
                                    3 )

            peak_pairs = estimate_peak_pairs( peaks, ladders, z )
            #print('Initial peak pairs ==>')
            #pprint.pprint(peak_pairs)

            last_rss = -1
            stage = 1

            while False:
            
                z, rss = dp.estimate_z(    [ x[1] for x in peak_pairs ],
                                        [ x[0] for x in peak_pairs ] )

                peak_pairs = estimate_peak_pairs( peaks, ladders, z )
                print('Stage %d ==>' % stage, rss)
                #pprint.pprint(peak_pairs)

                if last_rss < 0:
                    last_rss = rss
                elif last_rss - rss <= 1:
                    break

                last_rss = rss
                stage += 1
   
            # estimate initial Z based on peak_assignment

            peak_pairs.sort()
            z, rss = dp.estimate_z( * zip( *peak_pairs ) )
            #pprint.pprint(peak_pairs)


            # iterate using Dynamic Programming to get best RSS and DP score
            dp_score, dp_rss, dp_z, dp_peaks, S, D = dp.align_peaks( ladders, peaks, z, rss )
            print('==>', dp_score, dp_rss, len(dp_peaks))
            #pprint.pprint(dp_peaks)

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    #for dp_result in dp_results:
    #    print(' =>', dp_result[0], dp_result[1], len(dp_result[3]) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]

    return dp_score, dp_rss, dp_z, dp_peaks


def estimate_peak_pairs( peaks, ladders, z ):

    # generate synthetic standard peaks
    standard_peaks = []
    f = np.poly1d(z)
    for ladder in ladders:
        standard_peaks.append( ( f(ladder), ladder ) )

    #pprint.pprint( standard_peaks )

    # align peaks, must be done twice

    peak_pairs = []

    for (rtime, ladder) in standard_peaks:
        min_cost = abs( rtime - peaks[0].rtime )
        aln_peak = peaks[0]
        for p in peaks[1:]:
            c = abs( rtime - p.rtime )
            if c < min_cost:
                min_cost = c
                aln_peak = p
        peak_pairs.append( (aln_peak.rtime, ladder) )

    peak_pairs.sort()
    return peak_pairs

greedy_1 = greedy_1b
