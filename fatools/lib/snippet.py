def adaptive_align_naive( trace, peaks, avg_height, ladders ):

    ## TRY 2 -> use simple linear regression

    data = trace

    from matplotlib import pylab as plt

    linear_alignments = []
    for i in range(0, 3):
        for j in range(1, 3):
            
            # estimate z for degree = 1
            z, rss = estimate_z( [ ladders[0], ladders[-1] ], # => should adapt to N of peaks
                                [ peaks[i].rtime, peaks[-j].rtime ],
                                 1 )

            # generate standard peaks
            standard_peaks = []
            f = np.poly1d(z)
            for ladder in ladders:
                standard_peaks.append( ( f(ladder), ladder ) )

            pprint.pprint( standard_peaks )

            # calculate costs
            cost = 0.0
            for (rtime, _) in standard_peaks:
                min_cost = abs( rtime - peaks[0].rtime )
                for p in peaks[1:]:
                    c = abs( rtime - p.rtime )
                    if c < min_cost:
                        min_cost = c
                cost += min_cost

            linear_alignments.append( (cost, z, rss) )

            plt.plot(data)
            for p in standard_peaks:
                plt.plot( [ p[0], p[0] ], [ 0, avg_height ], 'g' )

            plt.show()


    linear_alignments.sort()

    pprint.pprint(linear_alignments)
    # iterate using Dynamic Programming to get best RSS and DP score
    for alignment in linear_alignments:
        (_, z, rss) = alignment
        dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders, peaks, z, rss )
        print(" => ",dp_score, dp_rss, len(dp_peaks))

    pprint.pprint(dp_peaks)

    for (std_size, peak) in dp_peaks:
        peak.size = std_size
        peak.type = peaktype.ladder

    return dp_score, dp_rss, dp_z, dp_peaks


def adaptive_align_naive_2( trace, peaks, avg_height, ladders ):

    # XXX: need to think when peak number < ladders (ie. missing ladder peaks) !!

    data = trace

    from matplotlib import pylab as plt

    pprint.pprint(peaks)

    dp_results = []
    for i in range(0, 3):
        for j in range(1, 3):
            
            # estimate z for degree = 1
            z, rss = estimate_z(    [ ladders[0], ladders[1], ladders[-2], ladders[-1] ],
                                    # => should adapt to N of peaks
                                    [ peaks[i].rtime, peaks[i+1].rtime,
                                      peaks[-j-1].rtime, peaks[-j].rtime ],
                                    3 )

            peak_pairs = estimate_peak_pairs( peaks, ladders, z )
            print('Initial peak pairs ==>')
            pprint.pprint(peak_pairs)

            last_rss = -1
            stage = 1

            while True:
            
                z, rss = estimate_z(    [ x[1] for x in peak_pairs ],
                                        [ x[0] for x in peak_pairs ] )

                peak_pairs = estimate_peak_pairs( peaks, ladders, z )
                print('Stage %d ==>' % stage, rss)
                pprint.pprint(peak_pairs)

                if last_rss < 0:
                    last_rss = rss
                elif last_rss - rss <= 1:
                    break

                last_rss = rss
                stage += 1
   
            # estimate initial Z based on peak_assignment

            peak_pairs.sort()
            z, rss = estimate_z( * zip( *peak_pairs ) )
            pprint.pprint(peak_pairs)


            # iterate using Dynamic Programming to get best RSS and DP score
            dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders, peaks, z, rss )
            print('==>', dp_score, dp_rss, len(dp_peaks))
            #pprint.pprint(dp_peaks)

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    #for dp_result in dp_results:
    #    print(' =>', dp_result[0], dp_result[1], len(dp_result[3]) )

    dp_results.sort( reverse=True )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]

    for (std_size, peak) in dp_peaks:
        peak.size = std_size
        peak.type = peaktype.ladder

    return dp_score, dp_rss, dp_z, dp_peaks


def adaptive_align_naive_3( trace, peaks, avg_height, ladders ):

    # XXX: when peak number is less than ladders, only use ladders - N

    data = trace

    dp_results = []
    for i in range(0, 2):
        for j in range(1, 3):
            
            # estimate z for degree = 3

            peak_pairs1 = list( zip( ladders, peaks[i:] ) )
            peak_pairs2 = list( zip( reversed(ladders), reversed(peaks[:-j]) ) )
            peak_pairs = list( sorted( peak_pairs1 + peak_pairs2,
                                key = lambda x: (x[0], x[1].rtime) ) )


            z, rss = estimate_z(    [ x[0] for x in peak_pairs ],
                                    # => should adapt to N of peaks
                                    [ x[1].rtime for x in peak_pairs ],
                                    3 )

            peak_pairs = estimate_peak_pairs( peaks, ladders, z )
            #print('Initial peak pairs ==>')
            #pprint.pprint(peak_pairs)

            last_rss = -1
            stage = 1

            while True:
            
                z, rss = estimate_z(    [ x[1] for x in peak_pairs ],
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
   
            # estimate initial Z based on peak_assignment

            peak_pairs.sort()
            z, rss = estimate_z( * zip( *peak_pairs ) )
            #pprint.pprint(peak_pairs)


            # iterate using Dynamic Programming to get best RSS and DP score
            dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders, peaks, z, rss )
            print('==>', dp_score, dp_rss, len(dp_peaks))
            #pprint.pprint(dp_peaks)

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    #for dp_result in dp_results:
    #    print(' =>', dp_result[0], dp_result[1], len(dp_result[3]) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]


    for (std_size, peak) in dp_peaks:
        peak.size = std_size
        peak.type = peaktype.ladder

    return dp_score, dp_rss, dp_z, dp_peaks


def adaptive_align_dtw( trace, peaks, avg_height, ladders ):

    data = trace

    ## TRY 1 -> use DTW

    from mlpy import dtw_std
    from matplotlib import pylab as plt
    from dtw import dtw

    peak_corr = {}
    for p in peaks:
        peak_corr[p.rtime] = []

    dtw_list = []
    for i in range(0, 3):
        for j in range(1, 3):

            standard_peaks, peak_index = generate_peaks( ladders, avg_height,
                                                peaks[i].rtime, peaks[-j].rtime )

            dist, cost, path = dtw_std( standard_peaks, data, dist_only=False, squared=True)
            plot_path( standard_peaks, data, path, [ p[0] for p in peak_index ] )
            #dist, cost, path = dtw( standard_peaks, data )

            # fill peak correlation based on path

            for map_x, map_y in zip( path[0], path[1] ):
                if map_y in peak_corr:
                    if standard_peaks[map_x] < avg_height/2:
                        peak_corr[map_y].append( -1 )
                    else:
                        peak_corr[map_y].append( search_peak_index( map_x, peak_index ) )

    peak_assignment = score_peak_correlation( peak_corr )
    dpscore, rss, z, aligned_peaks = adaptive_peak_alignment( peak_assignment, peaks, ladders )

    return (dpscore, rss, z, aligned_peaks)



def adaptive_align_naive_4( trace, peaks, avg_height, ladders ):

    # XXX: when peak number is less than ladders, only use ladders - N

    # when peak number == ladder number, just assume that every peak matches with each ladder
    if len(peaks) == len(ladders):
        i_max = 1
    else:
        i_max = 3

    data = trace
    reversed_peaks = list(reversed(peaks))
    reversed_ladders = list(reversed(ladders))


    dp_results = []
    for i in range(0, i_max):
            
            # estimate z for degree = 3


            peak_pairs = list( sorted( zip( reversed_ladders, reversed_peaks[i:] ),
                                key = lambda x: (x[0], x[1].rtime) ) )


            z, rss = estimate_z(    [ x[0] for x in peak_pairs ],
                                    # => should adapt to N of peaks
                                    [ x[1].rtime for x in peak_pairs ],
                                    3 )

            peak_pairs = estimate_peak_pairs( peaks, ladders, z )
            #print('Initial peak pairs ==>')
            #pprint.pprint(peak_pairs)

            last_rss = -1
            stage = 1

            while True:
            
                z, rss = estimate_z(    [ x[1] for x in peak_pairs ],
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
            z, rss = estimate_z( * zip( *peak_pairs ) )
            #pprint.pprint(peak_pairs)


            # iterate using Dynamic Programming to get best RSS and DP score
            dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders, peaks, z, rss )
            print('==>', dp_score, dp_rss, len(dp_peaks))
            #pprint.pprint(dp_peaks)

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    #for dp_result in dp_results:
    #    print(' =>', dp_result[0], dp_result[1], len(dp_result[3]) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]


    for (std_size, peak) in dp_peaks:
        peak.size = std_size
        peak.type = peaktype.ladder

    return dp_score, dp_rss, dp_z, dp_peaks


def generate_peak_pairs( ladders, peaks, start_pos, end_pos ):
    pass

def adaptive_align_naive_5( trace, peaks, avg_height, ladders ):
    # similar to align_naive_4, but with initial peak_pairs taken from the first and last of
    # peaks

    # XXX: when peak number is less than ladders, only use ladders - N

    data = trace

    # when peak number == ladder number, just assume that every peak matches with each ladder
    if len(peaks) == len(ladders):
        i_max, j_max = 1, 2
    else:
        i_max, j_max = 2, 3

    dp_results = []
    for i in range(0, 2):
        for j in range(1, 3):
            
            # estimate z for degree = 3

            peak_pairs = list( sorted( zip( ladders, peaks[i:] ),
                                key = lambda x: (x[0], x[1].rtime) ) )


            z, rss = estimate_z(    [ x[0] for x in peak_pairs ],
                                    # => should adapt to N of peaks
                                    [ x[1].rtime for x in peak_pairs ],
                                    3 )

            peak_pairs = estimate_peak_pairs( peaks, ladders, z )
            #print('Initial peak pairs ==>')
            #pprint.pprint(peak_pairs)

            last_rss = -1
            stage = 1

            while True:
            
                z, rss = estimate_z(    [ x[1] for x in peak_pairs ],
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
   
            # estimate initial Z based on peak_assignment

            peak_pairs.sort()
            z, rss = estimate_z( * zip( *peak_pairs ) )
            #pprint.pprint(peak_pairs)


            # iterate using Dynamic Programming to get best RSS and DP score
            dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders, peaks, z, rss )
            print('==>', dp_score, dp_rss, len(dp_peaks))
            #pprint.pprint(dp_peaks)

            dp_results.append( (dp_score, dp_rss, dp_z, dp_peaks, S, D) )

    #for dp_result in dp_results:
    #    print(' =>', dp_result[0], dp_result[1], len(dp_result[3]) )

    dp_results.sort( reverse=True, key = lambda x: (x[0], x[1]) )

    dp_score, dp_rss, dp_z, dp_peaks, S, D = dp_results[0]


    for (std_size, peak) in dp_peaks:
        peak.size = std_size
        peak.type = peaktype.ladder

    return dp_score, dp_rss, dp_z, dp_peaks


def generate_peaks( peaks, height, min_range, max_range ):

    longest_peak = max(peaks)
    shortest_peak = min(peaks)
    ratio = (max_range - min_range) / (longest_peak - shortest_peak)

    offset = round(min_range - shortest_peak * ratio)
    N = max_range + 50 * ratio
    x = np.linspace(0, N, N)
    print(N, ratio)
    
    from matplotlib import pylab as plt

    arrays = []
    peak_index = []
    for p in peaks:
        idx = round(p * ratio) + offset
        peak_index.append( (idx, p) )
        a = peak_func(x, height, idx, 1)
        arrays.append(a)

    gram = sum(arrays)
    #plt.plot(p)
    #plt.show()

    return gram, peak_index

    from matplotlib import pylab as plt
    import sys

    plt.plot(p)
    plt.show()
    sys.exit(1)


def peak_func(x, a, x0, sigma):
    return a*np.exp(-(x-x0)**2/(2*sigma**2))


def search_peak_index( idx, peak_list ):

    p = peak_list[0][1]
    d = abs( idx - peak_list[0][0] )

    for pos, std in peak_list[1:]:
        _d = abs( idx - pos )
        if _d < d:
            p = std
            d = _d

    return p


def adaptive_peak_alignment( peak_assignment, peaks, ladders ):
    """ return dpscore, rss, z, aligned_peaks """

    # estimate initial Z based on peak_assignment
    peak_pairs = []
    pprint.pprint( peak_assignment )
    for rtime, val in peak_assignment.items():
        if val[1] < 0.9: continue
        peak_pairs.append( (rtime, val[0]) )
    peak_pairs.sort()
    z, rss = estimate_z( * zip( *peak_pairs ) )

    # iterate using Dynamic Programming to get best RSS and DP score
    dp_score, dp_rss, dp_z, dp_peaks, S, D = align_peaks( ladders, peaks, z, rss )

    pprint.pprint(dp_peaks)

    for (std_size, peak) in dp_peaks:
        peak.size = std_size
        peak.type = peaktype.ladder

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

    return peak_pairs


## FUNCTIONS BELOW ARE USEFUL FOR ALGORITHM DEBUGGING

def plot_path( standard_peaks, data, path, peaks ):

    from matplotlib import pylab as plt

    plt.plot(data, 'r')
    plt.plot(standard_peaks, 'y')

    #print(path)
    #rtimes = [ x.rtime for x in peaks ]
    for [map_x, map_y] in zip( path[0], path[1] ):
        if map_x in peaks:
            plt.plot( [map_x, map_y], [ standard_peaks[map_x], data[map_y] ], 'g' )
    plt.show()


def simple_pca( a1 ):

    M = np.zeros( (len(a1), len(a1)) )

    for i in range( len(a1) ):
        for j in range( i, len(a1) ):
            M[i,j] = M[j,i] = (a1[i]-a1[j]) ** 2

    import mdp

    return mdp.pca( M, output_dim = 2 )



def plot_pca( comps, labels ):

    from matplotlib import pylab as plt

    plt.scatter( comps[:,0], comps[:,1] )
    #plt.show()
    #return
    
    ytext = -50
    for label, x, y in zip(labels, comps[:, 0], comps[:, 1]):
        plt.annotate(
            label, 
            xy = (x, y), xytext = (0, ytext),
            textcoords = 'offset points', ha = 'right', va = 'bottom',
            bbox = dict(boxstyle = 'round,pad=0.5', fc = 'yellow', alpha = 0.5),
            arrowprops = dict(arrowstyle = '->', connectionstyle = 'arc3,rad=0'))
        ytext = ytext * -1

    plt.show()

