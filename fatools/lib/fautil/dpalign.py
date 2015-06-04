
import numpy as np
import math
import pprint

from matplotlib import pylab as plt

# dynamic programming for peak alignment
#
# dp: y -> sizes
#     x -> retention_time
#


def estimate_z( x, y, degree = 3 ):
    """ estimate z and rss
            x: peak rtime
            y: standard sizes
        return (z, rss)

        y ~ f(x) where f = poly1d(z)
        rss ~ SUM( (f(x) - y)**2 ) for all (x,y)
    """
    z = np.polyfit( x, y, degree )
    p = np.poly1d( z )
    y_p = p(x)
    rss = ( (y_p - y) ** 2 ).sum()

    return z, rss


def annotate(M, ladders, peaks):
    M_ = np.insert( M, 0, values = np.array( [x.rtime for x in peaks] ), axis = 0 )
    return np.insert(M_, 0, values = np.array( [0] + ladders ), axis = 1 )


def align_peaks( ladders, peaks, z, rss ):
    """ align ladders with peaks using dynamic programming (global alignment)
        return (dpscore, RSS, Z, ladder_aligned_peaks)
    """

    # set initial peak alignment, just align from the biggest peaks
    ladders = list(sorted(ladders, reverse=True))
    peaks = list(sorted( peaks, key = lambda x: x.rtime, reverse = True))

    #print('initial RSS: %4.2f' % rss)
    dpscore = -1

    while True:

        S = generate_scores( ladders, peaks, np.poly1d(z))

        #pprint.pprint(S)

        result = dp(S, -5e-3)

        #pprint.pprint(result)

        cur_dpscore = result['D'][-1][-1]
        matches = result['matches']

        aligned_peaks = [ (ladders[i], peaks[j]) for i, j in matches ]

        # realign

        std_size, peak_sizes = zip(*aligned_peaks)
        cur_z, cur_rss = estimate_z( [x.rtime for x in peak_sizes], std_size )
        #print('current DP score: %3.3f RSS: %4.2f' % (cur_dpscore, cur_rss))
        #if cur_rss > rss:
        if cur_dpscore < dpscore:
            print('WARNING: algorithm did not converge!!')
            break

        #if cur_rss == rss:
        if cur_dpscore == dpscore:
            break

        z = cur_z
        rss = cur_rss
        dpscore = cur_dpscore
        sized_peaks = aligned_peaks

    return dpscore, rss, z, sized_peaks, annotate(S, ladders, peaks), result['D']


def create_scoring_function(A, B):
    def _scoring(ladder,peak,height):
        return ((ladder * A + B) - peak) * 1/height
    return _scoring


def generate_scores(ladders, peaks, func, tolerance = 4):
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
    M = np.zeros( (len(ladders), len(peaks)), dtype='d' ) 

    _TOL = 0.001
    cutoff = tolerance * math.sqrt( -1.0 * math.log(_TOL))
    highest_peak = max( [ x.height for x in peaks ] )
    ladder_N = len(ladders)

    row = 0
    col = 0

    for ladder in ladders:
        for peak in peaks:
            size = func(peak.rtime)
            similarity = (np.log10( peak.height/ highest_peak ) + ladder_N) / ladder_N
            M[row][col] = similarity * math.exp( - ((size - ladder)/(tolerance))**2  / 2 )
            #print(ladder, peak)
            #print('  =>', '%3.2f' % size, '%1.2f' % similarity, M[row][col] )
            col += 1
        row += 1
        col = 0

    return M


def plot_z(peaks, ladders, z):

    
    #x = np.linspace(0, peaks[-1].rtime + 100)
    
    x = np.linspace(0, max(ladders))
    p = np.poly1d( z )
    y_p = p(x)

    plt.plot(x, y_p)
    plt.show()


def plot_S(S, path=None):

    from matplotlib import pylab as plt

    im1 = plt.imshow(S, interpolation='nearest', cmap='Reds')
    #im2 = plt.plot(path[0], path[1], 'k')
    plt.gca().invert_yaxis()
    plt.xlabel("STD")
    plt.ylabel("PEAK")
    plt.grid()
    plt.colorbar()
    plt.show()


"""
Dynamic Programming routine
"""

 #############################################################################
 #                                                                           #
 #    PyMS software for processing of metabolomic mass-spectrometry data     #
 #    Copyright (C) 2005-2011 Vladimir Likic                                 #
 #                                                                           #
 #    This program is free software; you can redistribute it and/or modify   #
 #    it under the terms of the GNU General Public License version 2 as      #
 #    published by the Free Software Foundation.                             #
 #                                                                           #
 #    This program is distributed in the hope that it will be useful,        #
 #    but WITHOUT ANY WARRANTY; without even the implied warranty of         #
 #    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the          #
 #    GNU General Public License for more details.                           #
 #                                                                           #
 #    You should have received a copy of the GNU General Public License      #
 #    along with this program; if not, write to the Free Software            #
 #    Foundation, Inc., 675 Mass Ave, Cambridge, MA 02139, USA.              #
 #                                                                           #
 #############################################################################

import numpy

def dp(S, gap_penalty, peak_penalty = 0):
   
    """
    @summary: Solves optimal path in score matrix based on global sequence
    alignment

    @param S: Score matrix
    @type S: numpy.
    @param gap_penalty: Gap penalty
    @type gap_penalty: FloatType

    @return: A dictionary of results
    @rtype: DictType

    @author: Tim Erwin
    """
   
    row_length = len(S[:,0])
    col_length = len(S[0,:])

    # free for 25% peaks
    missing_20 = row_length * 0.25

    #D contains the score of the optimal alignment
    D = numpy.zeros((row_length+1,col_length+1), dtype='d')
    for i in range(1, row_length+1):
        #D[i,0] = 0 if i < missing_20 else gap_penalty*(i-missing_20)
        # missing ladders
        D[i,0] = 0.25 * gap_penalty
    for j in range(1, col_length+1):
        # missing peaks
        D[0,j] = 0 #gap_penalty*j
    D[0,0] = 0.0
    D[1:(row_length+1), 1:(col_length+1)] = S.copy();

    # Directions for trace
    # 0 - match               (move diagonal)
    # 1 - peaks1 has no match (move up)
    # 2 - peaks2 has no match (move left)
    # 3 - stop
    trace_matrix = numpy.zeros((row_length+1,col_length+1))
    trace_matrix[:,0] = 1;
    trace_matrix[0,:] = 2;
    trace_matrix[0,0] = 3;
   
    for i in range(1,row_length+1):
        for j in range(1,col_length+1):
 
            #
            # Needleman-Wunsch Algorithm assuming a score function S(x,x)=0
            #
            #              | D[i-1,j-1] + S(i,j)
            # D[i,j] = min | D(i-1,j] + gap
            #              | D[i,j-1] + gap
            #

            #darray = [D[i-1,j-1]+S[i-1,j-1], D[i-1,j]+gap_penalty, D[i,j-1]+gap_penalty]
            if j == col_length:
                penalty = 0.25 * gap_penalty
            else:
                penalty = gap_penalty
            darray = [D[i-1,j-1]+S[i-1,j-1], D[i-1,j]+penalty, D[i,j-1] + peak_penalty]
            D[i,j] = max(darray)
            #Store direction in trace matrix
            trace_matrix[i,j] = darray.index(D[i,j])

    # Trace back from bottom right
    trace = []
    matches = []
    i = row_length
    j = col_length
    direction = trace_matrix[i,j]
    p = [row_length-1]
    q = [col_length-1]
   
    while direction != 3:
       
        if direction == 0: #Match
            i = i-1
            j = j-1
            matches.append([i,j])
        elif direction == 1: #peaks1 has no match
            i = i-1
        elif direction == 2: #peaks2 has no match
            j = j-1
        p.append(i-1)
        q.append(j-1)
        trace.append(direction)
        direction=trace_matrix[i,j]

    #remove 'stop' entry
    p.pop()
    q.pop()
    # reverse the trace back
    p.reverse()
    q.reverse()
    trace.reverse()
    matches.reverse()

    return {'p':p, 'q':q, 'trace':trace, 'matches':matches, 'D':D, 'phi':trace_matrix}


