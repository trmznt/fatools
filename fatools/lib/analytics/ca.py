#
# Component Analysis
#
# PCoA (Principal Coordinate Analysis)
#   PCA with distance matrix, using Python MDP library
# MCA (Multiple Correspondence Analysis)
#   PCA with categorical data, using R FactoMineR library
#

import numpy as np
import matplotlib.pyplot as plt
import mdp

def jitters( data ):
    """ returning normal distribution noise for data jittering"""
    # need to calculate s based on the smallest and longest distance of the points
    d_min = d_max = np.abs(data[0] - data[1])
    for i in range(len(data)):
        for j in range(len(data)):
            d = np.abs( data[i] - data[j] )
            if d > 0:
                if d > d_max:
                    d_max = d
                elif d < d_min:
                    d_min = d

    s = d_min / d_max * len(data) / 2
    print("Jitters:", s)
    return s * np.random.randn( len(data) )


def pcoa( distance_matrix, dim = 2 ):

    pcan = mdp.nodes.PCANode( output_dim = dim )
    pcar = pcan.execute( distance_matrix.M )

    for i in range(dim):
        pcar[:, i] += jitters( pcar[:, i])

    return (pcar, pcan.d)


def plot_pca( pca_result, distance_matrix, pc1, pc2, filename=None ):
    """ plot PCA result using matplotlib """

    if not filename:
        raise RuntimeError('plot_pca() needs filename argument')

    fig = plt.figure()
    ax = fig.add_subplot(111)

    pca_matrix = pca_result[0]
    pca_var = pca_result[1]

    for hs, s, e in distance_matrix.S:
        ax.scatter( pca_matrix[s:s+e, pc1],
                    pca_matrix[s:s+e, pc2],
                    c = hs.colour,
                    edgecolor = hs.colour,
                    label = hs.label,
                    alpha = 0.75,
                    marker='+',
                    s = 30
        )

    if pca_var is not None:
        ax.set_xlabel('PC%d (%.3f%%)' % (pc1 + 1, pca_var[pc1]))
        ax.set_ylabel('PC%d (%.3f%%)' % (pc2 + 1, pca_var[pc2]))
    else:
        ax.set_xlabel('PC%d' % (pc1 + 1))
        ax.set_ylabel('PC%d' % (pc2 + 1))
    leg = ax.legend(loc='upper left', scatterpoints=1, fontsize='x-small', fancybox=True,
            bbox_to_anchor=(1,1))
    #leg.get_frame().set_alpha(0.5)
    fig.savefig( filename, bbox_extra_artists=(leg,), bbox_inches='tight' )
    plt.close()
    return filename


def mca( distance_matrix, dim = 2 ):
    """ calculate MCA matrix using R's FactorMineR """

    # build up haplotype dataframe

    from rpy2 import robjects
    from rpy2.robjects import pandas2ri

    pandas2ri.activate()

    r_df = pandas2ri.py2ri(distance_matrix.H)
    robjects.globalenv['haplo_data'] = r_df
    marker_len = len(distance_matrix.H.columns)
    arguments = ','.join('as.factor(haplo_data[,%d])' % x
                    for x in range(1, marker_len + 1))
    robjects.r('haplo_df <- data.frame(%s)' % arguments)
    robjects.r('library(FactoMineR)')
    mca_res = robjects.r('MCA(haplo_df, graph=FALSE)')

    # get the individual coordinate
    coord = pandas2ri.ri2py(mca_res.rx('ind')[0].rx('coord')[0])

    return (coord, None)

    #mca_result = mca_func(r_df, ** { 'quanti.sup': [1,2,3] })
    #

    raise RuntimeError


def pca( distance_matrix, dim = 2 ):
    """ calculate PCA matrix using in-house algorithm, treating variables
        as categorical data
    """

    D = {}

    for c in distance_matrix.H.columns:
        D[c] = pandas.get_dummies( distance_matrix.H[c])

    #H = pandas.concat()

    # this function hasn't been finished yet
    raise NotImplementedError
