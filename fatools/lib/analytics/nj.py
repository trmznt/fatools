# Neighbor-joining tree from distance matrix

from subprocess import call
from fatools.lib.utils import random_string


def plot_nj( distance_matrix, tmp_dir, fmt='pdf', label_callback=None ):
    """ NJ uses R's ape library
        R will be called as a separate process instead as a embedded library
        in order to utilize paralel processing in multiple processor
    """

    # create matrix distance & color distance

    file_id = random_string(3)

    matrix_file = '%s/matrix-distance-%s.txt' % (tmp_dir, file_id)
    colors_file = '%s/colors-distance-%s.txt' % (tmp_dir, file_id)
    script_file = '%s/njtree-%s.r' % (tmp_dir, file_id)
    njtree_file = '%s/njtree-%s.%s' % (tmp_dir, file_id, fmt)

    with open(matrix_file, 'w') as out_m, open(colors_file, 'w') as out_c:

        out_m.write( '\t'.join( str(x) for x in distance_matrix.sample_ids ) )
        out_m.write( '\n')
        for name, vals in zip( distance_matrix.sample_ids, distance_matrix.M ):
            out_m.write( '%s\t%s\n' % (str(name), '\t'.join( ['%2.3f' % x for x in vals] ) ))

        out_c.write('\n'.join( distance_matrix.C ) )

    with open(script_file, 'w') as scriptout:
        if fmt == 'pdf':
            cmd = 'pdf("%s", width = 11.2, height=7)' % njtree_file
        elif fmt == 'png':
            cmd = 'png("%s", width = 1024, height = 640)' % njtree_file
        scriptout.write("""
library(ape)
M <- as.matrix( read.table("%s", sep='\\t', header=T) )
C <- as.vector( read.table("%s", sep='\\t', header=F, comment.char = '')[,1] )
tree <- nj( M )
%s
plot(tree, "fan", tip.color = C, font=1, cex=0.7, label.offset = 0.009)
legend('topright', inset=c(0,0), c(%s), col = c(%s), lty=1, cex=0.85, xpd=T)
""" % (matrix_file, colors_file, cmd,
        ",".join( '"%s"' % hs.label for (hs,_,_) in distance_matrix.S),
        ",".join( '"%s"' % hs.colour for (hs,_,_) in distance_matrix.S) )
    )

    ok = call( ['Rscript', script_file] )

    if ok != 0:
        raise RuntimeError("Rscript run unsucessfully")

    return njtree_file