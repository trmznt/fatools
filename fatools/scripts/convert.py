

import sys
import argparse

from fatools.lib.traceio import read_abif_stream
from pprint import pprint

def init_argparser():

    parser = argparse.ArgumentParser( 'convert' )

    parser.add_argument('infiles', nargs='+')

    return parser


def main(args):

    for infile in args.infiles:
        with open( infile, 'rb' ) as instream:
            t = read_abif_stream( instream )
        channels = t.get_channels()
        names = [ '"' + c + '"' for c in channels ]
        print("Dyes: %s" % ' '.join( channels ))
        with open( infile + '.raw.tab', 'wt') as out:
            out.write( '\t'.join( names ) )
            out.write( '\n' )
            for p in zip( * [ channels[c].raw for c in channels ] ):
                out.write( '\t'.join( str(x) for x in p ) )
                out.write( '\n' )
        with open( infile + '.base.tab', 'wt') as out:
            out.write( '\t'.join( names ) )
            out.write( '\n' )
            for p in zip( * [ channels[c].smooth() for c in channels ] ):
                out.write( '\t'.join( str(x) for x in p ) )
                out.write( '\n' )


