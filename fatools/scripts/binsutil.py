
import sys, argparse

def init_argparser(parser=None):

    if parser:
        p = parser
    else:
        p = argparse.ArgumentParser('binsutil')

    p.add_argument('--infile', required=True)

    p.add_argument('--outfile')

    p.add_argument('--marker')

    p.add_argument('--repeats', type=int)

    p.add_argument('--min', type=int)

    p.add_argument('--max', type=int)

    p.add_argument('--anchor', type=int)

    return p


def main(args):

    do_binsutil(args)


def do_binsutil(args):
    from fatools.lib.fautil import binsutil
    binsutil.do_binsutil(args)
