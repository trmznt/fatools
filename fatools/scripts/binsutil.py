
import sys, argparse

# get the helpers

from fatools.scripts.dbmgr import get_assay_list

def init_argparser(parser=None):

    if parser:
        p = parser
    else:
        p = argparse.ArgumentParser('binsutil')

    p.add_argument('--sqldb', default=False,
        help = 'SQLite3 database filename')

    p.add_argument('--batch', default=False,
        help = 'Batch code')

    p.add_argument('--infile')

    p.add_argument('--outfile')

    p.add_argument('--marker')

    p.add_argument('--commit', default=False, action='store_true',
        help = 'commit to database')

    # commands

    p.add_argument('--init', default=False,
        help = 'create bins for a particular marker / batch')

    p.add_argument('--show', default=False,
        help = 'show bins for a particular marker / batch')

    p.add_argument('--optimize', default=False, action='store_true',
        help = 'optimize bins for a particular marker / batch')

    p.add_argument('--summarize', default=False, action='store_true',
        help = 'summarize data frame')

    p.add_argument('--repeats', type=int)

    p.add_argument('--min', type=int)

    p.add_argument('--max', type=int)

    p.add_argument('--anchor', type=int)

    p.add_argument('--shift', type=float, default=0)



    return p


def main(args):

    do_binsutil(args)


def do_binsutil(args, dbh=None):
    from fatools.lib.fautil import binsutil
    binsutil.do_binsutil(args)
