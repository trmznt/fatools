
import sys, os
import argparse
import importlib

def greet():
    print('fatools - Python-based DNA fragment-analysis tools')

def usage():
    print('Usage:')
    print('\t%s command [options]' % sys.argv[0])
    sys.exit(0)

def main():

    greet()

    command = sys.argv[1]
    opt_args = sys.argv[2:]

    print('Running command: %s' % command)

    try:
        M = importlib.import_module('fatools.scripts.' + command)
    except ImportError:
        print('Cannot import script name: %s' % command)
        raise

    parser = M.init_argparser()
    args = parser.parse_args(opt_args)
    M.main( args )



