
import sys

def cout(s, nl=True, flush=False):
    sys.stdout.write(s)
    if nl: sys.stdout.write('\n')
    if flush: sys.stdout.flush()

def cerr(s, nl=True, flush=False):
    sys.stderr.write(s)
    if nl: sys.stderr.write('\n')
    if flush: sys.stderr.flush()

def cexit(s, code=1):
    cerr(s)
    sys.exit(code)


def get_dbhandler(args, initial=False):
    """ return suitable handler from args """

    if args.sqldb:
        
        from fatools.lib.sqlmodels.handler import SQLHandler

        return SQLHandler(args.sqldb, initial)

    elif args.fsdb is not False:
        # we use filesystem-based database system
        raise NotImplementedError()


    cerr('ERR: Please specify database system using --sqldb or --fsdb options!')
    sys.exit(1)


def tokenize(options, converter=None):
    """ return { 'A': '1,2,3', 'B': True } for options 'A=1,2,3;B' """
    opt_dict = {}
    for token in options.split(';'):
        keys = token.split('=',1)
        if len(keys) == 1:
            opt_dict[keys[0].strip()] = True
        else:
            opt_dict[keys[0].strip()] = keys[1].strip()

    return opt_dict

    
