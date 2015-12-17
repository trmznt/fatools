
import sys, base64, os, math

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


random_string = lambda n: base64.b64encode(os.urandom(int(math.ceil(0.75*n))), b'-_')[:n].decode('UTF-8')


_R_lock_ = None

def acquire_R():

    global _R_lock_
    if _R_lock_ is None:

        # initialize rpy2 and set thread lock

        from rpy2 import robjects
        from rpy2.robjects import pandas2ri
        import threading

        pandas2ri.activate()
        _R_lock_ = threading.Lock()

    _R_lock_.acquire()


def release_R():

    global _R_lock_

    _R_lock_.release()
