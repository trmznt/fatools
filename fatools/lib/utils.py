
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


_VERBOSITY_ = 0
def set_verbosity(value):
    global _VERBOSITY_
    _VERBOSITY_ = value


def is_verbosity(value):
    global _VERBOSITY_
    return (_VERBOSITY_ >= value)

def cverr(value, txt, nl=True, flush=False):
    global _VERBOSITY_
    if _VERBOSITY_ >= value:
        cerr(txt, nl, flush)

def get_dbhandler(args, initial=False):
    """ return suitable handler from args """

    if args.sqldb:

        if args.schema == 1:

            from fatools.lib.sqlmodels.handler import SQLHandler

            return SQLHandler(args.sqldb, initial)

        elif args.schema == 2:

            from fatools.lib.sqlmodels2.handler import SQLHandler

            return SQLHandler(args.sqldb, initial)

        else:
            cexit('Please provide correct sql schema version!')

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


# utility to deal with tab or comma delimited text buffer

def detect_buffer( buf ):
    """ return (buf, delimiter) """

    # find our new line character, this is for Mac Excel blunder

    n_count = buf.count('\n')
    r_count = buf.count('\r')

    if n_count == 0 and r_count > 0:
        # Mac Excel
        buf = buf.replace('\r', '\n')
        n_count = r_count
    elif r_count > n_count:
        raise RuntimeError('Invalid text content')

    # we detect delimiter
    tab_count = buf.count('\t')
    comma_count = buf.count(',')

    if comma_count > tab_count and comma_count > n_count:
        return (buf, ',')
    return (buf, '\t')
