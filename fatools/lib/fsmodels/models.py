# basic classes that mimic msdb classes
# this is filesystem-based database

import .fautils2 as fautils

from .mixin import *



class sessionmgr(object):

    def __init__(self):
        self.rootdir = None

    def set_rootdir(self, rootdir):
        if self.rootdir is not None:
            raise RuntimeError('db is alredy open for rootdir: %s' % rootdir)

        self.rootdir = rootdir

    def get_rootdir(self):
        if self.rootdir is None:
            raise RuntimeError('db has not been set')


dbsession = sessionmgr()


class base(object):

    def __init__(self):
        pass

    def query(self):
        


class Sample(base):

    def __init__(self):
        self.assays = None


class Assay(base, AssayMixIn):

    def __init__(self):
        self.channels = None
        self.rss = -1
        self.dp = -1
        self.score = -1

        self.sample = None
        self.panel = None


    def set_directory(self, dirname):
        self._dirname = dirname
        self._assayname = None


    def load(self, with_trace=True):
        """ load assay from files into memory """

        # try to load channels from yaml if exists
        channel_yaml = "%s/channels.yaml" % self._dirname
        if os.path.exists( channel_yaml ):
            # load yaml
            pass


        # load smooth data, if not exists
        trace_file = "%s/traces.tab" % self._dirname
        if with_trace:
            if os.path.exists( trace_file ):
                pass
            else:
                pass


    def save(self):

        pass


class Channel(base, ChannelMixIn):

    def __init__(self):
        self.raw_data = None    # raw data from ABI
        self.data = None        # smoothed data using savitzky-golay and baseline correction
        self.marker = None
        self.dye = None
        self.alleles = None

        self._assay = None


    def get_allele_class(self):
        return Allele

    def get_raw_data(self):
        """ lazy loading of raw data """

        if self.raw_data is None:
            # load from yaml file
            pass

        return self.raw_data



class Allele(base):

    """ follow the structure of msaf Allele database schema
    """

    def __init__(self, bin=-1, asize=-1, aheight=-1, size=-1,
                rtime=-1, brtime=-1, ertime=-1, wrtime=-1, srtime=-1,
                height=-1, area=-1, beta=-1, delta=0,
                    type=None, method=None, marker=None):
        self.bin = bin
        self.asize = asize          # adjusted size from reference
        self.aheight = aheight      # adjusted height
        self.size = size            # real size
        self.rtime = rtime          # retention time
        self.brtime = brtime        # begin retention time
        self.ertime = ertime        # end retention time
        self.wrtime = wrtime        # width of peak by retention time
        self.srtime = srtime        # symmetrical of peak by retention time
        self.height = height        # real height
        self.area = area            # area
        self.beta = beta            # beta of peak, area / height
        self.delta = delta          # deviation from bin point size
        self.type = type            # type of peak
        self.method = method
        self.marker = marker

    def __repr__(self):
        return '<Allele rtime: %d height: %d>' % (self.rtime, self.height)


class Marker(base):

    """ Marker information
    """

    def __init__(self, code, min_size, max_size, repeats, bins):
        self.code = code
        self.min_size = min_size
        self.max_size = max_size
        self.repeats = repeats
        self.bins = bins
        # bins is [ [pos, tag], [pos, tag], ... ]




# the filesystem-based database for sample and allele
#
# root/
#   meta.yaml
#   sample.tab
#   assay.tab
#   panels.yaml
#       
#   data/
#       sample_1/
#           assay_1/assay_1.fsa         <- the assay file,
#                   meta.yaml           <- containing metadata about assay
#                   channels/
#                       6-FAM/
#                           traces.yaml
#                           data.yaml
#                       NED/
#                           traces.yaml
#                           data.yaml
#
#           assay_2/
#                   channels.yaml

def load_sample_manifest( filename ):
    """
    loading sample manifest
    """

    pass


def load_assay_manifest( filename ):
    """
    loading assay manifest
    """

    pass


def load_channel_manifest( filename ):
    """
    loading channel manifest
    """

    pass


def load_assay( dirname ):
    """ return Assay instance """

    # get the last name

    assay_name = 
    assay_file = assay_name + '.fsa'
    assay_path = "%s/%s" % (assay_name, assay_file)
    

    # 


def load_channels( yaml_file ):
    """
    load yaml containing channels and alleles
    """
    pass


class fsdb(object):
    
    def __init__(self, rootdir):
        self.rootdir = rootdir

    def 
