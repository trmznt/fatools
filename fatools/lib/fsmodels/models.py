# basic classes that mimic msdb classes
# this is filesystem-based database


from fatools.lib.fautil.mixin import (AssayMixIn, ChannelMixIn,
            MarkerMixIn, PanelMixIn, AlleleSetMixIn)


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
        pass        


class Sample(base):

    def __init__(self):
        self.assays = None


class Assay(base, AssayMixIn):

    def __init__(self):
        self.channels = []
        self.rss = -1
        self.dp = -1
        self.score = -1
        self.size_standard = None

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


    def new_channel(self, raw_data, data, dye, wavelen, status,
            median, mean, max_height, min_height, std_dev):
        channel = Channel()
        channel.raw_data = raw_data
        channel.data = data
        channel.dye = dye
        channel.wavelen = wavelen
        channel.status = status
        channel.median = median
        channel.mean = mean
        channel.max_height = max_height
        channel.min_height = min_height
        channel.std_dev = std_dev
        channel._assay = self

        # need to assign initial marker
        channel.marker = undefined_marker

        self.channels.append( channel )
        return channel



class Channel(base, ChannelMixIn, AlleleSetMixIn):

    def __init__(self):
        self.raw_data = None    # raw data from ABI
        self.data = None        # smoothed data using savitzky-golay and baseline correction
        self.marker = None
        self.dye = None

        self._assay = None


    def get_allele_class(self):
        return Allele

    def get_raw_data(self):
        """ lazy loading of raw data """

        if self.raw_data is None:
            # load from yaml file
            pass

        return self.raw_data


    def new_alleleset(self):
        # note: we don't really have AlleleSet, just return ourselves
        return self

    def get_latest_alleleset(self):
        return self


    def new_allele(self, rtime, height, area, brtime, ertime, wrtime, srtime, beta, theta,
                    type, method):
        allele = Allele( rtime = rtime, height = height, area = area,
                    brtime = brtime, ertime = ertime, wrtime = wrtime, srtime = srtime,
                    beta = beta, theta = theta, type = type, method = method )
        allele.alleleset = self

        return allele



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


class Marker(base, MarkerMixIn):

    """ Marker information
    """

    def __init__(self, code, min_size, max_size, repeats, bins):
        self.code = code
        self.species = 'x'
        self.min_size = min_size
        self.max_size = max_size
        self.repeats = repeats
        self.bins = bins
        # bins is [ [pos, tag], [pos, tag], ... ]


undefined_marker = Marker('undefined', 10, 600, 0, [])

class Panel(base, PanelMixIn):

    """ Panel information
    """

    def __init__(self, code, data ):
        self.code = code
        self.data = data
        self._dyes = {}

    def get_marker(self, code):
        print('creating marker: %s' % code)
        m = Marker(code, 10, 600, 0, None)
        return m

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

    assay_name = None
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

