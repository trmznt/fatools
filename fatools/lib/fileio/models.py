# fileio/models.py
#
# models here are used to ensure the integrity and consistency of each inherited model
#


from fatools.lib.utils import cout, cerr
from fatools.lib.fautil.mixin2 import MarkerMixIn, PanelMixIn, ChannelMixIn, FSAMixIn, AlleleMixIn

import os, pickle


class Marker(MarkerMixIn):

    container = {}

    @classmethod
    def upload(cls, d):
        for (k,v) in d.items():
            cls.container[k] = cls.from_dict(v)


    @classmethod
    def get_marker(cls, marker_code, species='x'):
        if '/' not in marker_code:
            marker_code = species + '/' + marker_code
        return cls.container[marker_code]


class Panel(PanelMixIn):

    container = {}

    Marker = Marker


    @classmethod
    def upload(cls, d):
        for (k,v) in d.items():
            cls.container[k] = cls.from_dict(v)


    @classmethod
    def get_panel(cls, panel_code):
        return cls.container[panel_code]


class Allele(AlleleMixIn):

    def __init__(self, rtime, rfu, area, brtime, ertime, wrtime, srtime,
                    beta, theta, omega):
        self.rtime = rtime
        self.rfu = rfu
        self.area = area
        self.brtime = brtime
        self.ertime = ertime
        self.wrtime = wrtime
        self.srtime = srtime
        self.beta = beta
        self.theta = theta
        self.omega = omega

        self.size = -1
        self.bin = -1
        self.dev = -1


class Channel(ChannelMixIn):

    Allele = Allele

    def __init__(self, data, dye, wavelen, status, fsa):
        self.data = data
        self.dye = dye
        self.wavelen = wavelen
        self.status = status
        self.fsa = fsa

        self.alleles = []

        self.assign()


    def add_allele(self, allele):
        self.alleles.append(allele)
        return allele



class FSA(FSAMixIn):

    __slots__ = [ '_fhdl', '_trace' ]

    Channel = Channel

    def __init__(self):
        self.channels = []
        self.excluded_markers = []

    def get_data_stream(self):
        return self._fhdl

    def add_channel(self, channel):
        self.channels.append(channel)
        return channel

    @classmethod
    def from_file(cls, fsa_filename, panel, excluded_markers=None, cache=True):
        fsa = cls()
        fsa.filename = os.path.basename(fsa_filename)
        fsa._fhdl = open(fsa_filename, 'rb')
        fsa.set_panel(panel, excluded_markers)

        # with fileio, we need to prepare channels everytime or seek from cache
        cache_file = '.fatools_caches/channels/%s' % fsa.filename
        if cache and os.path.exists(cache_file):
            if os.stat(fsa_filename).st_mtime < os.stat(cache_file).st_mtime:
                cerr('I: uploading channel cache for %s' % fsa_filename)
                fsa.channels = pickle.load( open(cache_file, 'rb') )
                for c in fsa.channels:
                    c.fsa = fsa
        else:
            fsa.create_channels()
            if cache and os.path.exists('.fatools_caches/channels'):
                for c in fsa.channels: c.fsa = None
                pickle.dump(fsa.channels, open(cache_file, 'wb'))
                for c in fsa.channels: c.fsa = fsa
        return fsa



