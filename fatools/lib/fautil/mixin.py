
from fatools.lib.fautil import traceio, traceutils
from fatools.lib.utils import cout, cerr
from fatools.lib.const import (peaktype, channelstatus, assaystatus, dyes, ladders,
                                    allelemethod, alignmethod, binningmethod, scanningmethod)
from fatools.lib.fautil import algo
from sortedcontainers import SortedListWithKey

import io, numpy as np
from copy import copy
from functools import lru_cache
import pprint, sys, time


class PanelMixIn(object):
    """ contains Panel methods """

    def update(self, obj):
        raise NotImplementedError('PROG/ERR - child class must provide this method')

    def _update(self, obj):

        if type(obj) == dict:
            if self.id is not None:
                raise RuntimeError('ERR: can only update from dictionary on new instance')
            self.code = obj['code']
            self.data = obj['data']
        else:
            if self.code != obj.code:
                raise RuntimeError('ERR: attempting to update Panel with different code!')
            if  obj.data is not None:
                self.data = obj.data


    def get_marker_codes(self):
        """ return a list of marker codes """
        if self.data:
            return list( self.data['markers'].keys() )
        return []


    def get_markers(self):
        """ return a list of markers used in this panel """
        return [ self.get_marker(m) for m in self.get_marker_codes() ]


    def get_marker(self, marker_code):
        """ return marker instance with marker_code """
        raise NotImplementedError('PROG/ERR - child class must provide this method')


    def get_ladder_code(self):
        if self.data:
            return self.data['ladder']
        return ''


    def get_marker_by_dye(self, dye):
        if True:
            for m in self.get_markers():
                self._dyes[ self.data['markers'][m.label]['dye'] ] = m
        return self._dyes[dye]



class MarkerMixIn(object):
    """ contains Marker method """

    def update(self, obj):
        raise NotImplementedError('PROG/ERR - child class must provide this method')

    def _update(self, obj):
        """ supplanted by child's update() method, if necessary """

        # check if obj is a dict type or not

        if type(obj) == dict:
            if self.id is not None:
                raise RuntimeError('ERR: can only update from dictionary on new instance')
            if 'code' in obj:
                self.code = obj['code']
            if 'species' in obj:
                self.species = obj['species']
            if 'min_size' in obj:
                self.min_size = obj['min_size']
            if 'max_size' in obj:
                self.max_size = obj['max_size']
            if 'repeats' in obj:
                self.repeats = obj['repeats']
            #if 'bins' in obj:
            #    self.bins = obj['bins']
            if 'z_params' in obj:
                self.z_params = obj['z_params']

            # field 'related_to' must be handled by implementation object

        else:
            if self.code != obj.code:
                raise RuntimeError('ERR: attempting to update Marker with different code')
            if obj.min_size is not None:
                self.min_size = obj.min_size
            if obj.max_size is not None:
                self.max_size = obj.max_size
            if obj.repeats is not None:
                self.repeats = obj.repeats
            #if obj.bins is not None:
            #    self.bins = obj.bins
            if obj.related_to is not None:
                self.related_to = obj.related_to
            if obj.z_params is not None:
                self.z_params = obj.z_params

    @property
    def label(self):
        return self.species + '/' + self.code

    #@property
    #def sortedbins(self):
    #    return SortedListWithKey(self.bins, key = lambda b: b[1])

    @lru_cache(maxsize=32)
    def get_sortedbins(self, batch):
        # get Bin from this batch
        bin = self.get_bin(batch)
        return bin.sortedbins

    def initbins(self, start_range, end_range, batch):
        # check whether we already have bin data, otherwise create new one
        bin = self.get_bin(batch, recursive=False)
        if bin is None:
            bin = self.new_bin( batch = batch )
        bin.initbins( start_range, end_range, self.repeats )


    def new_bin(self, batch):
        raise NotImplementedError('PROG/ERR - need to implement this on child class')



class BinMixIn(object):

    def initbins(self, start_range, end_range, repeats):
        self.bins = []
        for size in range(start_range, end_range, repeats):
            self.bins.append( [size, float(size), float(size-1), float(size+1)] )
            # the real bins are defined by [ bin_value, mean, 25percentile, 75percentile ]

    def adjustbins(self, updated_bins):
        bins = copy(self.bins)
        for idx in range(len(bins)):
            tag = bins[idx][1]
            if tag in updated_bins:
                bins[idx] = updated_bins[tag]
                cerr('updating bin: %d' % tag)

        # force db to update
        self.bins = bins

    @property
    def sortedbins(self):
        return SortedListWithKey(self.bins, key = lambda b: b[1])



class BatchMixIn(object):
    """ contains Batch methods """

    def get_panel(self, panel_code):
        """ shortcut to get single Panel instance, otherwise throws exception """
        raise NotImplementedError('PROG/ERR - child class must override this method')


    def get_marker(self, marker_code, species_code = None):
        """ shortcut to get single Marker instance, otherwise throws exception """
        raise NotImplementedError('PROG/ERR - child class must override this method')

    @property
    def sample_ids(self):
        """ this has to override by sub-class to provide more efficient implementation
            based on the underlying database architecture
        """
        return list([ x.id for s in self.samples])


    def update(self, obj):
        raise NotImplementedError('PROG/ERR - child class must provide this method')

    def _update(self, obj):

        if type(obj) == dict:
            # updating from dictionary (eg. YAML )
            if 'code' in obj:
                self.code = obj['code']
            if 'description' in obj:
                self.description = obj['description']
            if 'remark' in obj:
                self.remark = obj['remark']

        else:
            raise NotImplementedError('PROG/ERR - not implemented yet')


    def remove_assays(self):
        for sample in self.samples:
            sample.remove_assays()



class SampleMixIn(object):
    """ implement general Sample methods """

    def add_fsa_assay(self, trace, filename, panel_code, options = None, species = None,
                    dbhandler = None, dry_run=False):

        assert dbhandler, 'Please provide dbhandler function'

        # check panel

        batch = self.batch

        panels = [ dbhandler.get_panel(c.strip()) for c in panel_code.split(',') ]

        # parsing options

        excluded_markers = []
        if options:
            for key, val in options.items():
                if key.startswith('exclude'):
                    for marker_code in val.split(','):
                        if '/' not in marker_code:
                            if not species:
                                raise RuntimeError('ERR - need to specifiy species')
                            marker_code = species + '/' + marker_code.strip()
                        else:
                            marker_code = marker_code.strip()
                        excluded_markers.append(marker_code.upper())


        # Processing options
        if excluded_markers:
            # check excluded markers
            panel_markers = []
            for panel in panels:
                panel_markers.extend(
                    ('%s/%s' % (species, m) if '/' not in m else m).upper()
                    for m in panel.get_marker_codes()
                )
            print(str(excluded_markers), str(panel_markers))
            unknown_markers = set(excluded_markers) - set(panel_markers)
            if unknown_markers:
                raise RuntimeError('ERR - assay %s does not have exluded marker(s): %s'
                            % (filename, ','.join( unknown_markers )))

        if dry_run:
            # we need to check whether trace is valid
            trace_instance = traceio.read_abif_stream( io.BytesIO( trace ) )
            # things to be checked:
            # - dyes in the panel exist in DB
            # - duplicates of filename, panel_id, sample_id in the current DB
            # - duplicates of filename, panel_id, sample_id in the assay list -> check in caller
            ## XXX: fixme
            return None

        # creating assay
        for panel in panels:
            assay = self.new_fsa_assay(raw_data = trace, filename = filename,
                        status = assaystatus.uploaded, panel = panel)
            assay.status = assaystatus.uploaded
            assay.method = alignmethod.notapplicable
            assay.runtime = assay.get_trace().get_run_start_time()
            assay.create_channels()
            assay.assign_channels( excluded_markers )

        return assay


    def new_fsa_assay(self, trace, filename, panel):

        raise NotImplementedError('PROG/ERR - child class must override this method!')


    def _update(self, obj):

        if type(obj) == dict:
            if 'code' in obj:
                self.code = obj['code']
            if 'type' in obj:
                self.type = obj['type']
            if 'altcode' in obj:
                self.altcode = obj['altcode']
            if 'category' in obj:
                self.category = obj['category']
            if 'remark' in obj:
                self.remark = obj['remark']

        else:

            raise NotImplementedError('PROG/ERR - not implemented yet')


    def remove_fsas(self):
        raise NotImplementedError



class ChannelMixIn(object):
    """ contains Channel methods """


    def reset(self):
        """ reset this channel, ie. set all peaks into peak-scanned """
        raise NotImplementedError()


    def clear(self):
        """ clear this channel, ie. remove all peaks by removing all alleleset """
        raise NotImplementedError()


    def preprocess(self, params):

        raise NotImplementedError()


    def scan(self, params, peakdb = None):
        """ scan using params """

        #print('SCANNING: %s' % self.dye)

        alleleset = self.new_alleleset()
        alleleset.scanning_method = scanningmethod.notapplicable
        alleleset.calling_method = allelemethod.uncalled
        alleleset.binning_method = binningmethod.notavailable

        # first, check whether we are ladder or not
        if self.marker.code == 'ladder':

            ladder_code = self.assay.size_standard
            sizes = ladders[ladder_code]['sizes']
            params.ladder.max_peak_number = len(sizes) * 2
            params.ladder.expected_peak_number = len(sizes)

            alleles = algo.scan_peaks(self, params.ladder, peakdb)
            cerr('ladder: %d; ' % len(alleles), nl=False)
            alleleset.scanning_method = params.ladder.method


        else:

            alleles = algo.scan_peaks(self, params.nonladder, peakdb)
            cerr('%s: %d; ' % (self.marker.label, len(alleles)), nl=False)
            alleleset.scanning_method = params.nonladder.method


    def preannotate(self, params):
        """ preannotate must be conducted within assay, eg. need all channels """
        raise NotImplementedError()


    def alignladder(self, params):

        if self.marker.code != 'ladder':
            # sanity checking
            raise RuntimeError("ERR - can't align ladder on non-ladder channel")

        ladder_code = self.assay.size_standard
        ladder = ladders[ladder_code]
        ladder_sizes = ladder['sizes']
        ladder_qc_func = algo.generate_scoring_function( ladder['strict'], ladder['relax'] )

        # reset all calculated values
        for p in self.alleles:
            p.size = -1
            p.bin = -1

        start_time = time.process_time()
        (qcscore, remarks, results, method) = algo.size_peaks(self, params, ladder_sizes,
                                                ladder_qc_func)
        stop_time = time.process_time()
        (dpscore, rss, z, aligned_peaks) = results
        #qcscore, remarks = algo.score_ladder(rss, len(aligned_peak), len(ladder_sizes))

        f = np.poly1d(z)
        for (std_size, peak) in aligned_peaks:
            peak.size = std_size
            peak.type = peaktype.ladder
            peak.deviation = (f(peak.rtime) - std_size)**2


        # set channel & assay instance
        self.status = channelstatus.aligned
        assay = self.assay
        assay.dp = dpscore
        assay.score = qcscore
        assay.rss = rss
        assay.z = z
        assay.ladder_peaks = len(aligned_peaks)
        assay.process_time += int((stop_time - start_time) * 1000)  # in miliseconds
        assay.method = method
        if remarks:
            if assay.report:
                assay.report = assay.report + '//' + '|'.join(remarks)
            else:
                assay.report = '|'.join(remarks)

        return (dpscore, rss, len(aligned_peaks), len(ladder_sizes), qcscore, remarks, method)


    def call(self, params, func, min_rtime, max_rtime):

        # self sanity check
        if self.marker.code == 'ladder':
            return

        algo.call_peaks( self, params, func, min_rtime, max_rtime )


    def size(self, params):

        raise NotImplementedError()


    def bin(self, params):

        # self sanity check
        if self.marker.code == 'ladder':
            return

        if self.marker.code == 'undefined':
            return

        if self.marker.code == 'combined':
            # this is specific for combined marker
            # create allelesets as many as markers
            for marker in self.markers:
                alleleset = self.allelesets[0].clone()
                alleleset.marker = marker
                algo.bin_peaks( alleleset, params, marker )

        else:
            algo.bin_peaks( self.allelesets[0], params, self.marker )


    def postannotate(self, params):

        # self sanity check
        if self.marker.code == 'ladder':
            return

        if self.marker.code == 'undefined':
            return

        algo.postannotate_peaks( self.allelesets[0], params )
        #raise NotImplementedError()


    def tag(self):

        return '%s|%s|%s|%s|%s' % (self.assay.sample.batch.code,
                        self.assay.sample.code, self.assay.filename,
                        self.assay.runtime, self.dye.upper())


    def get_latest_alleleset(self):

        raise NotImplementedError('PROG/ERR - child class must override this method')


    @property
    def alleles(self):
        return self.get_latest_alleleset().alleles


    def new_allele(self, rtime, height, area, brtime, ertime, wrtime, srtime, beta, theta):
        alleleset = self.get_latest_alleleset()
        return alleleset.new_allele( rtime = rtime, height = height, area = area,
                    brtime = brtime, ertime = ertime, wrtime = wrtime, srtime = srtime,
                    beta = beta, theta = theta,
                    type = peaktype.scanned, method = binningmethod.notavailable )


    def showladderpca(self):

        import mdp
        from matplotlib import pylab as plt
        import pprint
        import math

        cerr('calculating PCA & plotting')
        peak_sizes = sorted(list([ x.rtime for x in self.alleles ]))
        #peak_sizes = sorted( peak_sizes )[:-5]
        #pprint.pprint(peak_sizes)
        #comps = algo.simple_pca( peak_sizes )
        #algo.plot_pca(comps, peak_sizes)

        from fatools.lib import const
        std_sizes = const.ladders['LIZ600']['sizes']

        x = std_sizes
        y = [ x * 0.1 for x in peak_sizes ]

        D = np.zeros( (len(y), len(x)) )
        for i in range(len(y)):
            for j in range(len(x)):
                D[i,j] = math.exp( ((x[j] - y[i]) * 0.001) ** 2 )

        pprint.pprint(D)
        im = plt.imshow(D, interpolation='nearest', cmap='Reds')
        plt.gca().invert_yaxis()
        plt.xlabel("STD")
        plt.ylabel("PEAK")
        plt.grid()
        plt.colorbar()
        plt.show()


class AssayMixIn(object):
    """ contains Assay processing method """


    def reset(self):
        """ reset all channels """
        for c in self.channels:
            c.reset()


    def clear(self):
        """ clear all channels """
        for c in self.channels:
            c.clear()
        self.status = assaystatus.assigned


    def preprocess(self, params):
        """ preprocess each channel, eg. smooth & normalize """
        for c in self.channels:
            c.preprocess(params)


    def scan(self, params, peakdb = None):
        """ scan for peaks """
        for c in self.channels:
            c.scan(params, peakdb)
        self.status = assaystatus.scanned
        cerr('')


    def preannotate(self, params):
        """ annotate peaks for broad, rtime-based stutter & overlapping peaks """
        channels =  list(self.channels)
        start_time = time.process_time()
        algo.preannotate_channels(channels, params.nonladder)
        stop_time = time.process_time()
        self.status = assaystatus.preannotated
        self.process_time = int((stop_time - start_time) * 1000) # in miliseconds


    def alignladder(self, excluded_peaks, force_mode=False):
        """ align ladder
            return ( score, rss, no of peak, no of ladder)
        """

        for c in self.channels:
            if c.marker.code != 'ladder': continue
            result = c.alignladder(excluded_peaks)
            break
        self.status = assaystatus.aligned
        return result


    def size(self, params):
        """ match ladder size with channel peaks """
        for c in self.channels:
            c.size(params)
        self.status = assaystatus.aligned


    def call(self, params, method=allelemethod.localsouthern):
        """ determine size of each peaks using ladder channel"""

        ladders = [ p for p in self.ladder.alleles if p.size > 0 ]

        # check the method
        if method == allelemethod.leastsquare:
            func = algo.least_square( ladders, self.z)
        elif method == allelemethod.cubicspline:
            func = algo.cubic_spline( ladders )
        elif method == allelemethod.localsouthern:
            func = algo.local_southern( ladders )
        else:
            raise RuntimeError
        min_rtime = ladders[1].rtime
        max_rtime = ladders[-2].rtime

        for c in self.channels:
            if c == self.ladder: continue
            c.call(params, func, min_rtime, max_rtime)
        self.status = assaystatus.called


    def postannotate(self, params, markers = None):
        """ annotate peaks for size-based stutter, etc """
        for c in self.channels:
            if c == self.ladder: continue
            if markers and c.marker not in markers:
                continue
            c.postannotate(params.nonladder)
        self.status = assaystatus.annotated


    def bin(self, params, markers = None):
        """ bin peaks to marker sizes """
        for c in self.channels:
            if c == self.ladder: continue
            if markers and c.marker not in markers:
                continue
            c.bin(params)
        self.status = assaystatus.binned


    def create_channels(self):
        """ create new channel based on current trace """
        #t = traceio.read_abif_stream( io.BytesIO( self.raw_data ) )
        t = self.get_trace()

        channels = traceutils.separate_channels( t )
        for tc in channels:
            # check n (dye name)
            if tc.dye_name not in dyes:
                raise RuntimeError('ERR - dye %s is unknown!' % n)
            c = self.new_channel(   raw_data = tc.raw_channel, data = tc.smooth_channel,
                                    dye = tc.dye_name, wavelen = tc.dye_wavelength,
                                    status = channelstatus.reseted,
                                    median=int(tc.median), mean=float(tc.mean),
                                    max_height=int(tc.max_height), min_height=int(tc.min_height),
                                    std_dev = float(tc.sd) )


    def assign_channels(self, excluded_markers=None):
        """ assign channel & ladder based on panel """

        # check panel
        panel = self.panel
        ladder_code = panel.get_ladder_code()
        ladder_dye = ladders[ladder_code]['dye']

        # check excluded_markers

        has_ladder = False
        marker_count = 0
        excluded = []

        cerr('Assay: %s' % self.filename)
        cerr('Dyes: ', nl=False)
        for channel in self.channels:
            if channel.dye.upper() == ladder_dye.upper():
                channel.marker = panel.get_marker('ladder')
                has_ladder = True
                channel.status = channelstatus.ladder
                self.size_standard = ladder_code
                self.ladder = channel
                continue

            try:
                marker = panel.get_marker_by_dye( channel.dye.upper() )
            except KeyError:
                channel.status = channelstatus.unused
                cerr('%s => Unused; ' % channel.dye, nl=False)
                continue

            if excluded_markers != None and marker.label.upper() in excluded_markers:
                channel.status = channelstatus.unassigned
                excluded.append( marker.label.upper() )
                cerr('%s => Unassigned; ' % channel.dye, nl=False)
                continue

            channel.marker = marker
            channel.status = channelstatus.assigned
            marker_count += 1
            cerr('%s => %s; ' % (channel.dye, marker.code), nl=False)

        cerr('')
        if not has_ladder:
            raise RuntimeError('ERR - sample %s assay %s does not have ladder!' %
                            (self.sample.code, self.filename))

        if excluded:
            self.exclude = ','.join( excluded )

        if self.status in [assaystatus.unassigned, assaystatus.uploaded]:
            self.status = assaystatus.assigned


    def showladderpca(self):
        for c in self.channels:
            if c.marker.code == 'ladder':
                c.showladderpca()


    def get_trace(self):
        if not hasattr(self, '_trace'):
            self._trace = traceio.read_abif_stream( io.BytesIO( self.raw_data ) )
        return self._trace



class AlleleSetMixIn(object):

    def new_allele(self, rtime, height, area, brtime, ertime, wrtime, srtime, beta):
        raise NotImplementedError('PROG/ERR - child class must override this method!')

    @property
    def batch(self):
        return self.sample.batch



class AlleleMixIn(object):

    def _update(self, obj):

        if type(obj) == dict:

            self.type_id = obj['type_id']


    def __repr__(self):
        return '<Allele [%3d] %5d %6d>' % (self.size, self.rtime, self.height)

    def __str__(self):
        return '<A [%3d %6.2f] %5d %6d | %5.1f  %+3.1f  %6.1f  %2.1f  %5d  %5.2f| %s>' % (
                self.bin, self.size, self.rtime, self.height, self.beta, self.srtime,
                self.theta, self.qscore, round(self.theta * self.beta), self.deviation,
                self.type)


class NoteMixIn(object):

    pass

class BatchNoteMixIn(object):

    pass


class SampleNoteMixIn(object):

    pass


class MarkerNoteMixIn(object):

    pass

class PanelNoteMixIn(object):

    pass

class AssayNoteMixIn(object):

    pass

class ChannelNoteMixIn(object):

    pass

class AlleleSetNoteMixIn(object):

    pass
