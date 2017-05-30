import numpy as np
import math

from fatools.lib.utils import cverr
from fatools.lib import const
from fatools.lib.fautil.hcalign import align_hc
from fatools.lib.fautil.gmalign import align_gm, align_sh, align_de

from scipy import signal, ndimage
from peakutils import indexes

import attr

@attr.s(repr=False)
class Peak(object):
    rtime = attr.ib(default=-1)
    rfu = attr.ib(default=-1)
    area = attr.ib(default=-1)
    brtime = attr.ib(default=-1)
    ertime = attr.ib(default=-1)
    srtime = attr.ib(default=-1)
    beta = attr.ib(default=-1)
    theta = attr.ib(default=-1)
    omega = attr.ib(default=-1)

    size = attr.ib(default=-1)
    bin = attr.ib(default=-1)

    def __repr__(self):
        return "<P: %4d | %4d | %5d | %2d | %+3.2f | b %4.1f | t %4.2f | o %3d>" % (
            self.rtime, self.rfu, self.area, self.ertime - self.brtime, self.srtime,
            self.beta, self.theta, self.omega)

@attr.s
class Channel(object):
    data = attr.ib()
    marker = attr.ib()
    alleles = attr.ib(default=list)

    fsa = attr.ib(default=None)


    def scan(self, params, offset=0):

        if self.is_ladder():
            alleles = scan_peaks(self, params.ladder)
        else:
            alleles = scan_peaks(self, params.ladder, offset)

        cverr(1, "# scanning %s: %d peak(s)" % (self.marker, len(alleles)))

        return alleles


    def align(self, params):

        if not self.is_ladder():
            raise RuntimeError('ERR: cannot align non-ladder channel')

        ladders, qcfunc = self.fsa.get_ladder_parameter()

        result = align_peaks(self, params.ladder, ladders, qcfunc)





def scan_peaks(channel, params, offset=0):
    """
    """

    # check if channel is ladder channel, and adjust expected_peak_number accordingly
    expected_peak_number = params.expected_peak_number
    if channel.is_ladder():
        expected_peak_number = len(channel.fsa.panel.get_ladder()['sizes'])

    initial_peaks = find_peaks(channel.data, params, offset, expected_peak_number)

    # create alleles based on these peaks
    alleles = []
    for p in initial_peaks:
        allele = channel.Allele(
            rtime = p.rtime,
            rfu = p.rfu,
            area = p.area,
            brtime = p.brtime,
            ertime = p.ertime,
            wrtime = p.wrtime,
            srtime = p.srtime,
            beta = p.beta,
            theta = p.theta,
            omega = p.omega,
        )
        allele.type = const.peaktype.scanned
        allele.method = const.binningmethod.notavailable
        allele.marker = channel.marker
        channel.add_allele( allele )
        alleles.append( allele )

    channel.status = const.channelstatus.scanned
    return alleles


def align_peaks(channel, params, ladder, anchor_pairs=None):
    """
    returns (score, rss, dp, aligned_peak_number)
    """

    alleles = channel.get_alleles()

    # reset all peaks first
    for p in channel.get_alleles():
        p.size = -1
        p.status = const.peaktype.scanned
        print(p)

    #anchor_pairs = pairs
    from fatools.lib.fautil import pmalign
    return pmalign.align_pm( alleles, ladder )

    if anchor_pairs:
        return align_gm( alleles, ladder, anchor_pairs)

    result = align_hc( alleles, ladder )

    if result.score > 0.9:
        return result

    if result.initial_pairs:
        result = align_gm( alleles, ladder, result.initial_pairs )
        if result.score > 0.75:
            return result

    result = align_sh( alleles, ladder )
    if result.score > 0.75:
            return result

    # perform differential evolution
    return align_de( alleles, ladder )

    raise RuntimeError


    result = hclust_align( alleles, ladder )

    # add relevant info to peaks
    aligned_peaks = result[2][3]
    f = np.poly1d( result[2][2] )
    for (size, p) in aligned_peaks:
        p.dev = abs( f(p.rtime) - size)
        p.size = size
    z, rss = estimate_z( [p[1].rtime for p in aligned_peaks],
                          [p[0] for p in aligned_peaks], 3)
    print('>>> RSS:', rss)
    #import pprint; pprint.pprint( aligned_peaks )
    return result





# helper functions

def find_raw_peaks(data, params, offset, expected_peak_number=0):
    """
    params.min_dist
    params.norm_thres
    params.min_rfu
    params.max_peak_number
    """
    print("expected:", expected_peak_number)
    # cut and pad data to overcome peaks at the end of array
    obs_data = np.append(data[offset:], [0,0,0])
    if expected_peak_number:
        min_dist = params.min_dist
        indices = []
        norm_threshold = params.norm_thres
        expected_peak_number = expected_peak_number * 1.5
        while len(indices) <= expected_peak_number and norm_threshold > 1e-7:
            indices = indexes( obs_data, norm_threshold, min_dist)
            print(len(indices), norm_threshold)
            norm_threshold *= 0.5
    else:
        indices = indexes( obs_data, params.norm_thres, params.min_dist)
    cverr(3, '## raw indices: %d' % len(indices))

    if len(indices) == 0:
        return []

    # normalize indices
    if offset > 0:
        indices += offset

    # filter peaks by minimum rfu, and by maximum peak number after sorted by rfu
    peaks = [ Peak( int(i), int(data[i]) ) for i in indices
            if (data[i] > params.min_rfu and params.min_rtime < i < params.max_rtime) ]
    #peaks = sorted( peaks, key = lambda x: x.rfu )[:params.max_peak_number * 2]

    cverr(3, '## peak above min rfu: %d' % len(peaks))

    return peaks


def find_peaks(data, params, offset=0, expected_peak_number=0):

    peaks = find_raw_peaks(data, params, offset, expected_peak_number)

    # check for any peaks
    if not peaks:
        return peaks

    # measure peaks parameters
    measure_peaks(peaks, data, offset)

    #import pprint; pprint.pprint(peaks)

    # filter artefact peaks
    non_artifact_peaks = filter_for_artifact(peaks, params, expected_peak_number)

    # for ladder, special filtering is applied
    if params.expected_peak_number:
        peaks = filter_for_ladder(non_artifact_peaks, params)
    else:
        peaks = non_artifact_peaks

    return peaks


def measure_peaks(peaks, data, offset=0):

    (q50, q70) = np.percentile( data[offset:], [50, 75] )
    for p in peaks:
        p.area, p.brtime, p.ertime, p.srtime, ls, rs = calculate_area( data,
                                    p.rtime, 5e-2, q50 )
        p.wrtime = p.ertime - p.brtime
        p.beta = p.area / p.rfu
        if p.wrtime == 0:
            p.theta = 0
            p.omega = 0
        else:
            p.theta = p.rfu / p.wrtime
            p.omega = p.area / p.wrtime


def calculate_area(y, t, threshold, baseline):
    """ return (area, brtime, ertime, srtime)
        area: area
        brtime: begin rtime
        ertime: end rtime
    """

    # right area
    data = y[t:]
    r_area, ertime, r_shared = half_area(data, threshold, baseline)

    # left area
    data = y[:t+1][::-1]
    l_area, brtime, l_shared = half_area(data, threshold, baseline)


    return ( l_area + r_area - y[t], t - brtime, ertime + t, math.log2(r_area / l_area),
                l_shared, r_shared )


def half_area(y, threshold, baseline):
    """ return (area, ertime, shared_status)
    """

    winsize = 3
    threshold = threshold/2
    shared = False
    area = y[0]
    edge = float(np.sum(y[0:winsize]))/winsize
    old_edge = 2 * edge

    index = 1
    limit = len(y)

    while ( edge > area * threshold and edge < old_edge and
            index < limit and y[index] >= baseline ):
        old_edge = edge
        area += y[index]
        edge = float(np.sum(y[index:index+winsize]))/winsize
        index += 1
    if edge >= old_edge:
        shared = True
    index -= 1

    return area, index, shared


def filter_for_artifact(peaks, params, expected_peak_number = 0):
    """
    params.max_peak_number
    params.artifact_ratio
    params.artifact_dist ~ 5
    """

    # if the first peak is small, usually it is a noise
    if len(peaks) > expected_peak_number and peaks[0].theta < 5 and peaks[0].omega < 30:
        peaks = peaks[1:]

    if len(peaks) == expected_peak_number:
        return peaks

    # we need to adapt to the noise level of current channel
    if expected_peak_number > 0 and len(peaks) > expected_peak_number:
        #  gather min_theta and min_omega
        theta = sorted( [p.theta for p in peaks], reverse = True )[ expected_peak_number]
        omega = sorted( [p.omega for p in peaks], reverse = True )[ expected_peak_number]
        # test for avery low quality FSA as we need to adapt
        if theta < 4:
            theta = 0.20 * theta
        if omega < 50:
            omega = 0.25 * omega

        min_theta = min(theta, params.min_theta)
        min_omega = min(omega, 100)
        min_theta_omega = min(500, 0.25 * theta * omega)
        print('min theta & omega:', min_theta, min_omega, min_theta_omega)
    else:
        min_theta = 0
        min_omega = 0
        min_theta_omega = 0

    # filter for too sharp/thin peaks
    filtered_peaks = []
    for p in peaks:
        #filtered_peaks.append(p); continue
        print(p)
        if min_theta and min_omega and p.omega < min_omega and p.theta < min_theta:
            print('! omega & theta')
            continue
        if min_theta_omega and p.theta * p.omega < min_theta_omega:
            print('! theta_omega')
            continue
        #if p.theta < 1.5 and p.area < 100:
        #    continue
        if p.wrtime < 3:
            continue
        if p.rfu >= 25 and p.beta * p.theta < 6:
            continue
        if p.rfu < 25 and p.beta * p.theta < 3:
            continue
        #if p.omega < 50:
        #    continue
        #if p.omega < 100 and p.theta < 5:
        #    continue
        if ( params.max_beta and min_theta and
                (p.beta > params.max_beta and p.theta < min_theta) ):
            print('! max_beta')
            continue
        print('->')
        filtered_peaks.append(p)

    #import pprint; pprint.pprint(filtered_peaks)

    # filter for distance between peaks and their rfu ratio
    peaks = sorted(filtered_peaks, key = lambda x: x.rtime)
    non_artifact_peaks = []
    for idx in range(len(peaks)):
        p = peaks[idx]

        if idx > 0:
            prev_p = peaks[idx-1]
            if ( p.brtime - prev_p.ertime < params.artifact_dist
                    and p.rfu < params.artifact_ratio * prev_p.rfu ):
                # we are artifact, just skip
                continue

        if idx < len(peaks)-1:
            next_p = peaks[idx+1]
            if ( next_p.brtime - p.ertime < params.artifact_dist
                    and p.rfu < params.artifact_ratio * next_p.rfu ):
                # we are artifact, just skip
                continue

        non_artifact_peaks.append( p )

    #import pprint; pprint.pprint(non_artifact_peaks)
    #print(len(non_artifact_peaks))

    #sorted_peaks = sorted( non_artifact_peaks, key = lambda x: (x.rtime, x.beta * x.theta),
    #                    reverse=True )#[:params.max_peak_number]
    #peaks = sorted( sorted_peaks, key = lambda x: x.rfu )
    peaks = non_artifact_peaks

    cverr(3, '## non artifact peaks: %d' % len(peaks))

    return peaks


def filter_for_ladder(peaks, params):
    """
    we need to obtaine enough peaks for ladder alignment purpose, but not too much to avoid
    excessive alignment process and potentially incorrect alignment

    peaks must in rtime ascending order
    """

    epn = params.expected_peak_number   # this is the number of ladder peaks

    #
    return peaks


def baseline_als(y, lam, p, niter=10):
    pass


@attr.s
class NormalizedTrace(object):
    signal = attr.ib()
    baseline = attr.ib()

    def get_qc(self):
        """ return tuple of qcfunc
        """
        return tuple()


def normalize_baseline( raw, medwinsize=399, savgol_size=11, savgol_order=5,
                tophat_factor = 0.01 ):
    """
    params.medwin_size
    params.savgol_order
    params.savgol_size
    """

    median_line = signal.medfilt(raw, [medwinsize])
    baseline = signal.savgol_filter( median_line, medwinsize, savgol_order)
    corrected_baseline = raw - baseline
    np.maximum(corrected_baseline, 0, out=corrected_baseline)
    savgol = signal.savgol_filter(corrected_baseline, savgol_size, savgol_order)
    smooth = ndimage.white_tophat(savgol, None,
                    np.repeat([1], int(round(raw.size * tophat_factor))))

    return NormalizedTrace( signal=smooth, baseline = baseline )


@attr.s
class TraceChannel(object):
    dye_name = attr.ib()
    dye_wavelength = attr.ib()
    raw_channel = attr.ib()
    smooth_channel = attr.ib()


def b(txt):
    """ return a binary string aka bytes """
    return txt.encode('UTF-8')


def separate_channels( trace ):
    # return a list of [ 'dye name', dye_wavelength, numpy_array, numpy_smooth_baseline ]

    results = []
    for (idx, data_idx) in [ (1,1), (2,2), (3,3), (4,4), (5,105) ]:
        try:
            dye_name = trace.get_data(b('DyeN%d' % idx)).decode('UTF-8')
            dye_wavelength = trace.get_data(b('DyeW%d' % idx))
            raw_channel = np.array( trace.get_data(b('DATA%d' % data_idx)) )
            nt = normalize_baseline( raw_channel )

            results.append(
                TraceChannel(dye_name, dye_wavelength, raw_channel, nt.signal)
            )
        except KeyError:
            pass

    return results


def generate_scoring_function( strict_params, relax_params ):

    def _scoring_func( dp_result, method ):
        # alignment_result is (dp_score, dp_rss, dp_z, dp_peaks)
        dp_score = dp_result.dpscore
        dp_rss = dp_result.rss
        dp_peaks = dp_result.sized_peaks

        if method == 'strict':
            if ( dp_score >= strict_params['min_dpscore'] and
                    dp_rss <= strict_params['max_rss'] and
                    len(dp_peaks) >= strict_params['min_sizes'] ):
                return (1, None)
            return (0, None)
        elif method == 'relax':
            msg = []
            # scoring based on parts of results

            # score based on DP score compared to minimum DP score
            delta_score = relax_params['min_dpscore'] - dp_score
            if delta_score <= 0:
                dp_score_part = 1
            else:
                dp_score_part = 1e-2 ** (1e-2 * delta_score)

            # score based on RSS compared to the maximum allowed RSS
            delta_rss = dp_rss - relax_params['max_rss']
            if delta_rss <= 0:
                dp_rss_part = 1
            else:
                dp_rss_part = 1e-2 ** ( 1e-3 * delta_rss )
                msg.append( 'RSS > %d' % ( relax_params['max_rss'] ) )

            # score based on how many peaks we might miss compared to minimum number of peaks
            delta_peaks = relax_params['min_sizes'] - len(dp_peaks)
            if delta_peaks <= 0:
                dp_peaks_part = 1
            else:
                dp_peaks_part = max( 0, - delta_peaks / 0.5 * relax_params['min_sizes'] - 1)
                msg.append( 'Missing peaks = %d' % delta_peaks )

            # total overall score
            score = 0.3 * dp_score_part + 0.5 * dp_rss_part + 0.2 * dp_peaks_part
            return (score, msg)

        raise RuntimeError("Shouldn't be here!")


    return _scoring_func

## this is a new algorithm and steps to perform peak analysis
##
## fsa = import_fsa()
## ladder_channel = fsa.ladder_channel()
## alleles = scan_peaks(ladder_channel, params)
## alleles = preannotate_peaks(ladder_channel, params)
## result = align_ladder(ladder_channel, params, size_standards)
##
## for channel in fsa.non_ladder_channel():
##     scan_peaks(channel, params)
##     preannotate_peaks(channel, params)
##     call_peaks(channel, params)
##     bin_peaks(channel, params)
##     postannotate_peaks(channel, params)

## the high level methods
##
##  fsa = import_fsa()
##  fsa.align_ladder(params.ladder)
##  fsa.scan_peaks(params.nonladder, marker=None)
##  fsa.preannotate_peaks(params.nonladder, marker=None)
##  fsa.call_peaks(params.nonladder, marker=None)
##  fsa.bin_peaks(params.nonladder, marker=None)