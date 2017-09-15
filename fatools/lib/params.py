
import numpy as np


class ScanningParameter(object):

    def __init__(self):

        self.method = 'pd'
        # 'mlpy' is fast, 'cwt' is accurate, 'relmax' is so-so, 'pd' is peak detect

        self.min_height = 1
        self.min_size = 100
        self.max_size = 600
        #self.stutter_threshold = 1.25
        self.overlap_threshold = 1.10
        self.bin_relative_ratio = 0.5
        self.widths = np.arange( 5, 15) #3, 10)
        #self.widths = [5, 7, 10, 15 ]
        self.min_snr = 3.0  # used to be 2.5
        self.min_relative_ratio = 0
        self.max_relative_ratio = 0
        self.min_height_ratio = 0
        self.min_rtime = 750
        self.max_rtime = 12500
        self.max_peak_number = 30
        self.max_beta = 15
        self.max_gradient_threshold = 0
        self.overlap_height_threshold = 0.75
        self.stutter_rtime_threshold = 10
        self.stutter_height_threshold = 0.5
        #self.stutter_size_threshold = 1.25
        self.height = -1    # if this is > 0, then the peaks will be filtered based on this peak
        self.ignoredpeaks = None
        self.width_ratio = 2000
        self.expected_peak_number = 0
        self.stutter_ratio = 0.95
        self.stutter_range = 3.5
        self.artifact_ratio = 0.8

        self.norm_thres = 0.05
        self.min_dist = 10
        self.min_rfu = 20
        self.artifact_dist = 20

        self.min_theta = 0
        self.max_beta = 0


class LadderScanningParameter(ScanningParameter):

    def __init__(self):
        super().__init__()
        self.min_rtime = 1
        self.min_relative_ratio = 0
        self.max_relative_ratio = 0
        self.min_height_ratio = 0
        self.min_peak_number = 0.85
        self.ladder_height_win = range( 1, 500)
        self.max_gradient_threhold = 20
        self.widths = np.arange( 2, 8) # (2, 8) #(3, 10 )
        #self.min_snr = 2.25    # for RELMAX
        #self.min_snr = 2.25    # for CWT
        self.min_snr = 1.25
        self.max_peak_number = 0
        self.init_separation_time = -1
        self.min_height = 6 # 50 is preferable
        self.width_ratio = 5000
        self.expected_peak_number = 36

        self.norm_thres = 0.1
        self.min_rfu = 3
        self.min_dist = 12

        self.min_theta = 4
        self.max_beta = 18

        self.artifact_dist = 15
        self.artifact_ratio = 0.5


class Params(object):

    ladder = LadderScanningParameter()
    nonladder = ScanningParameter()


default_panels = {
    'GS600LIZ': {
        'code': 'GS600LIZ',
        'data': {
            'ladder': 'LIZ600',
            'markers': {
                'x/VIC': dict(dye='VIC', filter='G'),
                'x/PET': dict(dye='PET', filter='R'),
                'x/NED': dict(dye='NED', filter='Y'),
                'x/6-FAM': dict(dye='6-FAM', filter='B'),
            }
        }
    },

    'GS500LIZ': {
        'code': 'GS500LIZ',
        'data': {
            'ladder': 'LIZ500',
            'markers': {
                'x/VIC': dict(dye='VIC', filter='G'),
                'x/PET': dict(dye='PET', filter='R'),
                'x/NED': dict(dye='NED', filter='Y'),
                'x/6-FAM': dict(dye='6-FAM', filter='B'),
            }
        }
    },

}

default_markers = {

    'x/ladder': dict( code='ladder', species='x', min_size=-1, max_size=-1 ),

    'x/VIC': dict( code='VIC', species='x', min_size=80, max_size=500 ),

    'x/PET': dict( code='PET', species='x', min_size=80, max_size=500 ),

    'x/NED': dict( code='NED', species='x', min_size=80, max_size=500 ),

    'x/6-FAM': dict( code='6-FAM', species='x', min_size=80, max_size=500 ),
}
