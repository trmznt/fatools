

class peaktype(object):
    scanned = 'peak-scanned'
    broad = 'peak-broad'
    noise = 'peak-noise'
    overlap = 'peak-overlap'
    unassigned = 'peak-unassigned'
    ladder = 'peak-ladder'
    called = 'peak-called'
    stutter = 'peak-stutter'
    artifact = 'peak-artifact'
    bin = 'peak-bin'
    ignored = 'peak-ignored'


class channelstatus(object):
    unassigned = 'channel-unassigned'
    assigned = 'channel-assigned'
    unused = 'channel-unused'
    noisy = 'channel-noisy'
    empty = 'channel-empty'
    scanned = 'channel-scanned'
    aligned = 'channel-aligned'         # ladder peaks has been aligned to standard size
    called = 'channel-called'
    binned = 'channel-binned'
    ladder = 'channel-ladder'


class assaystatus(object):
    uploaded = 'assay-uploaded'
    assigned = 'assay-assigned'
    scanned = 'assay-scanned'
    preannotated = 'assay-preannotated'
    aligned = 'assay-aligned'
    called = 'assay-called'
    binned = 'assay-binned'
    annotated = 'assay-annotated'


class allelemethod(object):
    uncalled = 'allele-uncalled'


class binningmethod(object):
    unavailable = 'binning-unavailable'
    auto = 'binning-auto'
    semiauto = 'binning-semiauto'


dyes = [ '6-FAM', 'NED', 'VIC', 'PET', 'LIZ' ]

ladders = { 'LIZ600': { 'dye': 'LIZ',
                        'sizes': [ 20.0, 40.0, 60.0, 80.0, 100.0, 114.0, 120.0, 140.0, 160.0,
                        180.0, 200.0, 214.0, 220.0, 240.0, 250.0, 260.0, 280.0, 300.0,
                        314.0, 320.0, 340.0, 360.0, 380.0, 400.0, 414.0, 420.0, 440.0,
                        460.0, 480.0, 500.0, 514.0, 520.0, 540.0, 560.0, 580.0, 600.0 ],
                        'strict': {
                            'max_rss': 40.0,
                            'min_dpscore': 34.0,
                            'min_sizes': 36
                        },
                        'relax': {
                            'max_rss': 56.25,
                            'min_dpscore': 33.0,
                            'min_sizes': 36
                        }
                    },
            'LIZ500': { 'dye': 'LIZ',
                        'sizes': [ 35, 50, 75, 100, 139, 150, 160, 200, 250, 300, 340, 350,
                        400, 450, 490, 500 ],
                        'strict': {
                            'max_rss': 17.5,
                            'min_dpscore': 14.0,
                            'min_sizes': 16
                        },
                        'relax': {
                            'max_rss': 25.0,
                            'min_dpscore': 13.0,
                            'min_size': 16
                        }                   
    }

