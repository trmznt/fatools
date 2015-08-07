

class peaktype(object):
    scanned = 'scanned'
    broad = 'broad'
    noise = 'noise'
    overlap = 'overlap'
    unassigned = 'unassigned'
    ladder = 'ladder'
    called = 'called'
    stutter = 'stutter'
    artifact = 'artifact'
    bin = 'bin'
    ignored = 'ignored'


class channelstatus(object):
    unassigned = 'unassigned'
    assigned = 'assigned'
    unused = 'unused'
    noisy = 'noisy'
    empty = 'empty'
    reseted = 'reseted'
    scanned = 'scanned'
    preannotated = 'preannotated'
    aligned = 'aligned'         # ladder peaks has been aligned to standard size
    called = 'called'
    binned = 'binned'
    annotated = 'annotated'
    ladder = 'ladder'


class assaystatus(object):
    uploaded = 'uploaded'
    assigned = 'assigned'
    scanned = 'scanned'
    preannotated = 'preannotated'
    aligned = 'aligned'
    called = 'called'
    binned = 'binned'
    annotated = 'annotated'


class alignmethod(object):
    notapplicable = 'notapplicable'
    fast_hq = 'fast|highqual'
    fast_mq = 'fast|medqual'
    fast_hqr = 'fast|highqual-relax'
    fast_mqr = 'fast|medqual-relax'
    greedy_filtered = 'greedy|filtered'
    greedy_shifted = 'greedy|shifted'
    greedy_scored = 'greedy|scored'


class scanningmethod(object):
    notapplicable = 'notapplicable'
    cwt = 'cwt'


class allelemethod(object):
    uncalled = 'uncalled'
    leastsquare = 'leastsquare'
    cubicspline = 'cubicspline'
    localsouthern = 'localsouthern'


class binningmethod(object):
    unavailable = 'unavailable'
    auto = 'auto'
    semiauto = 'semiauto'


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
                            'min_sizes': 16
                        }                   
                    },
}

