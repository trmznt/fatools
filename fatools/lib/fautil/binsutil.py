import pandas, attr, yaml
import numpy as np
from fatools.lib.fautil.mixin import BinMixIn
from fatools.lib.utils import cout, cerr
from collections import defaultdict

from IPython import embed

def do_binsutil(args):

    if args.optimize:
        do_optimize(args)

    elif args.summarize:
        do_summarize(args)

    else:
        cerr('ERR: operation undefined!')


def do_summarize(args):

    d = pandas.read_table(args.infile)
    peaks = d[ d['MARKER'] == args.marker ]
    stats = bin_stats(peaks)
    for s in sorted(stats.keys()):
        cout(stats[s].repr())


def do_optimize(args):

    d = pandas.read_table(args.infile)

    p = d[ d['MARKER'] == args.marker ].copy()

    counts = np.round( p['SIZE'] )

    c = defaultdict(int)
    for i in counts:
        c[int(i)] += 1

    bins = sorted(c.items(), key = lambda x: x[1])
    print(bins)

    tbin = Bin()
    tbin.initbins(args.anchor, args.repeats, args.min, args.max, args.shift)

    for i in range(0, 30):
        print('<<ITER:  %d >>' % i)
        call_peaks(tbin, p)
        stat = bin_stats(p)
        for s in sorted(stat.keys()):
            print(stat[s].repr())
        if i % 10 == 0:
            tbin.adjust_bins(stat, reset=True, repeats=args.repeats)
        else:
            tbin.adjust_bins(stat)

    if args.outfile:
        with open(args.outfile, 'w') as f:
            yaml.dump({ args.marker: { 'label': args.marker, 'bins': tbin.bins}}, f)


class Bin(BinMixIn):

    def __init__(self):
        self.bins = None

    def initbins(self, anchor, repeats, min_range, max_range, shift=0):

        mod = anchor % repeats
        print('mod', mod)
        print('min_range', min_range)
        min_range = (min_range // repeats -1) * repeats + mod
        print('min_range', min_range)

        #super().initbins(min_range, max_range, repeats)
        self.bins = []
        for i in range(min_range, max_range, repeats):
            self.bins.append([i, float(i) + shift, i - 0.5 + shift, i + 0.5 + shift])
        print(self.bins)

    def adjust_bins(self, containers, reset=False, repeats=-1):
        adjust_bins(self.bins, containers, reset, repeats)


@attr.s
class BinContainer(object):

    size = attr.ib(default=-1)
    values = attr.ib(default=attr.Factory(list))

    def repr(self):
        return "<Bin: %d / %5.4f / %5.4f / %5.4f - %5.4f d: %5.4f f: %d>" % (
            self.size, np.mean(self.values), np.median(self.values), min(self.values), max(self.values), max(self.values) - min(self.values), len(self.values)
        )

    def percentile(self, q):
        return np.percentile(self.values, q)

    def d(self):
        return max(self.values) - min(self.values)

    def f(self):
        return len(self.values)


def call_peaks(bins, peaks):

    sortedbins = bins.sortedbins

    #embed()

    for i in peaks.index:
        size = peaks.loc[i, 'SIZE']
        idx = sortedbins.bisect_key_right( size )

        if idx == 0:
            curr_bin = sortedbins[0]
        elif idx == len(sortedbins):
            curr_bin = sortedbins[-1]
        else:
            left_bin = sortedbins[idx-1]
            right_bin = sortedbins[idx]

            if size - left_bin[3] < right_bin[2] - size:
                # belongs tp left_bin
                curr_bin = left_bin
            else:
                curr_bin = right_bin

        peaks.loc[i, 'BIN'] = curr_bin[0]


def bin_stats(peaks):

    b = {}
    for i in peaks.index:
        bin_idx = int(peaks.loc[i, 'BIN'])
        size = peaks.loc[i, 'SIZE']
        try:
            b[bin_idx].values.append(size)
        except KeyError:
            b[bin_idx] = BinContainer(size=bin_idx, values = [ size ])

    return b


def adjust_bins(bins, stat, reset=False, repeats=-1):

    if reset:
        reset_bins(bins, stat, repeats)
        return

    for idx in range(len(bins)):
        b = bins[idx]
        value = b[0]
        if value in stat:
            #print('update bins', value)
            s = stat[value]
            values = sorted(s.values)
            size = len(values)
            if size % 2 == 0:
                idx = int(size/2) - 1
                med = (values[idx] + values[idx+1])/2
                q1 = values[:idx]
                q2 = values[idx:]
            else:
                idx = int(size//2)
                med = values[idx]
                q1 = values[:idx]
                q2 = values[idx+1:]

            b[1] = float(med)
            percentiles = np.percentile( s.values, [10, 50, 90])
            if not reset and s.d() > 0.5:
                #b[2] = float(np.mean(q1))
                #b[3] = float(np.mean(q2))
                b[2] = float(percentiles[0])
                b[3] = float(percentiles[2])
            else:
                b[2] = b[1] - 0.5
                b[3] = b[1] + 0.5

            #percentiles = np.percentile( s.values, [20, 50, 80] )
            #b[1] = float(percentiles[1])
            #if not reset and s.d() > 1.0:
            #    b[2] = float(percentiles[0])
            #    b[3] = float(percentiles[2])
            #else:
            #    if len(s.values) > 1:
            #        avg = float(b[1] + np.mean(s.values)) / 2
            #    else:
            #        avg = b[1]
            #    b[2] = avg - 1.0
            #    b[3] = avg + 1.0
                #b[2] = b[1] - 1.0
                #b[3] = b[1] + 1.0

def reset_bin_item(a_bin):
    a_bin[2] = a_bin[1] - 0.5
    a_bin[3] = a_bin[1] + 0.5


def reset_bins(bins, stat, repeats):

    anchor_bin = None
    for s in stat.values():
        if anchor_bin is None:
            anchor_bin = s
            continue
        if s.f() > anchor_bin.f():
            anchor_bin = s
        elif s.f() == anchor_bin.f() and s.d() < anchor_bin.d():
            anchor_bin = s

    # create a dict of bins from original bins
    bin_d = {}
    for b in bins:
        bin_d[b[0]] = b

    # set anchor bin
    u_bin = bin_d[anchor_bin.size]
    reset_bin_item(u_bin)

    # going up for each bin
    base_size = u_bin[1]
    for i in range(anchor_bin.size + repeats, bins[-1][0]+1, repeats):
        b = bin_d[i]
        if b[1] - base_size < repeats - 0.75:
            b[1] = base_size + repeats - 0.5
            reset_bin_item(b)
        base_size = b[1]

    # going down for each bin
    base_size = u_bin[1]
    for i in range(anchor_bin.size - repeats, bins[0][0]-1, -repeats):
        b = bin_d[i]
        if base_size - b[1] < repeats - 0.75:
            b[1] = base_size - repeats + 0.5
            reset_bin_item(b)
        base_size = b[1]

