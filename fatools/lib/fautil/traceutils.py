
from math import factorial
import numpy as np
import attr
from scipy import ndimage, signal, optimize


_TOPHAT_FACTOR = 0.01 #025   #05
_MEDWINSIZE = 299
_MEDMMSIZE = 1991

def b(txt):
    """ return a binary string aka bytes """
    return txt.encode('UTF-8')


def smooth_signal( raw_signal ):
    """ smooth signal using savitzky_golay algorithm """
    return savitzky_golay( raw_signal, 11, 7 )


def correct_baseline( signal ):
    """ use tophat morphological transform to correct for baseline """

    return ndimage.white_tophat(signal, None,
                np.repeat([1], int(round(signal.size * _TOPHAT_FACTOR)))
            )


@attr.s
class NormalizedTrace(object):
    signal = attr.ib()
    baseline = attr.ib()
    mma = attr.ib()
    mmb = attr.ib()


def func_mm(x, a, b):
    """ Michaelis Menten kinetics equation -
        this has a nice property of passing initial value and plateau
    """
    return a*x/(b+x)

def normalize_baseline( raw ):
    """ return mean, median, sd and smooth signal """

    median_line = signal.medfilt(raw, [_MEDWINSIZE])
    baseline = signal.savgol_filter( median_line, _MEDWINSIZE, 7)
    corrected_baseline = raw - baseline
    np.maximum(corrected_baseline, 0, out=corrected_baseline)
    smooth = correct_baseline( signal.savgol_filter(corrected_baseline, 11, 7) )

    # perform michaelis-menten equation for baseline assessment
    mm_line = signal.medfilt(raw, [ _MEDMMSIZE ])
    xx = np.linspace(0, len(raw) )
    try:
        popt, pcov = optimize.curve_fit(func_mm, xx, mm_line)
    except:
        # michaelis menten are not appropriate for this scale
        popt = pcov = [0, 0]

    return NormalizedTrace( signal=smooth, baseline = baseline, mma = popt[0], mmb = popt[1] )


def search_peaks( signal, cwt_widths, min_snr ):
    """ returns [ (peak, height, area) ], ... ] """

    # find all peaks by cwt-based algorithm
    indices = find_peaks_cwt( signal, cwt_widths, min_snr=min_snr )

    if not indices:
        return []

    # find absolute heights

    raw_peaks = []
    for idx in indices:
        for i in range(3, -1, -1):
            try:
                height, index = max( [ (signal[i], i) for i in range(idx-3, idx+3) ] )
            except IndexError:
                continue
            break
        if height < params.min_height:
            continue
        if index < 0:
            continue
        raw_peaks.append( (index, height) )

    if not raw_peaks:
        return []


    # calculate area
    peaks = []
    for (peak, height) in raw_peaks:
        area, brtime, ertime = calculate_area( signal, peak, 5e-2 )
        peaks.append( (peak, height, area, brtime, ertime) )


    return peaks


def calculate_area(y, t, threshold):
    """ return (area, begin_rtime, end_rtime)
    """

    # right area
    data = y[t:]
    r_area, ertime, r_shared = half_area(data, threshold)

    # left area
    data = y[:t+1][::-1]
    l_area, brtime, l_shared= half_area(data, threshold)

    return ( l_area + r_area - y[t], t - brtime, ertime + t )



def half_area(y, threshold):
    """ return (area, ertime)
    """

    winsize = 3
    threshold = threshold/2
    shared = False
    area = y[0]
    edge = float(np.sum(y[0:winsize]))/winsize
    old_edge = 2 * edge

    index = 1
    limit = len(y)

    while edge > area * threshold and edge < old_edge and index < limit:
        old_edge = edge
        area += y[index]
        edge = float(np.sum(y[index:index+winsize]))/winsize
        index += 1
    if edge >= old_edge:
        shared = True
    index -= 1

    return area, index, shared


@attr.s
class TraceChannel(object):
    dye_name = attr.ib()
    dye_wavelength = attr.ib()
    raw_channel = attr.ib()
    smooth_channel = attr.ib()
    median = attr.ib()
    mean = attr.ib()
    sd = attr.ib()
    max_height = attr.ib()
    min_height = attr.ib()


def separate_channels( trace ):
    # return a list of [ 'dye name', dye_wavelength, numpy_array, numpy_smooth_baseline ]

    from fatools.lib.fautil.traceio import WAVELENGTH

    results = []
    for (idx, data_idx) in [ (1,1), (2,2), (3,3), (4,4), (5,105) ]:
        try:
            dye_name = trace.get_data(b('DyeN%d' % idx)).decode('UTF-8')
            #dye_wavelength = trace.get_data(b('DyeW%d' % idx))

            # below is to workaround on some strange dye names
            if dye_name == '6FAM': dye_name = '6-FAM'
            elif dye_name == 'PAT': dye_name = 'PET'
            elif dye_name == 'Bn Joda': dye_name = 'LIZ'

            try:
                dye_wavelength = trace.get_data(b('DyeW%d' % idx))
            except KeyError:
                dye_wavelength = WAVELENGTH[dye_name]

            raw_channel = np.array( trace.get_data(b('DATA%d' % data_idx)) )
            nt = normalize_baseline( raw_channel )

            results.append( TraceChannel(dye_name, dye_wavelength, raw_channel,
                            nt.signal, np.median(nt.baseline), np.mean(nt.baseline),
                            np.std(nt.baseline), np.max(nt.baseline), np.min(nt.baseline))
            )
        except KeyError:
            pass

    return results


def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError as msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')


def smooth(x,window_len=11,window='hanning'):
    """smooth the data using a window with requested size.

    This method is based on the convolution of a scaled window with the signal.
    The signal is prepared by introducing reflected copies of the signal
    (with the window size) in both ends so that transient parts are minimized
    in the begining and end part of the output signal.

    input:
        x: the input signal
        window_len: the dimension of the smoothing window; should be an odd integer
        window: the type of window from 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'
            flat window will produce a moving average smoothing.

    output:
        the smoothed signal

    example:

    t=linspace(-2,2,0.1)
    x=sin(t)+randn(len(t))*0.1
    y=smooth(x)

    see also:

    numpy.hanning, numpy.hamming, numpy.bartlett, numpy.blackman, numpy.convolve
    scipy.signal.lfilter

    TODO: the window parameter could be the window itself if an array instead of a string
    NOTE: length(output) != length(input), to correct this: return y[(window_len/2-1):-(window_len/2)] instead of just y.
    """

    if x.ndim != 1:
        raise ValueError( "smooth only accepts 1 dimension arrays." )

    if x.size < window_len:
        raise ValueError( "Input vector needs to be bigger than window size." )


    if window_len<3:
        return x


    if not window in ['flat', 'hanning', 'hamming', 'bartlett', 'blackman']:
        raise ValueError( "Window is on of 'flat', 'hanning', 'hamming', 'bartlett', 'blackman'")


    s=np.r_[x[window_len-1:0:-1],x,x[-1:-window_len:-1]]
    #print(len(s))
    if window == 'flat': #moving average
        w=np.ones(window_len,'d')
    else:
        w=eval('np.'+window+'(window_len)')

    y=np.convolve(w/w.sum(),s,mode='valid')
    return y

