"""
Collection of functions to do assay plotting using matplotlib.
"""
import matplotlib.pyplot as plt


def align_fsa(fsa):
    """
    Align fsa to prepare for size and retention time extraction from each allele.

    Input

    fsa: class of fsa
    import Params() from fatools.lib.params for parameter in fsa alignment

    Output

    fsa that has been aligned
    """
    from fatools.lib import params

    fsa.align(params.Params())


def determine_number_of_subplots(channels):
    """
    Prepare the subplots needed by getting number of channels.
    The subplots are x and y-axis bound, allowing to be zoomed at the same time.

    Input

    channels: attribute from fsa class containing list of fsa channel

    Output

    matplotlib.figure containing N rows of axes in a column.
    """
    return plt.subplots(len(channels), 1, 'all')


def colorize_wavelength(wavelength):
    """
    Find dye color by using wavelen2rgb.

    Input

    wavelength: int of a dye wavelength

    Output

    RGB value in 3-tuple divided by 100: (R, G, B)

    The division by 100 is necessary because matplotlib color parameter
    only accepts value from 0-1.
    """
    from fatools.lib.fautil.wavelen2rgb import wavelen2rgb

    return tuple([color / 100 for color in wavelen2rgb(wavelength)])


def get_size_rtime(channel):
    """
    Get size and retention time from the align method of fsa.

    Input

    channel: a channel class from one of the channels attribute in fsa class

    Output

    size_rtime: pairs of size and rtime from alleles in channel

    Size with value '-1' are not included in the collection.
    """
    alleles = channel.alleles
    size_rtime = []

    if alleles == []:
        return size_rtime

    for allele in alleles:
        if allele.size == -1:
            continue
        size_rtime.append((allele.size, allele.rtime, allele.rfu))

    return size_rtime


def prepare_second_x_axis(channel_axis, size_rtime):
    """
    Create a second x-axis to indicate the size of alleles.

    Input

    channel_axis: the channel axis to be marked
    size_rtime: the data for marking the second x-axis

    Output

    channel_axis that have size markings (if available)
    """
    sizes = []
    rtimes = []
    for size, rtime, rfu in size_rtime:
        sizes.append(int(size))
        rtimes.append(rtime)

    second_x_axis = channel_axis.twiny()
    second_x_axis.set_xlim(channel_axis.get_xlim())
    second_x_axis.set_xticks(rtimes)
    second_x_axis.set_xticklabels(sizes, rotation='vertical', fontsize=8)


def do_split_plot(fsa, plot_file=None):
    """
    Plot an assay dye, in every subplot.

    Input

    fsa: class of fsa
    plot_file: path for saving plot to file

    Output

    If plot_file is None, then show plot to user.
    If plot_file is supplied, then the plot is saved to file.
    """
    align_fsa(fsa)
    channels = fsa.channels
    whole_fig, whole_axes = determine_number_of_subplots(channels)

    for channel_axis_num, channel in enumerate(channels):
        color = colorize_wavelength(channel.wavelen)
        channel_axis = whole_axes[channel_axis_num]
        channel_axis.plot(channel.data, color=color, label=channel.dye)
        channel_axis.legend(framealpha=0.5)

        size_rtime = get_size_rtime(channel)
        prepare_second_x_axis(channel_axis, size_rtime)
        if len(size_rtime) > 0:
            max_rfu = max( p[2] for p in size_rtime ) * 1.2
            channel_axis.set_ylim((0, max_rfu))

    plt.suptitle(fsa.filename)
    plt.tight_layout()

    if plot_file is not None:
        # HACK: to reduce complexity, the figure size is set to
        # approximately 22" monitor (16:19)
        whole_fig.set_size_inches(19, 11)
        plt.savefig(plot_file, dpi=150, bbox_inches='tight')
    else:
        plt.show()


def file_handler(fsa_list):
    """
    Generate fsa file from list of fsa to plot.

    Input

    fsa_list: list containing tuples of fsa and index

    Output

    A generator object returning fsa file
    """
    for fsa, index in fsa_list:
        yield fsa


def split_plot(args, fsa_list, dbh=None):
    """
    The main function to handle all arguments given.

    Input

    args: arguments namespace from argparse
    fsa_list: list containing tuples of fsa and index
    dbh: *reserved for database handling*

    Output

    Calling do_split_plot function for every fsa file passed
    """
    fsas = file_handler(fsa_list)
    for fsa in fsas:
        do_split_plot(fsa, args.plot_file)
