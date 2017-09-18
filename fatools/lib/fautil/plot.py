"""
Collection of functions to do assay plotting using matplotlib.
"""
import matplotlib.pyplot as plt

from fatools.lib.utils import cerr


def align_fsa(fsa):
    """
    Align fsa to prepare for size and retention time extraction from each allele.

    Input
    -----
    fsa: class of fsa
    import Params() from fatools.lib.params for parameter in fsa alignment

    Output
    ------
    fsa that has been aligned
    """
    from fatools.lib import params

    fsa.align(params.Params())


def determine_number_of_subplots(channels):
    """
    Prepare the subplots needed by getting number of channels.
    The subplots are x and y-axis bound, allowing to be zoomed at the same time.

    Input
    -----
    channels: attribute from fsa class containing list of fsa channel

    Output
    ------
    matplotlib.figure containing N rows of axes in a column.
    """
    return plt.subplots(len(channels), 1, 'all')


def colorize_wavelength(wavelength):
    """
    Find dye color by using wavelen2rgb.

    Input
    -----
    wavelength: int of a dye wavelength

    Output
    ------
    RGB value in 3-tuple divided by 100: (R, G, B)

    The division by 100 is necessary because matplotlib color parameter
    only accepts value from 0-1.
    """
    from fatools.lib.fautil.wavelen2rgb import wavelen2rgb

    return tuple([color / 100 for color in wavelen2rgb(wavelength)])


def get_size_rtime_rfu(channel):
    """
    Get size, retention time, and RFU from the align method of fsa.

    Input
    -----
    channel: a channel class from one of the channels attribute in fsa class

    Output
    ------
    size_rtime_rfu: 3-tuples of size, rtime, and RFU from alleles in channel

    Size with value '-1' are not included in the collection.
    """
    alleles = channel.alleles
    size_rtime_rfu = []

    if alleles == []:
        return size_rtime_rfu

    for allele in alleles:
        if allele.size == -1:
            continue
        size_rtime_rfu.append((allele.size, allele.rtime, allele.rfu))

    return size_rtime_rfu


def prepare_second_x_axis(channel_axes, size_rtime_rfu):
    """
    Create a second x-axis to indicate the size of alleles.

    Input
    -----
    channel_axes: the channel axis to be marked
    size_rtime_rfu: the data for marking the second x-axis

    Output
    ------
    channel_axes that have size markings (if available)
    """
    sizes = []
    rtimes = []
    for size, rtime, rfu in size_rtime_rfu:
        sizes.append(int(size))
        rtimes.append(rtime)

    second_x_axis = channel_axes.twiny()
    second_x_axis.set_xlim(channel_axes.get_xlim())
    second_x_axis.set_xticks(rtimes)
    second_x_axis.set_xticklabels(sizes, rotation='vertical', fontsize=8)

    return second_x_axis


def save_or_show(figure, plot_file):
    """
    Determine if the plot is to be saved or shown.

    Input
    -----
    figure: class of figure from matplotlib
    plot_file: location and file name for saving plot

    Output
    ------
    If plot_file is None, then show plot to user.
    If plot_file is supplied, then the plot is saved to file.
    If plot_file is PdfPages object, save to PdfPages object.
    """
    saving_params = {'dpi': 150}
    plt.tight_layout()
    # HACK: to reduce complexity, the figure size is set to
    # approximately 22" monitor (21.3:9)
    figure.set_size_inches(20, 9)
    try:
        plot_file.savefig(**saving_params)
    except AttributeError:
        if plot_file is not None:
            plt.savefig(plot_file, **saving_params)
        else:
            plt.show()
    finally:
        plt.close()


def do_plot(fsa, plot_file=None):
    """
    Plot an assay in a plot.

    Input
    -----
    fsa: class of fsa
    plot_file: path for saving plot to file

    Output
    ------
    a figure class ready to be saved/shown
    """
    channels = fsa.channels
    fig = plt.figure()

    for channel in channels:
        color = colorize_wavelength(channel.wavelen)
        plt.plot(channel.data, color=color, label=channel.dye)

    plt.legend(framealpha=0.5)
    plt.title(fsa.filename)
    save_or_show(fig, plot_file)


def do_split_plot(fsa, plot_file=None):
    """
    Plot an assay dye, in every subplot.

    Input
    -----
    fsa: class of fsa
    plot_file: path for saving plot to file

    Output
    ------
    a figure class ready to be saved/shown
    """
    align_fsa(fsa)
    channels = fsa.channels
    whole_fig, whole_axes = determine_number_of_subplots(channels)
    twiny_axes = []

    for channel_axes_num, channel in enumerate(channels):
        color = colorize_wavelength(channel.wavelen)
        channel_axes = whole_axes[channel_axes_num]
        channel_axes.plot(channel.data, color=color, label=channel.dye)
        channel_axes.legend(framealpha=0.5)

        size_rtime_rfu = get_size_rtime_rfu(channel)
        if size_rtime_rfu:
            max_rfu = max(p[2] for p in size_rtime_rfu) * 1.2
        else:
            max_rfu = max(channel.data) * 1.2
        channel_axes.set_ylim((0, max_rfu))

        second_x_axis = prepare_second_x_axis(channel_axes, size_rtime_rfu)
        twiny_axes.append(second_x_axis)

    for axes in whole_axes:
        axes.get_shared_x_axes().join(*twiny_axes)

    plt.suptitle(fsa.filename)
    save_or_show(whole_fig, plot_file)


def ladder_plot(args, fsa_list, dbh=None):

    # filter FSAs
    plotted_fsas = []
    for (fsa, idx) in fsa_list:
        score, rss, nladder = fsa.align()
        if score < args.score:
            plotted_fsas.append( (score, rss, nladder, fsa) )
        elif args.rss > 0 and rss > args.rss:
            plotted_fsas.append( (score, rss, nladder, fsa) )

    # sort FSAs by score (ascending) and rss (descending)
    plotted_fsas.sort( key = lambda k: (k[0], -k[1]) )

    # for all fsas in plotted_fsas, create a plot
    whole_fig, whole_axes = determine_number_of_subplots(plotted_fsas)
    for fsa_axis_num, fsa_item in enumerate(plotted_fsas):
        fsa = fsa_item[3]

        # continue to prepare the plot


    # save the plot


def file_handler(fsa_list):
    """
    Generate fsa file from list of fsa to plot.

    Input
    -----
    fsa_list: list containing tuples of fsa and index

    Output
    ------
    A generator object returning fsa file
    """
    for fsa, index in fsa_list:
        yield fsa


def prepare_multi_page_pdf(plot_file):
    """
    Preparing PdfPages object for plotting to multi page pdf.

    Input
    -----
    plot_file: string name of plot file

    Output
    ------
    pdf: PdfPages object with plot_file name
    """
    from matplotlib.backends.backend_pdf import PdfPages

    pdf = PdfPages(plot_file)

    return pdf


def command_block(args, fsas, plot_file):
    """
    Commands to be done.

    Input
    -----
    args: arguments namespace from argparse
    fsas: list of fsa files
    plot_file: PdfPages object, string, or None
    """
    for fsa in fsas:
        if args.plot:
            do_plot(fsa, plot_file)
        if args.split_plot:
            do_split_plot(fsa, plot_file)


def plot(args, fsa_list, dbh=None):
    """
    The main function to handle all plot arguments given.

    Input
    -----
    args: arguments namespace from argparse
    fsa_list: list containing tuples of fsa and index
    dbh: *reserved for database handling*

    Output
    ------
    Determine if a PdfPages object needs to be created, then
    passing the plot file to the commands.
    """
    if args.plot and args.split_plot and args.plot_file:
        cerr('W: --plot, --split-plot, and --plot-file are flagged')
        cerr('W: This will only save the --split-plot results if format is not pdf')

    plot_file = args.plot_file
    if plot_file is not None:
        plot_file_ext = plot_file.split('.')[-1]
        if plot_file_ext == 'pdf':
            plot_file = prepare_multi_page_pdf(plot_file)

    fsas = file_handler(fsa_list)
    try:
        with plot_file as pdf:
            command_block(args, fsas, pdf)
    except AttributeError:
        command_block(args, fsas, plot_file)
