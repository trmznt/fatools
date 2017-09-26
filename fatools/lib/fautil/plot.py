"""
Collection of functions to do assay plotting using matplotlib.
"""
import matplotlib.pyplot as plt

from fatools.lib.utils import cerr, cexit


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


def determine_figure_size(list_of_data):
    """
    Prepare the figure size needed by getting number of data.

    Input
    -----
    list_of_data: list for determining the data amount

    Output
    ------
    matplotlib.figure with size to accomodate subplots
    """
    # Every axes are given 2" in height
    height = 2 * len(list_of_data)
    figure = plt.figure()
    figure.set_size_inches(20, height)

    return figure


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
    for size, rtime, _ in size_rtime_rfu:
        sizes.append(int(size))
        rtimes.append(rtime)

    second_x_axis = channel_axes.twiny()
    second_x_axis.set_xlim(channel_axes.get_xlim())
    second_x_axis.set_xticks(rtimes)
    second_x_axis.set_xticklabels(sizes, rotation='vertical', fontsize=8)

    return second_x_axis


def save_or_show(plot_file):
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
    plt.tight_layout()
    try:
        plot_file.savefig(dpi=150)
    except AttributeError:
        if plot_file is not None:
            plt.savefig(plot_file, dpi=150)
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

    for channel in channels:
        color = colorize_wavelength(channel.wavelen)
        plt.plot(channel.data, color=color, label=channel.dye)

    plt.legend(framealpha=0.5)
    plt.title(fsa.filename)
    save_or_show(plot_file)


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
    figure = determine_figure_size(channels)
    fsa_subplots = []
    twiny_axes = []

    for channel_axes_num, channel in enumerate(channels):
        color = colorize_wavelength(channel.wavelen)
        channel_axes = figure.add_subplot(len(channels), 1, channel_axes_num + 1)
        channel_axes.plot(channel.data, color=color, label=channel.dye)
        channel_axes.legend(framealpha=0.5)
        fsa_subplots.append(channel_axes)

        size_rtime_rfu = get_size_rtime_rfu(channel)
        if size_rtime_rfu:
            max_rfu = max(p[2] for p in size_rtime_rfu) * 1.2
        else:
            max_rfu = max(channel.data) * 1.2
        channel_axes.set_ylim((0, max_rfu))

        second_x_axis = prepare_second_x_axis(channel_axes, size_rtime_rfu)
        twiny_axes.append(second_x_axis)

        if channel_axes_num == 0:
            channel_axes.set_title(fsa.filename, y=1.3)

    for axes in fsa_subplots:
        axes.get_shared_x_axes().join(*fsa_subplots)
        axes.get_shared_x_axes().join(*twiny_axes)

    save_or_show(plot_file)


def determine_included_fsa_to_plot(score, rss, fsas):
    """
    Separate based on score and RSS value.

    Input
    -----
    score: int/float score threshold for fsa exclusion
    rss: int/float RSS threshold for fsa exclusion
    fsas: list of fsa files

    Output
    ------
    included_fsas: list of fsa that is lower than the score threshold
                   set in --score parameter and higher than the RSS
                   threshold set in --rss parameter
    """
    included_fsas = []
    for fsa in fsas:
        align_fsa(fsa)
        if fsa.score < score:
            included_fsas.append(fsa)
        elif rss > 0 and fsa.rss > rss:
            included_fsas.append(fsa)

    # sort FSAs by score (ascending) and rss (descending)
    included_fsas.sort(key=lambda fsa: (fsa.score, -fsa.rss))
    # Limit to the top 100 worst fsa score & RSS
    included_fsas = included_fsas[:100]

    return included_fsas


def do_ladder_plot(fsas, plot_file):
    """
    Create a plot of the ladder channel from fsa files.

    Input
    -----
    args: arguments namespace from argparse
    fsas: list of fsa files
    plot_file: path for saving plot to file

    Output
    ------
    a figure class ready to be saved/shown
    """
    figure = determine_figure_size(fsas)

    for ladder_axes_num, fsa in enumerate(fsas):
        ladder = fsa.get_ladder_channel()
        color = colorize_wavelength(ladder.wavelen)
        ladder_axes = figure.add_subplot(len(fsas), 1, ladder_axes_num + 1)
        ladder_axes.plot(ladder.data, color=color, label=ladder.dye)
        ladder_axes.legend(framealpha=0.5)

        size_rtime_rfu = get_size_rtime_rfu(ladder)
        if size_rtime_rfu:
            max_rfu = max(p[2] for p in size_rtime_rfu) * 1.2
        else:
            max_rfu = max(ladder.data) * 1.2
        ladder_axes.set_ylim((0, max_rfu))

        prepare_second_x_axis(ladder_axes, size_rtime_rfu)

        title = '{filename} | Score: {score:.2f} | RSS: {rss:.2f}'.format(
            filename=fsa.filename,
            score=fsa.score,
            rss=fsa.rss)
        ladder_axes.set_title(title, y=1.3)

    save_or_show(plot_file)


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


def check_and_prepare_pdf(plot_file):
    """
    Check if format is supported by matplotlib, then determine if
    PdfPages object needs to be prepared for plotting to pdf.

    Input
    -----
    plot_file: string of plot file name and format

    Output
    ------
    plot_file: PdfPages object with plot_file name if format is '.pdf'
    """
    from os.path import splitext

    from matplotlib.backends.backend_pdf import PdfPages

    if plot_file is not None:
        plot_file_ext = splitext(plot_file)[-1]
        if plot_file_ext == '.pdf':
            plot_file = PdfPages(plot_file)
        else:
            try:
                plt.savefig(plot_file)
            except ValueError:
                cerr('E: Format {} is not supported!'.format(plot_file_ext))
                cexit('Exiting...')

    return plot_file


def command_block(args, fsas, plot_file):
    """
    Prepare the necessary data and hold the commands to be done.
    Give warnings if there are overlapping arguments.

    Input
    -----
    args: arguments namespace from argparse
    fsas: list of fsa files
    plot_file: PdfPages object, string, or None
    """
    if args.ladder_plot:
        if args.plot or args.split_plot:
            cerr('W: --plot, --split-plot and --ladder-plot are flagged')
            cerr('W: Only the --ladder-plot option will be done')

        included_fsas = determine_included_fsa_to_plot(args.score, args.rss, fsas)
        do_ladder_plot(included_fsas, plot_file)

        return

    elif args.plot or args.split_plot:
        if args.plot and args.split_plot and args.plot_file:
            cerr('W: --plot, --split-plot, and --plot-file are flagged')
            cerr('W: This will only save the --split-plot results if format is not in pdf')

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
    passing the args, fsas, and plot file to the commands.
    """
    plot_file = check_and_prepare_pdf(args.plot_file)

    try:
        with plot_file as pdf:
            command_block(args, file_handler(fsa_list), pdf)
    except AttributeError:
        command_block(args, file_handler(fsa_list), plot_file)
