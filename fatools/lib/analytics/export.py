
import sys
import csv
from itertools import zip_longest
from pandas import pivot_table


def export_major_tab(analytical_sets, dbh, outstream):

    output = []

    for analytical_set in analytical_sets:
        data, aux_data, assay_data = tabulate_data( analytical_set.allele_df.dominant_df,
                                            dbh )
        output.append( (analytical_set.label, data, aux_data, assay_data) )

    write_csv(output, outstream)

    return output


def export_tab(analytical_sets, dbh, outstream):

    output = []

    for analytical_set in analytical_sets:
        data, aux_data, assay_data = tabulate_data( analytical_set.allele_df.df, dbh )
        output.append( (analytical_set.label, data, aux_data, assay_data) )

    write_csv(output, outstream)

    return output


def export_major_r(analytical_sets, dbh, outstream):
    """ export to file suitable for loading into R
        the file will be tab-delimited and has header
    """

    output = []
    for analytical_set in analytical_sets:
        data, aux_data, assay_data = tabulate_data( analytical_set.allele_df.dominant_df,
                                            dbh )
        output.append( (analytical_set.label, data, aux_data, assay_data) )

    write_r(output, outstream)

    return output


def export_alleledf(analytical_sets, dbh, outstream):
    """ export allele dataframe to file suitable for loading
        into R or Python's pandas
    """

    # format: LABEL SAMPLE MARKER ALLELE SIZE HEIGHT AREA BETA THETA SYM SCORE TYPE

    df = None


def write_csv(output, outstream, delimiter='\t'):

    writer = csv.writer(outstream, delimiter=delimiter)
    for (label, rows, aux_rows, assay_rows) in output:
        writer.writerow( ('Label: %s' % label,))
        writer.writerows( rows )


export_format = {
    'major_tab': export_major_tab,
    'tab': export_tab,
    'major_r': export_major_r,
    'alleledf': export_alleledf,
}


def write_r(output, outstream, delimiter='\t'):

    # sanitity checking
    header = output[0][1][0]
    for (label, rows, aux_rows, assay_rows) in output[1:]:
        if header != rows[0]:
            raise RuntimeError('Headers between tabulated data do not match')

    writer = csv.writer(outstream, delimiter=delimiter)
    writer.writerow( ('Group',) + header )
    for (label, rows, aux_rows, assay_rows) in output:
        for row in rows[1:]:
            writer.writerow( (label, ) + row )


def export(analytical_sets, dbh, outfile, format='major_tab'):

    export_func = export_format[format]

    if outfile == '-':
        outstream = sys.stdout
    else:
        outstream = open(outfile, 'wt')

    export_func(analytical_sets, dbh, outstream)


def tabulate_data( allele_df, dbh ):

    buf = []
    buf2 = []
    buf3 = []

    table = pivot_table( allele_df,
                            index='sample_id', columns='marker_id', values='value',
                            aggfunc = lambda x: tuple(x) )

    heights = pivot_table( allele_df,
                            index='sample_id', columns='marker_id', values='height',
                            aggfunc = lambda x: tuple(x) )

    assay_ids = pivot_table( allele_df,
                            index='sample_id', columns='marker_id', values='assay_id',
                            aggfunc = lambda x: tuple(x) )

    buf.append( tuple( ['Sample', 'ID'] +
                    [ dbh.get_marker_by_id(x).code for x in table.columns ] ) )
    buf2.append( tuple( ['Sample', 'ID'] +
                    [ dbh.get_marker_by_id(x).code for x in heights.columns ] ) )
    buf3.append( tuple( ['Sample', 'ID'] +
                    [ dbh.get_marker_by_id(x).code for x in assay_ids.columns ] ) )


    empty = tuple()

    rows = [ ((dbh.get_sample_by_id(r[0]).code,), (r[0],)) + r[1:]
                for r in table.itertuples() ]
    rows.sort()

    height_rows = [ ((dbh.get_sample_by_id(r[0]).code,), (r[0],)) + r[1:]
                                for r in heights.itertuples() ]
    height_rows.sort()

    assayid_rows = [ ((dbh.get_sample_by_id(r[0]).code,), (r[0],)) + r[1:]
                                for r in assay_ids.itertuples() ]
    assayid_rows.sort()

    for row in rows:
        data = [ x if type(x) == tuple else empty for x in row ]
        for cols in zip_longest( *data, fillvalue='' ):
            buf.append( cols )

    for height_row in height_rows:
        data = [ x if type(x) == tuple else empty for x in height_row ]
        for cols in zip_longest( *data, fillvalue='' ):
            buf2.append( cols )

    for assayid_row in assayid_rows:
        data = [ x if type(x) == tuple else empty for x in assayid_row ]
        for cols in zip_longest( *data, fillvalue='' ):
            buf3.append( cols )


    return (buf, buf2, buf3)

