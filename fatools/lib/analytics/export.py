
import sys
import csv
from itertools import zip_longest
from pandas import pivot_table


class autostream(object):

    def __init__( self, outstream ):
        self.outstream = outstream
        if 'b' in self.outstream.mode:
            self.write = self.encoder
        else:
            self.write = self.outstream.write

    def encoder(self, buff):
        if type(buff) == str:
            self.outstream.write( buff.encode('ASCII') )
        else:
            self.outstream.write( buff )


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

    outstream.write('LABEL\tMARKER\tSAMPLE\tBIN\tSIZE\tHEIGHT\tRATIO\tRANK\n')

    for analytical_set in analytical_sets:

        allele_df = analytical_set.allele_df.df
        for t in allele_df.itertuples():
            (marker_id, sample_id, value, size, height, assay_id, ratio, rank) = t[1:]
            marker = dbh.get_marker_by_id(marker_id)
            sample = dbh.get_sample_by_id(sample_id)
            outstream.write('%s\t%s\t%s\t%d\t%f\t%d\t%f\t%d\n' %
                (analytical_set.label, marker.code, sample.code, value, size, height, ratio, rank))


def export_moidf(analytical_sets, dbh, outstream):
    """ export MoI dataframe to a file suitable for loading into
        R or Python's pandas
    """

    from fatools.lib.analytics.moi import calculate_moi

    # format: LABEL SAMPLE MOI MLOCI

    outstream.write('LABEL\tSAMPLE\tMOI\tMLOCI\n')

    for analytical_set in analytical_sets:
        moi_result = calculate_moi(analytical_set.allele_df)
        label = analytical_set.label
        for t in moi_result.sample_dist.itertuples():
            (sample_id, moi_number, mloci_number) = t
            if dbh:
                sample_code = dbh.get_sample_by_id(sample_id).code
            else:
                sample_code = str(sample_id)
            outstream.write('%s\t%s\t%d\t%d\n' %
                (label, sample_code, moi_number, mloci_number))


def export_arlequin(analytical_sets, dbh, outstream, recode=False):
    """ export MLGT to Arlequin format
        recode: whether to use population-spesific alleles
    """

    _ = [   '[Profile]',
            '  Title="MsAF exported data"',
            '  NbSamples=%d' % len(analytical_sets),
            '  DataType=MICROSAT',
            '  GenotypicData=0',
            '  GameticPhase=0',
            '  MissingData="?"',
            '  LocusSeparator=WHITESPACE',
            '',
    ]

    _ += [  '[Data]',
            '  [[Samples]]',
    ]

    for (idx, analytical_set) in enumerate(analytical_sets):

        _ += [  '    SampleName="%s"' % analytical_set.label,
                '    SampleSize=%d' % analytical_set.sample_set.N,
                '    SampleData={',
        ]

        for e in analytical_set.allele_df.mlgt.itertuples():
            if recode:
                _.append('    %d 1 ' % e[0] + ' '.join('%02d%03d' % (idx,x) for x in e[1:]))
            else:
                _.append('    %d 1 ' % e[0] + ' '.join('%03d' % x for x in e[1:]))

        _.append('    }')

    outstream.write( '\n'.join(_))


def export_demetics(analytical_sets, dbh, outstream):
    """ export genotype data for export_demetics
        individual population fragment.length locus
        (individual -> sample code, population -> label, fragment.length -> allele, locus -> marker)
        demetics deals better with non-space-containing label
    """

    outstream.write('individual\tpopulation\tfragment.length\tlocus\n')

    for analytical_set in analytical_sets:

        allele_df = analytical_set.allele_df.df
        for t in allele_df.itertuples():
            (marker_id, sample_id, value, size, height, assay_id, allele_id, ratio, rank) = t[1:]
            marker = dbh.get_marker_by_id(marker_id)
            sample = dbh.get_sample_by_id(sample_id)
            outstream.write('%s\t%s\t%d\t%s\n' %
                (sample.code, reformat_label(analytical_set.label), value, marker.code))


def export_flat(analytical_set, dbh, outstream):
    """ export MLGT from single analytical set to flat format, eg:
        SAMPLECODE ALLELE1 ALLELE2 ALLELE3 ALLELE4
    """

    out = autostream( outstream )

    for e in analytical_set.allele_df.mlgt.itertuples():
        row = [ str(e[0])]
        out.write(
            (' '.join( row + list( str(int(x)) for x in e[1:]) ) + '\n')
        )


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
    'moidf': export_moidf,
    'arlequin': export_arlequin,
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

def reformat_label(label):
    """ some software can't deal with space or non-alphanumeric characters """
    return label.replace(' ','').replace('|','_').replace('#','_').replace('/','_')

