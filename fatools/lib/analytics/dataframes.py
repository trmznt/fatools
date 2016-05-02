
from pandas import DataFrame, pivot_table



class AlleleDataFrame(object):

    def __init__(self, dbh, sample_ids, marker_ids, params):
        self.sample_ids = sample_ids
        self.marker_ids = marker_ids
        self.params = params
        self._df = dbh.get_allele_dataframe(self.sample_ids, self.marker_ids,
                        self.params)
        self._dominant_df = None

        # grouped dataframe based on [ 'marker_id', 'value']
        self._group_df = None
        self._group_dominant_df = None

        # distributions
        self._df_distribution = None
        self._dominant_df_distribution = None

        # genotypes and MLGT
        self._genotype_df = None
        self._mlgt_df = None
        self._unique_mlgt_df = None

        # moi
        self._allele_multiplicity = None
        self._sample_multiplicity = None
        self._locus_multiplicity = None


    @property
    def df(self):
        return self._df


    @property
    def dominant_df(self):
        """ return Pandas dataframe of (marker_id, sample_id, value, size, height) """
        if self._dominant_df is None:
            df = self.df
            idx = df.groupby(['marker_id','sample_id'])['height'].transform(max) == df['height']
            self._dominant_df = df[idx]
        return self._dominant_df


    @property
    def grouped_df(self):
        """ return Pandas dataframe grouped by ['marker_id', 'value']
            ie: [marker_id, value] = (marker_id, sample_id, value, size, height, assay_id)
        """
        if self._group_df is None:
            self._group_df = self.df.groupby( ['marker_id', 'value'] )
        return self._group_df


    @property
    def grouped_dominant_df(self):
        """ return Pandas dataframe for dominant alleles group by ['marker_id', 'value'] """
        if self._group_dominant_df is None:
            self._group_dominant_df = self.dominant_df.groupby(['marker_id', 'value'])
        return self._group_dominant_df


    @property
    def df_distribution(self):
        if self._df_distribution is None:
            self._df_distribution = pivot_table(self.df,
                    index = ['marker_id', 'value'],
                    values = 'sample_id',
                    aggfunc = len)
        return self._df_distribution


    @property
    def dominant_df_distribution(self):
        if self._df_distribution is None:
            self._dominant_df_distribution = pivot_table(self.dominant_df,
                    index = ['marker_id', 'value'],
                    values = 'sample_id',
                    aggfunc = len)
        return self._dominant_df_distribution


    @property
    def genotype_df(self):
        if self._genotype_df is None:
            self._genotype_df = pivot_table(self.df,
                    index = ['sample_id'],
                    columns = ['marker_id'],
                    values = ['value', 'height'],
                    aggfunc = tuple)
        return self._genotype_df


    @property
    def mlgt(self):
        """ return MLGT (Multi Locus GenoTyping) dataframe
            this will drop samples with NaN (or None) allele, hence please
            check the dataframe first
        """
        if self._mlgt_df is None:
            self._mltg_df = pivot_table(self.dominant_df,
                    index = 'sample_id',
                    columns = 'marker_id',
                    values = 'value',
                    dropna=False).dropna(how='any')
        return self._mltg_df


    @property
    def unique_mlgt(self):
        """ return unique MLGT dataframe (randomly select samples with
            identical MLGT)
        """
        if self._unique_mlgt_df is None:
            self._unique_mlgt_df = self.mlgt.drop_duplicates()
        return self._unique_mlgt_df


    @property
    def allele_multiplicity(self):
        if self._allele_multiplicity is None:
            self._allele_multiplicity = pivot_table(self.df,
                    index = ['sample_id'],
                    columns = 'marker_id',
                    values = 'value',
                    aggfunc = len,
                    fill_value=0)
        return self._allele_multiplicity


    @property
    def sample_multiplicity(self):
        if self._sample_multiplicity is None:
            # apply max to each row (axis=1)
            self._sample_multiplicity = self.allele_multiplicity.max(1)
        return self._sample_multiplicity


    @property
    def locus_multiplicity(self):
        if self._locus_multiplicity is None:
            self._locus_multiplicity = self.allele_multiplicity.applymap(lambda x: 1 if x > 1 else 0)
        return self._locus_multiplicity
