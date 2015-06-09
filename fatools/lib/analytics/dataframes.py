
from pandas import DataFrame, pivot_table



class AlleleDataFrame(object):

    def __init__(self, dbh, sample_ids, marker_ids, params):
        self.sample_ids = sample_ids
        self.marker_ids = marker_ids
        self.params = params
        self._df = None
        self._dominant_allele_df = None


    def get_dataframe(self):
        if self._df is None:
            self._df = self.dbh.get_allele_dataframe(self.sample_ids,
                        self.marker_ids, self.params)
        return self._df


    def get_dominant_alleles(self):
        """ return Pandas dataframe of (marker_id, sample_id, bin, size, height) """
        if self._dominant_allele_df is None:
            df = self.get_dataframe()
            idx = df.groupby(['marker_id','sample_id'])['height'].transform(max) == df['height']
            self._dominant_allele_df = df[idx]
        return self._dominant_allele_df



