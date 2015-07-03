
from collections import defaultdict
from fatools.lib.analytics.dataframes import AlleleDataFrame
from pandas import pivot_table
from pprint import pprint

class AnalyticalSet(object):
    """ AnalyticalSet

        _allele_df: dataframe of (marker_id, sample_id, value, size, height, assay_id)
        _marker_df:
    """

    def __init__(self, sample_set, params, marker_ids, dbh):
        assert sample_set and params and dbh

        self._sample_set = sample_set
        self._params = params

        if not marker_ids:
            marker_ids = params.get_marker_ids(dbh)
        self._marker_ids = marker_ids

        self._allele_df = AlleleDataFrame(dbh,
                    sample_ids = self.sample_ids,
                    marker_ids = marker_ids,
                    params = params
            )

        # placeholder

        self._marker_df = None

        self._filtered_sample_ids = None
        self._filtered_marker_ids = None
        self._sample_marker = None

        self._sample_genotyped_dist = None
        self._marker_genotyped_dist = None


    @property
    def label(self):
        return self._sample_set.label

    @property
    def sample_set(self):
        return self._sample_set

    @property
    def marker_ids(self):
        return self._allele_df.marker_ids

    @property
    def sample_ids(self):
        return self._sample_set.sample_ids

    @property
    def allele_df(self):
        return self._allele_df

    @property
    def marker_df(self):
        """ return a dataframe of:
                sample_id   marker_id1  marker_id2  marker_id3
                1           2           1           0

            where value = number of alleles in the marker for this particular sample
        """
        if self._marker_df is None:
            self._marker_df = pivot_table( self._allele_df.df,
                    index = 'sample_id', columns = 'marker_id', values='value', aggfunc = len )
        return self._marker_df

    @property
    def sample_marker(self):
        if self._sample_marker is None:
            self._sample_marker = {}
            for (idx, marker_id, sample_id, value, size, height, assay_id) in self.allele_df.df.itertuples():
                try:
                    self._sample_marker[sample_id].add( marker_id )
                except KeyError:
                    self._sample_marker[sample_id] = { marker_id }   # create initial set
        return self._sample_marker
    

    def get_filtered_sample_ids(self):
        """ get sample_ids that passed sample quality assessment """
        if self._filtered_sample_ids is None:
            self._filtered_sample_ids, self._sample_genotyped_dist = self.assess_sample_quality()
        return self._filtered_sample_ids


    def get_filtered_marker_ids(self):
        """ get marker_ids that passed marker quality assessment """
        if self._filtered_marker_ids is None:
            self._filtered_marker_ids, self._marker_genotyped_dist = self.assess_marker_quality()
        return self._filtered_marker_ids


    def get_sample_genotyped_distribution(self):
        if self._sample_genotyped_dist is None:
            self.get_filtered_sample_ids()
        return self._sample_genotyped_dist


    def get_marker_genotyped_distribution(self):
        if self._marker_genotyped_dist is None:
            self.get_filtered_marker_ids()
        return self._marker_genotyped_dist


    def assess_sample_quality(self, sample_qual_threshold = -1):
        """ assess sample based on successful genotyped markers
            param: sample_qual_threshold
        """

        sample_quality = [ (s_id, len(m)) for s_id, m in self.sample_marker.items() ]
        genotyped_dist = [ x[1] for x in sample_quality ]
        #n = max(genotyped_dist)
        n = len(self.marker_ids)
        if sample_qual_threshold < 0:
            sample_qual_threshold = self._params.sample_qual_threshold
        threshold = n * sample_qual_threshold
        passed_sample_ids = set([ int(x[0]) for x in sample_quality if x[1] >= threshold ])
        #failed_samples = len(sample_quality) - len(passed_sample_ids)

        return (passed_sample_ids, genotyped_dist)


    def assess_marker_quality(self, marker_qual_threshold = -1):
        """ assess marker based on successful genotyped samples, which must be done on
            all samples!!
            param: marker_qual_threshold
        """
        marker_genotyped = []
        for marker_id in self.marker_ids:
            # check of any x > 0 for marker_df[marker_id] = [ 2 1 0 0 0 ]
            genotyped = 0
            for m in self.marker_df[marker_id]:
                if m > 0:
                    genotyped += 1
            marker_genotyped.append( (marker_id, genotyped) )

        if marker_qual_threshold < 0:
            marker_qual_threshold = self._params.marker_qual_threshold
        threshold = len(self.sample_ids) * marker_qual_threshold
        passed_marker_ids = set([ x[0] for x in marker_genotyped if x[1] >= threshold ])
        return (passed_marker_ids, marker_genotyped)


    def get_filtered_analytical_set(self, sample_ids=None, marker_ids=None):

        if not (sample_ids or marker_ids):
            return None

        raise NotImplementedError()


class AnalyticalSetContainer(list):

    def __init__(self, sample_sets, params, marker_ids, dbh):
        super().__init__()
        self._sample_sets = sample_sets
        self._params = params
        for s in self._sample_sets:
            self.append( AnalyticalSet( s, params, marker_ids, dbh ) )


    def get_sample_sets(self):
        return self._sample_sets
            

    def get_filtered_sample_ids(self):
        """ return a filtered sample_ids, ie sample with minimum number of marker """
        sample_ids = set()
        for s in self:
            sample_ids.update( s.get_filtered_sample_ids() )
        return sample_ids


    def get_filtered_marker_ids(self):
        """ return a filtered marker_ids from total of all analytical sets """
        marker_counts = self.assess_marker_quality()
        threshold = self.total_samples * self._params.marker_qual_threshold
        marker_ids = [ x[0] for x in marker_counts.items() if x[1] > threshold ]
        return marker_ids


    def assess_marker_quality(self):
        """ assess marker quality """
        # collect all necessary data from each analyticalset
        marker_counts = defaultdict( int )
        for s in self:
            for (marker_id, count) in s.get_marker_genotyped_distribution():
                marker_counts[marker_id] += count

        return marker_counts
        

    def get_filtered_analytical_sets(self):
        """ return a filtered analytical set container """
        raise NotImplementedError


    @property
    def total_samples(self):
        return self._sample_sets.total_samples


def get_analytical_sets(dbh, sample_sets, params, marker_ids=None):

    assert sample_sets and params
    sets = AnalyticalSetContainer( sample_sets, params, marker_ids, dbh )

    return sets
        



