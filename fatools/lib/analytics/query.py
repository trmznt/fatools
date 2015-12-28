
import yaml
from fatools.lib.utils import cerr, cout
from fatools.lib.analytics.selector import Selector, Filter
from fatools.lib.analytics.analyticalset import get_analytical_sets
from fatools.lib.analytics.haploset import get_haplotype_sets
from pandas import DataFrame
from pprint import pprint

def load_yaml(yaml_text):

    d = yaml.load( yaml_text )
    instances = {}
    for k in d:
        if k == 'selector':
            instances['selector'] = Selector.from_dict( d[k] )
        elif k == 'filter':
            instances['filter'] = Filter.from_dict( d[k] )
        else:
            raise RuntimeError()

    return instances


class Query(object):

    def __init__(self, query_params, dbh):
        self._params = query_params
        self._dbh = dbh
        self._sample_sets = None
        self._analytical_sets = None
        self._filtered_sample_sets = None
        self._filtered_analytical_sets = None
        self._filtered_haplotype_sets = None


    def get_sample_sets(self, sample_ids = None):
        if self._sample_sets is None or sample_ids:
            selector = self._params['selector']
            self._sample_sets = selector.get_sample_sets(self._dbh, sample_ids)
        return self._sample_sets


    def get_analytical_sets(self, sample_ids = None):
        if self._analytical_sets is None or sample_ids:
            cerr('[query]: getting initial analytical sets')
            sample_sets = self.get_sample_sets( sample_ids )
            self._analytical_sets = get_analytical_sets( self._dbh, sample_sets,
                                        self._params['filter'] )
            cerr('[query]: initial total samples: %d' % self._analytical_sets.total_samples)
        return self._analytical_sets


    def get_filtered_sample_sets(self, sample_ids = None):
        if self._filtered_sample_sets is None or sample_ids:
            if not sample_ids:
                sample_ids = self.get_analytical_sets().get_filtered_sample_ids()
            self._filtered_sample_sets = self.get_sample_sets().filtered( sample_ids )
            cerr('[query]: filtered total samples: %d' %
                    self._filtered_sample_sets.total_samples)
        return self._filtered_sample_sets


    def get_filtered_analytical_sets(self, sample_ids = None):
        if self._filtered_analytical_sets is None or sample_ids:

            # get initial sample set, and filter the sample set by sample ids
            filtered_sample_sets = self.get_filtered_sample_sets(sample_ids)

            # create new analytical sets based on the filtered sample sets
            cerr('[query]: getting analytical sets with filtered sample sets')
            filtered_analytical_sets = get_analytical_sets( self._dbh, filtered_sample_sets,
                                        self._params['filter'] )
            cerr('[query]: filtered total samples: %d'
                        % filtered_analytical_sets.total_samples)

            # get filtered marker ids
            filtered_marker_ids = filtered_analytical_sets.get_filtered_marker_ids()
            cerr('[query]: filtered marker ids: %s' % str(filtered_marker_ids))

            # filter markers by retaining marker ids and removing others
            filtered_analytical_sets = get_analytical_sets( self._dbh, filtered_sample_sets,
                                        self._params['filter'],
                                        marker_ids = filtered_marker_ids )

            self._filtered_analytical_sets = filtered_analytical_sets

        return self._filtered_analytical_sets


    def get_filtered_haplotype_sets(self):
        if self._filtered_haplotype_sets is None:
            self._filtered_haplotype_sets = get_haplotype_sets(
                    self.get_filtered_analytical_sets())
        return self._filtered_haplotype_sets


    def get_sample_summary(self, mode='allele'):
        """ return a pandas dataframe containing (label, initial, filtered) headings
            since sample_summary does not involve analysis anymore, it is included here
        """

        if mode not in [ 'allele', 'mlgt' ]:
            raise RuntimeError('unknown get_sample_summary() mode: %s' % mode)

        # count initial sample set
        initial_samples = {}
        for sample_set in self.get_sample_sets():
            initial_samples[sample_set.label] = sample_set.N

        # count filtered sample set
        filtered_samples = {}
        for analytical_set in self.get_filtered_analytical_sets():
            filtered_samples[analytical_set.label] = analytical_set.N

        if not mode == 'mlgt':
            return DataFrame( { 'Initial Samples': initial_samples,
                            'Filtered Samples': filtered_samples },
                            columns = [ 'Initial Samples', 'Filtered Samples']
                )

        # count MLGs sample set
        mlg_samples = {}
        for haplotype_set in self.get_filtered_haplotype_sets():
            mlg_samples[haplotype_set.label] = haplotype_set.N

        return DataFrame( { 'Initial Samples': initial_samples,
                            'Filtered Samples': filtered_samples,
                            'MLG Samples': mlg_samples },
                            columns = [ 'Initial Samples', 'Filtered Samples', 'MLG Samples']
                )


