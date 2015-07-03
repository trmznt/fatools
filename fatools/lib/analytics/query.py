
import yaml
from fatools.lib.utils import cerr, cout
from fatools.lib.analytics.selector import Selector, Filter
from fatools.lib.analytics.analyticalset import get_analytical_sets
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

