
import yaml
from fatools.lib.analytics.selector import Selector, Filter


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


    def get_sample_sets(self, sample_ids = None):
        if self._sample_sets is None or sample_ids:
            selector = self._params['selector']
            self._sample_sets = selector.get_sample_sets(dbh, sample_ids)
        return self._sample_sets


    def get_analytical_sets(self, sample_ids = None):
        if self._analytical_sets is None or sample_ids:
            sample_sets = self.get_sample_sets(self._dbh, sample_ids)
            self._analytical_sets = get_filtered_analytical_sets( sample_sets,
                                        self._params['filter'] )
        return self._analytical_sets

