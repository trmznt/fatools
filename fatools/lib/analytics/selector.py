
from itertools import cycle
from fatools.lib.analytics.sampleset import SampleSet, SampleSetContainer
from fatools.lib.const import peaktype

#colour_list = [ 'r', 'g', 'b' ]
colour_list = [ '#1f78b4','#33a02c','#e31a1c','#ff7f00','#6a3d9a','#b15928',
                '#a6cee3','#b2df8a','#fb9a99','#fdbf6f','#cab2d6','#ffff99']


class Selector(object):

    def __init__(self, samples = []):
        self.samples = samples
        self.global_options = None
        self._sample_sets = None

    @classmethod
    def from_dict(cls, d):
        selector = cls()
        if '_:_' in d:
            global_options = d.pop('_:_')
        else:
            global_options = {}
        selector.samples = d
        selector.global_options = global_options
        return selector

    def to_dict(self):
        return { 'samples': self.samples }


    @classmethod
    def load(cls, yaml_text):
        d = yaml.load( yaml_text )
        selector = cls.from_dict( d )
        return selector

    def dump(self):
        d = self.to_dict()
        return yaml.dump( d )


    def get_sample_ids(self, db):
        """ return sample ids; db is SQLa dbsession handler """
        pass


    def spec_to_sample_ids(self, spec_list, dbh, sample_ids=None):

        global_ids = set()
        for spec in spec_list:

            # initial spec
            ids = set()

            if 'query' in spec:

                if '$' in spec['query']:
                    raise RuntimeError('cannot process differentiating query')

                if 'batch' in spec:
                    query = spec['batch'] + '[batch] & (' + spec['query'] + ')'

                ids.update( query2set( parse_querycmd( query ) ) )

            elif 'codes' in spec:

                batch = dbh.get_batch(spec['batch'])
                ids.update( set(batch.get_sample_ids_by_codes( spec['codes'] )) )

            elif 'ids' in spec:
                ids.update( set(spec['ids']) )

            elif 'batch' in spec:
                batch = dbh.get_batch(spec['batch'])
                ids.update( batch.sample_ids )

            elif 'batch_id' in spec:
                batch = dbh.get_batch_by_id(spec['batch_id'])
                ids.update( batch.sample_ids )

            else:
                raise RuntimeError(
                    'sample spec format is incorrect, mandatory fields must exist!')

            # filtering spec

            q = dbh.session().query(dbh.Sample.id).filter(dbh.Sample.id.in_(ids))

            if 'category' in spec:
                q = q.filter(dbh.Sample.category == int(spec['category']))
            if 'int1' in spec:
                q = q.filter(dbh.Sample.int1 == int(spec['int1']))
            if 'int2' in spec:
                q = q.filter(dbh.Sample.int2 == int(spec['int2']))

            q = self.filter_sample(spec, dbh, q)

            ids = set(x.id for x in q)

            global_ids.update( ids )

        if sample_ids is not None:
            assert type(sample_ids) is set, "Please provide sample_ids as set"
            global_ids = global_ids.intersection( sample_ids )

        return global_ids


    def get_sample_sets(self, dbh, sample_ids=None):

        if not self._sample_sets:

            assert dbh, "dbh must be specified"

            sample_sets = SampleSetContainer()

            if type(self.samples) == list:
                sample_sets.append(
                    SampleSet(label = '-', colour = 'blue',
                        sample_ids = self.spec_to_sample_ids(self.samples, dbh, sample_ids)
                    )
                )

            elif type(self.samples) == dict:

                colours = cycle( colour_list )

                for label in self.samples:
                    sample_sets.append(
                        SampleSet(label = label, colour = next(colours),
                            sample_ids = self.spec_to_sample_ids(self.samples[label],
                                                        dbh, sample_ids)
                        )
                    )

            self._sample_sets = sample_sets

        return self._sample_sets


    def filter_sample(self, spec_list, dbh, q):
        """ please override this method as necessary"""
        return q


class Filter(object):


    def __init__(self):
        self.markers = []
        self.marker_ids = None
        self.species = None
        self.abs_threshold = 0      # includes alelles above rfu
        self.rel_threshold = 0.0      # includes alleles above % of highest rfu
        self.rel_cutoff = 0.0         # excludes alleles above % of highest rfu [ rel_threshold < height < rel_cutoff ]
        self.sample_qual_threshold = 0.0    # includes samples with marker more than %
        self.marker_qual_threshold = 0.0    # includes markers with sample more than %
        self.peaktype = peaktype.bin
        self.sample_options = None
        self.stutter_ratio = 0.0
        self.stutter_range = 3.5
        self.stutter_baserange = 4.5
        self.stutter_baseratio = 0.1


    @staticmethod
    def from_dict(d):
        params = Filter()
        params.markers = d.get('markers', None)
        params.marker_ids = d.get('marker_ids', None)
        params.abs_threshold = int( d['abs_threshold'] )
        params.rel_threshold = float( d['rel_threshold'] )
        params.rel_cutoff = float( d['rel_cutoff'] )
        params.sample_qual_threshold = float( d['sample_qual_threshold'] )
        params.marker_qual_threshold = float( d['marker_qual_threshold'] )
        params.sample_filtering = d['sample_filtering']
        params.stutter_ratio = float(d.get('stutter_ratio', 0))
        params.stutter_range = float(d.get('stutter_range', 3.5))
        return params


    def get_marker_ids(self, dbh=None):
        """ return marker ids;  """
        # self.markers is name
        if (self.marker_ids is None and self.markers) and dbh:
            # only execute below if dbh is provided, marker_ids is empty and
            # markers is not empty
            markers = [ dbh.get_marker(name) for name in self.markers ]
            self.marker_ids = [ marker.id for marker in markers ]
        return self.marker_ids


    def to_dict(self):
        pass


    @staticmethod
    def load(yaml_text):
        pass

    def dump(self):
        pass


    def get_analytical_sets(self, sample_sets, marker_ids=None ):

        sets = []
        if marker_ids == None:
            marker_ids = self.get_marker_ids()
        for sample_set in sample_sets:
            sets.append( sample_set.get_analytical_set( marker_ids = marker_ids,
                                allele_absolute_threshold = self.abs_threshold,
                                allele_relative_threshold = self.rel_threshold,
                                allele_relative_cutoff = self.rel_cutoff,
                                unique = (self.sample_options == 'U') ) )

        return sets
