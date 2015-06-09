
class SampleSelector(object):

    def __init__(self, samples = []):
        self.samples = samples
        self._sample_sets = None

    @staticmethod
    def from_dict(d):
        selector = SampleSelector()
        selector.samples = d
        return selector

    def to_dict(self):
        return { 'samples': self.samples }


    @staticmethod
    def load(yaml_text):
        d = yaml.load( yaml_text )
        selector = Selector.from_dict( d )
        return selector

    def dump(self):
        d = self.to_dict()
        return yaml.dump( d )


    def get_sample_ids(self, db):
        """ return sample ids; db is SQLa dbsession handler """
        pass


    def spec_to_sample_ids(self, spec_list):

        sample_ids = []
        for spec in spec_list:

            if 'query' in spec:
                    
                if '$' in spec['query']:
                    raise RuntimeError('cannot process differentiating query')

                if 'batch' in spec:
                    query = spec['batch'] + '[batch] & (' + spec['query'] + ')'

                sample_ids += query2list( parse_querycmd( query ) )

            elif 'codes' in spec:

                batch = Batch.search(spec['batch'])
                q = dbsession.query( Sample.id ).join( Batch ).filter( Batch.id == batch.id).filter( Sample.code.in_( spec['codes'] ) ) 

                sample_ids += list( q )

            elif 'ids' in spec:
                sample_ids += spec['ids']

            else:
                raise RuntimeError('sample spec format is incorrect')

        return sample_ids


    def get_sample_sets(self, db=None, sample_ids=None):

        if not self._sample_sets:
            
            if not db:
                db = dbsession
        
            sample_set = []
            colours = cycle( colour_list )

            if type(self.samples) == list:
                sample_set.append( SampleSet(   location='', year=0,
                        label = '-',
                        colour = 'blue',
                        sample_ids = self.spec_to_sample_ids( self.samples ) ) )

            elif type(self.samples) == dict:

                for label in self.samples:
            
                    sample_set.append( SampleSet( location='', year=0,
                                    label = label, colour = next(colours),
                                    sample_ids = self.spec_to_sample_ids(self.samples[label])))
                        

            self._sample_sets = sample_set

        if sample_ids:
            # filter based on sample_ids
            pass

        return self._sample_sets



class Parameter(object):

    def __init__(self):
        self.markers = []
        self.marker_ids = None
        self.species = None
        self.abs_threshold = 0
        self.rel_threshold = 0
        self.rel_cutoff = 0
        self.sample_qual_threshold = 0
        self.marker_qual_threshold = 0
        self.sample_options = None


    @staticmethod
    def from_dict(d):
        params = Parameter()
        params.markers = d.get('markers', None)
        params.marker_ids = d.get('marker_ids', None)
        params.abs_threshold = int( d['abs_threshold'] )
        params.rel_threshold = float( d['rel_threshold'] )
        params.rel_cutoff = float( d['rel_cutoff'] )
        params.sample_qual_threshold = float( d['sample_qual_threshold'] )
        params.marker_qual_threshold = float( d['marker_qual_threshold'] )
        params.sample_option = d['sample_option']
        return params


    def get_marker_ids(self):
        """ return marker ids;  """
        # self.markers is name
        if self.marker_ids is None and self.markers:
            markers = [ Marker.search(name) for name in self.markers ]
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




class Differentiation(object):
    """ MODE:
        spatial ~  -1 no differentiation
                    0 country
                    1 admin level 1
                    2 admin level 2
                    3 admin level 3
    """

    def __init__(self):
        self.spatial = 0
        self.temporal = 0
        self.differentiation = 0

    @staticmethod
    def from_dict(d):
        differentiation = Differentiation()
        differentiation.spatial = d['spatial']
        differentiation.temporal = d['temporal']
        differentiation.detection = d['detection']
        return differentiation


    def to_dict(self):
        pass

    @staticmethod
    def load(yaml_text):
        pass

    def dump(self):
        pass


    def create_groups(self, sample_sets):
        """ return a new sample set based on original sample sets and differentiation """

        sets = []
        sample_dfs = None
        colours = cycle( colour_list )

        for sample_set in sample_sets:

            sample_df = SampleDF.get_by_sample_ids( sample_set.get_sample_ids(),
                            self.spatial )

            samples = {}   # group sets

            for idx, sample_id, location, year, month, passive_detection in sample_df.itertuples():
        
                if self.temporal == 1:
                    tag = (location, year)
                elif self.temporal == 3:
                    quarter = 'Q1'
                    if month >= 9:
                        quarter = 'Q4'
                    elif month >= 6:
                        quarter = 'Q3'
                    elif month >= 3:
                        quarter = 'Q2'
                    tag = (location, '%d %s' % (year, quarter)) 
                else:
                    tag = (location, None)

                if self.detection:
                    (location, year) = tag
                    tag = (location, year, passive_detection)

                try:
                    samples[tag].append( int(sample_id) )
                except KeyError:
                    samples[tag] = [ int(sample_id) ]


            for tag in sorted(samples.keys()):
                if self.detection:
                    (location, year, passive_detection) = tag
                else:
                    (location, year) = tag
                    passive_detection = None

                label = sample_set.get_label() if not (location or year) else ''

                if passive_detection is not None:
                    extra_label = 'PD' if passive_detection else 'AD'
                else:
                    extra_label = None

                sets.append( SampleSet( location, year, next(colours), samples[tag],
                    label = label, extra_label = extra_label ))



            if sample_dfs is None:
                sample_dfs = sample_df
            else:
                sample_dfs = sample_dfs.append( sample_df, ignore_index=True )

        sets.sort( key = lambda x: x.get_label() )
        return (sets, sample_dfs)

