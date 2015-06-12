

class SampleSet(object):

    def __init__(self, sample_ids, label=None, colour=None):
        """ sample_ids is list-like object """
        self.label = label
        self.sample_ids = sample_ids
        self.colour = colour

    def __repr__(self):
        return '<SampleSet: %s with %d sample(s)' % (self.label, len(self.sample_ids))

    @property
    def N(self):
        return len( self.sample_ids )




class SampleSetContainer(list):

    def append(self, sampleset):
        # check whether there's overlap sample ids
        for curr_set in self:
            if curr_set.sample_ids.intersection( sampleset.sample_ids ):
                raise RuntimeError( 'sample_id in container is not unique' )
        super().append( sampleset )

    @property
    def total_samples(self):
        N = 0
        for s in self:
            N += s.N
        return N


