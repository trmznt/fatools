

class SampleSet(object):

    def __init__(self, sample_ids, label=None, colour=None):
        """ sample_ids is list-like object """
        assert type(sample_ids) is set, "sample_ids must be a set type"
        assert all(isinstance(x, int) for x in sample_ids), "all sample_ids member must have type int"
        self.label = label
        self.sample_ids = sample_ids
        self.colour = colour


    def __repr__(self):
        return '<SampleSet: %s with %d sample(s)' % (self.label, len(self.sample_ids))

    def filtered(self, sample_ids):
        assert type(sample_ids) is set, "sample_ids must be a set type"
        filtered_sample_ids = self.sample_ids.intersection( sample_ids )
        #print('filtered_sample_ids:', filtered_sample_ids)
        s = self.__class__(sample_ids = filtered_sample_ids, label = self.label,
                            colour = self.colour)
        return s

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


    def filtered(self, sample_ids):
        """ return a new SampleSetContainer with only sample_ids
        """
        assert type(sample_ids) is set, "sample_ids mut be a set type"
        s = self.__class__()
        for sample_set in self:
            s.append( sample_set.filtered(sample_ids) )
        return s

    @property
    def total_samples(self):
        N = 0
        for s in self:
            N += s.N
        return N




