

class SampleSet(object):

    def __init__(self, sample_ids, label=None, colour=None):
        """ sample_ids is list-like object """
        self.label = label
        self.sample_ids = sample_ids
        self.colour = colour



class SampleSetContainer(list):

    def append(self, sampleset):
        # check whether there's overlap sample ids
        for curr_set in self:
            if curr_set.intersect_with( sampleset ):
                return False
        super().append( sampleset )


