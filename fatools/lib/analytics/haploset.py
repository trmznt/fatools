#
# haploset.py
#

class HaploSet(object):

    def __init__(self, analytical_set):
        self._analytical_set = analytical_set

        self._haplotype_df = analytical_set.allele_df.mlgt.dropna()
        self._sample_ids = set( self._haplotype_df.index )


    @property
    def haplotype_df(self):
        return self._haplotype_df


    @property
    def total_samples(self):
        return len(self._sample_ids)


    @property
    def sample_ids(self):
        return self._sample_ids



class HaploSetContainer(list):

    def __init__(self, analytical_sets):
        super().__init__()
        self._analytical_sets = analytical_sets
        for s in self._analytical_sets:
            self.append( HaploSet(s) )

    @property
    def total_samples(self):
        return sum( s.total_samples for s in self )


def get_haplotype_sets(analytical_sets):

    assert analytical_sets
    sets = HaploSetContainer(analytical_sets)

    return sets





