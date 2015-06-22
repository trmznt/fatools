
def correlate_alleles( set_1, set_2, marker_id ):

    sample_dict_1 = make_sample_allele_dict( set_1, marker_id )
    sample_dict_2 = make_sample_allele_dict( set_2, marker_id )

    sample_codes = set(sample_dict_1.keys()).intersection( set(sample_dict_2.keys()) )



def make_sample_allele_dict( a_set, marker_id ):

    pass
