
from pandas import DataFrame

class base_sqlhandler(object):
    """ base class for SQLAlchemy-friendly handler """

    def __init__(self):
        self.engine = None
        self.session = None


    ## constructor for root classes

    def new_panel(self):
        p = self.Panel()
        p._dbh_session_ = self.session

    def new_marker(self):
        p = self.Marker()
        p._dbh_session_ = self.session

    def new_batch(self):
        p = self.Batch()
        p._dbh_session_ = self.session


    ## getter for single root classes

    def get_panel(self, panel_code):
        assert panel_code
        return self.Panel.search(panel_code, self.session)

    def get_batch(self, batch_code):
        assert batch_code
        return self.Batch.search(batch_code, self.session)

    def get_marker(self, marker_code):
        assert marker_code
        return self.Marker.search(marker_code, self.session)

    def get_by_id(self, class_, id):
        return class_.get(id, self.session)


    ## getter for multi root classes

    def get_markers(self):
        return self.Marker.query(self.session)

    def get_panels(self):
        pass

    def get_batches(self):
        pass

    def get_by_ids(self):
        pass


    ## getter for data

    def get_allele_dataframe(self, sample_ids, marker_ids, params):
        """ return a Pandas dataframe with this columns
            ( marker_id, sample_id, bin, size, height, assay_id )
        """

        # params ->
        # allele_absolute_threshold
        # allele_relative_threshold
        # allele_relative_cutoff
        # peak_type

        assert sample_ids and marker_ids and params

        q = self.session.query( self.AlleleSet.sample_id, self.Channel.assay_id,
                self.Allele.marker_id, self.Allele.bin,
                self.Allele.size, self.Allele.height
            ).join(self.Allele).join(self.Channel)

        q = q.filter( self.AlleleSet.sample_id.in_( sample_ids ) )

        if type(params.peaktype) in [ list, tuple ]:
            q = q.filter( self.Allele.type.in_( params.peaktype  ) )
        else:
            q = q.filter( self.Allele.type == params.peaktype )

        # we order based on marker_id, sample_id and then descending height
        q = q.order_by( self.Allele.marker_id, self.AlleleSet.sample_id,
                        self.Allele.height.desc() )

        print('MARKER IDS:', marker_ids)

        if marker_ids:
            q = q.filter( self.AlleleSet.marker_id.in_( marker_ids ) )
            #q = q.outerjoin( Marker, Allele.marker_id == Marker.id )
            #q = q.filter( Marker.id.in_( marker_ids ) )

        if params.abs_threshold > 0:
            q = q.filter( self.Allele.height > params.abs_threshold )

        if params.rel_threshold == 0 and params.rel_cutoff == 0:
            df = DataFrame( [ (marker_id, sample_id, value, size, height, assay_id )
                    for ( sample_id, assay_id, marker_id, value, size, height ) in q ] )

        else:

            alleles = []

            max_height = 0
            last_marker_id = 0
            last_sample_id = 0
            skip_flag = False
            for ( sample_id, assay_id, marker_id, value, size, height ) in q:
                if sample_id == last_sample_id:
                    if last_marker_id == marker_id:
                        if skip_flag:
                            continue
                        ratio = height / max_height
                        if ratio < params.rel_threshold:
                            continue
                        if (    params.rel_cutoff > 0 and
                                ratio > params.rel_cutoff ):
                            # turn off this marker by skipping this sample_id & marker_id
                            skip_flag = True
                            # don't forget to remove the latest allele
                            del alleles[-1]
                            continue
                            
                else:
                    last_sample_id = sample_id
                    last_marker_id = marker_id
                    max_height = height
                    skip_flag = False

                alleles.append( (marker_id, sample_id, value, size, height, assay_id) )

            df = DataFrame( alleles )
        
        if len(df) == 0:
            return df

        df.columns = ( 'marker_id', 'sample_id', 'value', 'size', 'height', 'assay_id' )
        return df





