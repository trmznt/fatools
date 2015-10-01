
from pandas import DataFrame
from sqlalchemy.orm.exc import NoResultFound

class base_sqlhandler(object):
    """ base class for SQLAlchemy-friendly handler """

    Panel = None
    Marker = None
    Batch = None
    Sample = None
    Assay = None

    def __init__(self):
        self.engine = None
        self.session = None


    ## constructor for root classes

    def new_panel(self):
        i = self.Panel()
        i._dbh_session_ = self.session()
        return i

    def new_marker(self):
        i = self.Marker()
        i._dbh_session_ = self.session()
        return i

    def new_batch(self):
        i = self.Batch()
        i._dbh_session_ = self.session()
        return i


    ## getter for single root classes

    def get_panel(self, panel_code):
        try:
            return self.Panel.search(panel_code, self.session)
        except NoResultFound:
            raise RuntimeError('Panel code %s does not exist!' % panel_code)

    def get_batch(self, batch_code):
        assert batch_code
        return self.Batch.search(batch_code, self.session)

    def get_marker(self, marker_code):
        assert marker_code
        return self.Marker.search(marker_code, self.session)

    def get_marker_by_id(self, id):
        return self.get_by_id(self.Marker, id)

    def get_batch_by_id(self, id):
        return self.get_by_id(self.Batch, id)

    def get_sample_by_id(self, id):
        return self.get_by_id(self.Sample, id)

    def get_assay_by_id(self, id):
        return self.get_by_id(self.Assay, id)

    def get_by_id(self, class_, id):
        assert class_ and id
        return class_.get(id, self.session)


    ## getter for multi root classes
    ## this will return a query object that can be further filtered by the caller

    def get_markers(self):
        return self.Marker.query(self.session)

    def get_panels(self):
        return self.Panel.query(self.session)

    def get_batches(self):
        return self.Batch.query(self.session)

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

        assert sample_ids
        assert marker_ids
        assert params

        q = self.session.query( self.AlleleSet.sample_id, self.Channel.assay_id,
                self.Allele.marker_id, self.Allele.bin,
                self.Allele.size, self.Allele.height
            ).join(self.Allele).join(self.Channel)

        q = q.filter( self.AlleleSet.sample_id.in_( sample_ids ) )

        q = self.customize_filter(q, params)

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


    def customize_filter(self, q, params):

        if type(params.peaktype) in [ list, tuple ]:
            q = q.filter( self.Allele.type.in_( params.peaktype  ) )
        else:
            q = q.filter( self.Allele.type == params.peaktype )

        return q

