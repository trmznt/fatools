

## handler

#from fatools.lib.sqlmodels import schema
#from fatools.lib.utils import cout, cerr
#from pandas import DataFrame
import os, sys

class SQLHandler_XXX(object):


    #Marker = schema.Marker

    #Batch = schema.Batch

    def __init__(self, dbfile, initial=False):
        print("Opening db: %s" % dbfile)
        if not initial and not os.path.isfile(dbfile):
            cerr('ERR - sqlite db file not found: %s' % dbfile)
            sys.exit(1)
        if initial and os.path.isfile(dbfile):
            cerr('ERR - sqlite db file already exists: %s' % dbfile)
            sys.exit(1)
        self.dbfile = dbfile
        self.engine, self.session = schema.engine_from_file(dbfile)


    def initdb(self, create_table = True):
        if create_table:
            schema.Base.metadata.create_all(self.engine)
        from fatools.lib.sqlmodels.setup import setup
        setup( self.session )
        cout('Database at %s has been initialized.' % self.dbfile)


    def get_batch(self, batch_code = None):

        if not batch_code:
            cerr('ERR - batch code must be supplied!')
            sys.exit(1)

        return schema.Batch.search(batch_code, self.session)


    def get_marker(self, marker_code = None):

        if not marker_code:
            cerr('ERR - marker code must be supplied')
            sys.exit(1)

        return schema.Marker.search(marker_code, self.session)


    def get_code_pairs(self, sample_ids):
        """ return [(sample_id, code), ...] pairs """
        
        q = self.session.query( schema.Sample.id, schema.Sample.code )
        q = q.filter( schema.Sample.id.in_( sample_ids ))
        return list(q)


    def get_sample_by_id(self, sample_id):
        return schema.Sample.get(sample_id, self.session)


    def get_marker_by_id(self, marker_id):
        return schema.Marker.get(marker_id, self.session)

    def Panel(self):
        p = schema.Panel()
        p._dbh_session_ = self.session
        return p

    
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

        q = self.session.query( schema.AlleleSet.sample_id, schema.Channel.assay_id,
                schema.Allele.marker_id, schema.Allele.bin,
                schema.Allele.size, schema.Allele.height
            ).join(schema.Allele).join(schema.Channel)

        q = q.filter( schema.AlleleSet.sample_id.in_( sample_ids ) )

        if type(params.peaktype) in [ list, tuple ]:
            q = q.filter( schema.Allele.type.in_( params.peaktype  ) )
        else:
            q = q.filter( schema.Allele.type == params.peaktype )

        # we order based on marker_id, sample_id and then descending height
        q = q.order_by( schema.Allele.marker_id, schema.AlleleSet.sample_id,
                        schema.Allele.height.desc() )

        print('MARKER IDS:', marker_ids)

        if marker_ids:
            q = q.filter( schema.AlleleSet.marker_id.in_( marker_ids ) )
            #q = q.outerjoin( Marker, Allele.marker_id == Marker.id )
            #q = q.filter( Marker.id.in_( marker_ids ) )

        if params.abs_threshold > 0:
            q = q.filter( schema.Allele.height > params.abs_threshold )

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


    def get_sample_ids(self, sample_selector):
        pass
