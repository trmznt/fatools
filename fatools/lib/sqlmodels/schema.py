
from sqlalchemy.ext.declarative import declarative_base
from sqlalchemy.orm import scoped_session, sessionmaker, mapper
from sqlalchemy.engine import Engine
from sqlalchemy import event, create_engine, func

from sqlalchemy import and_, or_, schema, types, MetaData, Sequence, Column, ForeignKey, UniqueConstraint, Table
from sqlalchemy.orm import relationship, backref, dynamic_loader, deferred, reconstructor
from sqlalchemy.orm.collections import column_mapped_collection, attribute_mapped_collection
from sqlalchemy.orm.interfaces import MapperExtension
from sqlalchemy.orm.exc import NoResultFound
from sqlalchemy.orm.session import object_session
from sqlalchemy.exc import OperationalError, IntegrityError
from sqlalchemy.ext.associationproxy import association_proxy
from sqlalchemy.ext.declarative import declared_attr, declarative_base
from sqlalchemy.sql.functions import current_timestamp

from zope.sqlalchemy import ZopeTransactionExtension

from fatools.lib.fautil.mixin import ( PanelMixIn, AssayMixIn, ChannelMixIn, MarkerMixIn,
                BinMixIn, AlleleSetMixIn, AlleleMixIn, SampleMixIn, BatchMixIn,
                NoteMixIn, BatchNoteMixIn, SampleNoteMixIn, AssayNoteMixIn,
                ChannelNoteMixIn, AlleleSetNoteMixIn, PanelNoteMixIn, MarkerNoteMixIn )

import os, io, yaml

#__all__ = ['get_base', 'get_dbsession', 'set_datalogger']

# this is necessary for SQLite to use FOREIGN KEY support (as well as ON DELETE CASCADE)
@event.listens_for(Engine, 'connect')
def set_sqlite_pragma(dbapi_connection, connection_record):
    cursor = dbapi_connection.cursor()
    cursor.execute('PRAGMA foreign_keys=ON')
    cursor.close()


Base = declarative_base()

def _generic_query(cls, session):
    return session.query(cls)

Base.query = classmethod(_generic_query)

def _generic_get(cls, pkid, session):
    q = cls.query(session)
    return q.get(int(pkid))

Base.get = classmethod(_generic_get)

def _generic_lowername(cls):
    return cls.__name__.lower()

Base.lowername = classmethod(_generic_lowername)

def _generic_delete(cls, pkid, session):
    q = cls.query(session)
    return q.filter(cls.id == int(pkid)).delete()

Base.delete = classmethod(_generic_delete)


# create JSON column
# XXX: may be more appropriate to create a dict-based object that will serialize to JSON?

null = object()

class JSONCol(types.TypeDecorator):
    impl = types.Unicode

    def process_bind_param(self, value, dialect):
        if value is null:
            value = None
        return json.dumps(value)

    def process_result_value(self, value, dialect):
        if value is None:
            return None
        return json.loads(value)

    def copy_value(self, value):
        return copy.deepcopy(value)

# create YAML column
#

null = null

class YAMLCol(types.TypeDecorator):
    impl = types.Unicode

    def process_bind_param(self, value, dialect):
        if value is null:
            value = None
        return yaml.dump(value, default_flow_style=True)

    def process_result_value(self, value, dialect):
        if value is None:
            return None
        return yaml.load(value)

    def copy_value(self, value):
        return copy.deepcopy(value)

import numpy, copy

class NPArray(types.TypeDecorator):
    impl = types.LargeBinary

    def process_bind_param(self, value, dialect):
        if value is None:
            return None
        #buf = value.tostring()
        buf = io.BytesIO()
        numpy.save(buf, value)
        return buf.getvalue()

    def process_result_value(self, value, dialect):
        if value is None:
            return None
        buf = io.BytesIO(value)
        return numpy.load(buf)
        #return numpy.fromstring(value)

    def copy_value(self, value):
        return copy.deepcopy( value )



class Note(Base, NoteMixIn):

    __tablename__ = 'notes'
    id = Column(types.Integer, primary_key=True)
    text = Column(types.String(1024), nullable=False, default='')
    cat = Column(types.String(32), nullable=False, default='')
    stamp = Column(types.DateTime, nullable=False)



class Batch(Base, BatchMixIn):

    __tablename__ = 'batches'
    id = Column(types.Integer, primary_key=True)
    code = Column(types.String(16), nullable=False, unique=True)
    assay_provider = Column(types.String(32), nullable=False)
    species = Column(types.String(16), nullable=False)
    description = Column(types.String(1024), nullable=False, default='')
    remark = deferred(Column(types.String(1024), nullable=True))
    data = deferred(Column(YAMLCol(4096), nullable=False, default=''))
    bin_batch_id = Column(types.Integer, ForeignKey('batches.id'), nullable=True)




    def add_sample(self, sample_code):
        """ return a new Sample with sample_code """

        sample = Sample()
        sample.code = sample_code
        sample.batch = self

        return sample


    def search_sample(self, sample_code):
        """ return a single Sample from the current batch with sample_code """
        try:
            return self.samples.filter( func.lower(Sample.code) == func.lower(sample_code),
                Sample.batch_id == self.id ).one()
        except NoResultFound:
            return None


    @staticmethod
    def search(code, session):
        """ provide case-insensitive search for batch code """
        q = Batch.query(session).filter( func.lower(Batch.code) == func.lower(code) )
        return q.one()


    def get_panel(self, panel_code):
        return Panel.search(panel_code, object_session(self))


    def get_marker(self, marker_code, species=None):
        session = object_session(self)
        markers = None  #XXX: Fix me
        raise NotImplementedError()


    @property
    def sample_ids(self):
        """ faster implementation of getting sample ids """
        session = object_session(self)
        return [ x[0] for x in session.query( Sample.id ).filter( Sample.batch_id == self.id ) ]



class BatchNote(Base, BatchNoteMixIn):

    __tablename__ = 'batchnotes'
    id = Column(types.Integer, primary_key=True)
    batch_id = Column(types.Integer, ForeignKey('batches.id', ondelete='CASCADE'),
                nullable=False)
    note_id = Column(types.Integer, ForeignKey('notes.id', ondelete='CASCADE'),
                nullable=False)



class Sample(Base, SampleMixIn):

    __tablename__ = 'samples'
    id = Column(types.Integer, primary_key=True)
    code = Column(types.String(64), nullable=False)
    type = Column(types.String(1), default='S')
    altcode = Column(types.String(16), nullable=True)               # custom usage
    category = Column(types.Integer, nullable=False, default=0)     # custom usage
    batch_id = Column(types.Integer, ForeignKey('batches.id', ondelete='CASCADE'),
                nullable=False)
    int1 = Column(types.Integer, nullable=False, default=-1)        # custom usage
    int2 = Column(types.Integer, nullable=False, default=-1)        # custom usage
    string1 = Column(types.String(16), nullable=False, default='')  # custom usage
    string2 = Column(types.String(16), nullable=False, default='')  # custom usage
    batch = relationship(Batch, uselist=False, backref=backref('samples', lazy='dynamic'))
    remark = deferred(Column(types.String(1024), nullable=True))

    __table_args__ = (  UniqueConstraint( 'code', 'batch_id' ),
                        UniqueConstraint( 'altcode', 'batch_id')
                    )

    def new_assay(self, raw_data, filename, status, panel=None):
        assay = Assay( raw_data = raw_data, filename = filename )
        if panel is None:
            panel = Panel.search('undefined', object_session(self))
        assay.panel = panel
        assay.sample = self
        assay.status = status
        return assay



class SampleNote(Base, SampleNoteMixIn):

    __tablename__ = 'samplenotes'
    id = Column(types.Integer, primary_key=True)
    sample_id = Column(types.Integer, ForeignKey('samples.id', ondelete='CASCADE'),
                nullable=False)
    note_id = Column(types.Integer, ForeignKey('notes.id', ondelete='CASCADE'),
                nullable=False)



class Panel(Base, PanelMixIn):

    __tablename__ = 'panels'
    id = Column(types.Integer, primary_key=True)
    code = Column(types.String(8), nullable=False, unique=True)
    data = Column(YAMLCol(1024), nullable=False)
    remark = deferred(Column(types.String(1024), nullable=True))

    def update(self, obj):

        self._update(obj)

        # verify that each marker in data exists
        session = object_session(self) or self._dbh_session_
        for m_code in self.data['markers']:
            m = Marker.search(m_code, session)
            if m is None:
                cerr("ERR: can't find marker: %s" % m_code)
                sys.exit(1)


    def sync(self, session):
        """ sync'ing is getting the object in the database and update with our content """
        db_panel = Panel.search(self.code, session)
        db_panel.update( self )
        return db_panel


    def get_marker(self, marker_code):
        """ return marker instance """
        return Marker.search( marker_code, object_session(self))

    @reconstructor
    def init_data(self):
        if not hasattr(self, '_dyes'):
            self._dyes = {}


    @staticmethod
    def search(code, session):
        """ provide case-insensitive search for marker code """
        q = Panel.query(session).filter( func.lower(Panel.code) == func.lower(code) )
        return q.one()



class PanelNote(Base, PanelNoteMixIn):

    __tablename__ = 'panelnotes'
    id = Column(types.Integer, primary_key=True)
    panel_id = Column(types.Integer, ForeignKey('panels.id', ondelete='CASCADE'),
                nullable=False)
    note_id = Column(types.Integer, ForeignKey('notes.id', ondelete='CASCADE'),
                nullable=False)



class Marker(Base, MarkerMixIn):

    __tablename__ = 'markers'
    id = Column(types.Integer, primary_key=True)
    code = Column(types.String(64), nullable=False, unique=True)
    species = Column(types.String(32), nullable=False)
    repeats = Column(types.Integer, nullable=False, default=-1)

    min_size = Column(types.Integer, nullable=False, default=0)
    max_size = Column(types.Integer, nullable=False, default=0)
    """ range of allele size for this marker """

    bins = deferred(Column(YAMLCol(2048), nullable=False, default=''))
    """ sorted known bins for this markers """
    
    related_to_id = Column(types.Integer, ForeignKey('markers.id'),
                          nullable=True)
    related_to = relationship("Marker", uselist=False)
    """ points to related marker """
    
    z_params = Column(types.String(32), nullable=False, default='')
    """ mathematical expression correlating with the related_to marker """

    remark = deferred(Column(types.String(1024), nullable=True))

    __table_args__ = ( UniqueConstraint( 'code', 'species' ), )

    def update(self, obj):
        
        self._update( obj )
        if type(obj) == dict and 'related_to' in obj:
            related_marker = Marker.search( d['related_to'],
                    session = object_session(self) or self.__dbh_session )
            self.related_to = related_marker


    def sync(self, session):
        """ sync assume that the current instance is not attached to any session """
        db_marker = Marker.search(marker.code, session=session)
        db_marker.update( self )
        return db_marker


    @staticmethod
    def search(code, session):
        """ provide case-insensitive search for marker code """
        if '/' in code:
            species, code = code.split('/')
            q = Marker.query(session).filter( func.lower(Marker.code) == func.lower(code),
                                    func.lower(Marker.species) == func.lower(species) )
        else:
            q = Marker.query(session).filter( func.lower(Marker.code) == func.lower(code) )
        return q.one()



class MarkerNote(Base, MarkerNoteMixIn):

    __tablename__ = 'markernotes'
    id = Column(types.Integer, primary_key=True)
    marker_id = Column(types.Integer, ForeignKey('markers.id', ondelete='CASCADE'),
                nullable=False)
    note_id = Column(types.Integer, ForeignKey('notes.id', ondelete='CASCADE'),
                nullable=False)

class Bin(Base, BinMixIn):

    __tablename__ = 'bins'
    id = Column(types.Integer, primary_key=True)
    batch_id = Column(types.Integer, ForeignKey('batches.id'), nullable=False)
    marker_id = Column(types.Integer, ForeignKey('markers.id'), nullable=False)
    z = deferred(Column(NPArray))

    related_to_id = Column(types.Integer, ForeignKey('bins.id'), nullable=True)

    bins = deferred(Column(YAMLCol(2048), nullable=False, default=''))
    """ sorted known bins for this markers """

    meta = deferred(Column(YAMLCol(4096), nullable=False, default=''))
    """ metadata for this bin """

    remark = deferred(Column(types.String(512)))


    def search(self, batch_id, marker_id, session):
        q = Bin.query(session).filter(Bin.batch_id == batch_id, Bin.marker_id == marker_id)
        return q.one()


class Assay(Base, AssayMixIn):

    __tablename__ = 'assays'
    id = Column(types.Integer, primary_key=True)
    filename = Column(types.String(128), nullable=False, index=True)
    runtime = Column(types.DateTime, nullable=False)
    rss = Column(types.Float, nullable=False, default=-1)
    dp = Column(types.Float, nullable=False, default=-1)
    score = Column(types.Float, nullable=False, default=-1)
    z = deferred(Column(NPArray))
    ladder_peaks = Column(types.Integer, nullable=False, default=-1)

    size_standard = Column(types.String(32), nullable=False, default='')

    sample_id = Column(types.Integer, ForeignKey('samples.id', ondelete='CASCADE'),
                        nullable=False)
    sample = relationship(Sample, uselist=False, backref=backref('assays', lazy='dynamic'))

    panel_id = Column(types.Integer, ForeignKey('panels.id'), nullable=False)
    panel = relationship(Panel, uselist=False)


    ladder_id = Column(types.Integer,
                        ForeignKey('channels.id', use_alter=True, name = 'ladderchannel_fk'),
                        nullable=True)

    #channels = relationship('Channel', primaryjoin = "Assay.id == Channel.id", lazy='dynamic',
    #                    post_update = True,
    #                    backref = backref('assay', uselist=False))

    ladder = relationship('Channel', uselist=False,
                primaryjoin = "Assay.ladder_id == Channel.id")

    status = Column(types.String(32), nullable=False)
    method = deferred(Column(types.String(16), nullable=False))
    report = deferred(Column(types.String(512), nullable=False, default=''))
    remark = deferred(Column(types.String(1024), nullable=False, default=''))

    raw_data = deferred(Column(types.Binary(), nullable=False))
    """ raw data for this assay (FSA file content) """

    __table_args__ = (  UniqueConstraint( 'filename', 'panel_id', 'sample_id' ), )


    def new_channel(self, raw_data, data, dye, wavelen, status, median, mean,
            max_height, min_height, std_dev, initial_marker=None, initial_panel=None):
        """ create new channel and added to this assay """
        if not initial_marker:
            initial_marker = Marker.search('undefined', session = object_session(self))
        if not initial_panel:
            initial_panel = Panel.search('undefined', session = object_session(self))
        channel = Channel( raw_data = data, data = data, dye = dye, wavelen = wavelen,
                            status = status, median = median, mean = mean,
                            max_height = max_height, min_height = min_height,
                            std_dev = std_dev )
        channel.assay = self
        channel.marker = initial_marker
        channel.panel = initial_panel

        return channel


    def get_ladder(self):
        """ get ladder channel """
        assert self.ladder_id, "ERR/PROG -  pls make sure ladder_id is not null!"
        session = object_session(self)
        return Channel.get(self.ladder_id, session)



class AssayNote(Base, AssayNoteMixIn):

    __tablename__ = 'assaynotes'
    id = Column(types.Integer, primary_key=True)
    assay_id = Column(types.Integer, ForeignKey('assays.id', ondelete='CASCADE'),
                nullable=False)
    note_id = Column(types.Integer, ForeignKey('notes.id', ondelete='CASCADE'),
                nullable=False)



class Channel(Base, ChannelMixIn):

    __tablename__ = 'channels'
    id = Column(types.Integer, primary_key=True)

    assay_id = Column(types.Integer, ForeignKey('assays.id', ondelete='CASCADE'),
                        nullable=False)
    assay = relationship(Assay, uselist=False, primaryjoin = assay_id == Assay.id, 
                    backref=backref('channels', lazy='dynamic'))

    marker_id = Column(types.Integer, ForeignKey('markers.id'), nullable=False)
    marker = relationship(Marker, uselist=False, backref=backref('channels', lazy='dynamic'))

    dye = Column(types.String(32), nullable=False)

    markers = relationship(Marker, secondary='channels_markers', viewonly=True)

    raw_data = deferred(Column(NPArray, nullable=False))
    """ raw data from channel as numpy array, can have empty array to accomodate
        allele data from CSV uploading """

    status = Column(types.String(32), nullable=False)

    wavelen = Column(types.Integer, nullable=False, default=0)
    median = Column(types.Integer, nullable=False, default=0)
    mean = Column(types.Float, nullable=False, default=0.0)
    std_dev = Column(types.Float, nullable=False, default=0.0)
    max_height = Column(types.Integer, nullable=False, default=-1)
    min_height = Column(types.Integer, nullable=False, default=-1)
    """ basic descriptive statistics for data"""
        
    data = deferred(Column(NPArray, nullable=False))
    """ data after smoothed using savitzky-golay algorithm and baseline correction
        using top hat morphologic transform
    """

    remark = deferred(Column(types.String(1024), nullable=True))


    def new_alleleset(self, revision=-1):
        return AlleleSet( channel = self, sample = self.assay.sample,
                            marker = self.marker )


    def clear(self):
        for alleleset in self.allelesets:
            del alleleset


    def get_latest_alleleset(self):
        if self.allelesets.count() < 1:
            raise RuntimeError("ERR - channel does not have alleleset, probably hasn't been scanned!")
        return self.allelesets[-1]



class ChannelNote(Base, ChannelNoteMixIn):

    __tablename__ = 'channelnotes'
    id = Column(types.Integer, primary_key=True)
    channel_id = Column(types.Integer, ForeignKey('channels.id', ondelete='CASCADE'),
                nullable=False)
    note_id = Column(types.Integer, ForeignKey('notes.id', ondelete='CASCADE'),
                nullable=False)



channels_markers = Table(
    'channels_markers', Base.metadata,
    Column('channel_id', types.Integer, ForeignKey('channels.id')),
    Column('marker_id', types.Integer, ForeignKey('markers.id')),
    UniqueConstraint( 'channel_id', 'marker_id' )
    )
    # unique compound (channel_id, marker_id)



class AlleleSet(Base, AlleleSetMixIn):

    __tablename__ = 'allelesets'
    id = Column(types.Integer, primary_key=True)

    channel_id = Column(types.Integer, ForeignKey('channels.id', ondelete='CASCADE'),
                    nullable=False)
    channel = relationship(Channel, uselist=False,
                backref=backref('allelesets', lazy='dynamic', passive_deletes=True))
    # a channel can have several allele set for different revision numbers

    sample_id = Column(types.Integer, ForeignKey('samples.id', ondelete='CASCADE'),
                    nullable=False)
    sample = relationship(Sample, uselist=False,
                backref=backref('allelesets', lazy='dynamic', passive_deletes=True))
    """ link to sample """

    marker_id = Column(types.Integer, ForeignKey('markers.id'), nullable=False)
    marker = relationship(Marker, uselist=False,
                    backref=backref('allelesets', lazy='dynamic'))
    """ link to marker """

    scanning_method = deferred(Column(types.String(32), nullable=False))
    """ method used for scanning and generating this alleleset """

    calling_method = deferred(Column(types.String(32), nullable=False))
    """ method used for calling this alleleset """

    binning_method = deferred(Column(types.String(32), nullable=False))
    """ method used for binning this alleleset """


    def new_allele(self, rtime, height, area, brtime, ertime, wrtime, srtime, beta, theta,
                    type, method):
        allele = Allele( rtime = rtime, height = height, area = area,
                    brtime = brtime, ertime = ertime, wrtime = wrtime, srtime = srtime,
                    beta = beta, theta = theta, type = type, method = method )
        allele.alleleset = self

        return allele



class AlleleSetNote(Base, AlleleSetNoteMixIn):

    __tablename__ = 'allelesetnotes'
    id = Column(types.Integer, primary_key=True)
    alleleset_id = Column(types.Integer, ForeignKey('allelesets.id', ondelete='CASCADE'),
                nullable=False)
    note_id = Column(types.Integer, ForeignKey('notes.id', ondelete='CASCADE'),
                nullable=False)



class Allele(Base, AlleleMixIn):

    __tablename__ = 'alleles'
    id = Column(types.Integer, primary_key=True)

    alleleset_id = Column(types.Integer, ForeignKey('allelesets.id', ondelete='CASCADE'),
                nullable=False)
    alleleset = relationship(AlleleSet, uselist=False,
                backref=backref('alleles', cascade='all, delete-orphan',
                passive_deletes=True))

    marker_id = Column(types.Integer, ForeignKey('markers.id', ondelete='CASCADE'),
                nullable=False)
    marker = relationship(Marker, uselist=False,
                backref=backref('alleles', cascade='all, delete-orphan',
                    passive_deletes=True))

    abin = Column(types.Integer, nullable=False, default=-1)    # adjusted bin
    asize = Column(types.Float, nullable=False, default=-1)     # adjusted size
    aheight = Column(types.Float, nullable=False, default=-1)   # adjusted height

    bin = Column(types.Integer, nullable=False, default=-1)
    size = Column(types.Float, nullable=False, default=-1)
    deviation = Column(types.Float, nullable=False, default=-1)
    # deviation -> for ladder channel, this is ( z(rtime) - size )**2 or square of residual
    # for marker channel, this depends on the method
    # method cubic-spline, this is avg of deviation of the nearest peaks
    # for local southern, this is (size1 - size2) ** 2

    height = Column(types.Float, nullable=False, default=-1)
    area = Column(types.Float, nullable=False, default=-1)
    rtime = Column(types.Integer, nullable=False, default=-1)
    delta = Column(types.Float, nullable=False, default=-1)     # bin - actual size 
    beta = Column(types.Float, nullable=False, default=-1)      # area / height
    theta = Column(types.Float, nullable=False, default=-1)     # height / width

    type = Column(types.String(32), nullable=False)
    method = Column(types.String(32), nullable=False)   # binning method

    brtime = Column(types.Integer, nullable=False, default=-1)
    ertime = Column(types.Integer, nullable=False, default=-1)
    wrtime = Column(types.Integer, nullable=False, default=-1)
    srtime = Column(types.Float, nullable=False, default=-1)    # log2( right_area/left_area )
    """ begin, end, width, symmetrical retention time of this peak and peak quality"""

    qscore = Column(types.Float, nullable=False, default=-1)    # calculated in preannotate()
    qcall = Column(types.Float, nullable=False, default=-1)     # calculated in call()


    @property
    def channel(self):
        return self.alleleset.channel



def engine_from_file( dbfilename, bind=True ):
    """ return (engine, session) """

    # make absolute path
    if dbfilename != ':memory:':
        abspath = os.path.abspath( dbfilename )
        engine = create_engine('sqlite:///' + abspath)
    else:
        # use memory-based sqlite database
        engine = create_engine('sqlite://')
    session = scoped_session(sessionmaker(extension=ZopeTransactionExtension()))
    if bind:
        session.configure(bind = engine)
    return (engine, session)



