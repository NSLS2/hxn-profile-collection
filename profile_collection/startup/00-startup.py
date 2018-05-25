import pandas as pd

# Make ophyd listen to pyepics.
from ophyd import setup_ophyd
setup_ophyd()

# Set up a Broker.
# TODO clean this up
from databroker import Broker
from databroker.headersource.mongo import MDS
from databroker.assets.mongo import Registry

from datetime import timedelta, datetime, tzinfo

import pymongo
from pymongo import MongoClient

import uuid
from jsonschema import validate as js_validate
import six
from collections import deque

# DB1

db1_name = 'mdb01'
db1_addr = 'xf03id1-mdb01'

# db1_name = 'ca1'
# db1_addr = 'xf03id-ca1'

_mds_config_db1 = {'host': db1_addr,
                   'port': 27017,
                   'database': 'datastore-new',
                   'timezone': 'US/Eastern'}

_fs_config_db1 = {'host': db1_addr,
                  'port': 27017,
                  'database': 'filestore-new'}

# DB2

db2_addr = 'xf03id-ca1'

db2_name = 'ca1-1'
db2_datastore = 'datastore-1'
db2_filestore = 'filestore-1'

#db2_name = 'ca1'
#db2_datastore = 'datastore-new'
#db2_filestore = 'filestore-new'

# db2_name = 'mdb01'
# db2_addr = 'xf03id1-mdb01'

_mds_config_db2 = {'host': db2_addr,
                   'port': 27017,
                   'database': db2_datastore,
                   'timezone': 'US/Eastern'}

_fs_config_db2 = {'host': db2_addr,
                  'port': 27017,
                  'database': db2_filestore}

mongo_client = MongoClient(db2_addr, 27017)

# Benchmark file

f_benchmark = open("benchmark.out", "a+")

# Composite Repository

fs_db2 = mongo_client[db2_filestore]

def sanitize_np(val):
    "Convert any numpy objects into built-in Python types."
    if isinstance(val, (np.generic, np.ndarray)):
        if np.isscalar(val):
            return val.item()
        return val.tolist()
    return val

def apply_to_dict_recursively(d, f):
    for key, val in d.items():
        if hasattr(val, 'items'):
            d[key] = apply_to_dict_recursively(val, f)
        d[key] = f(val)

def _write_to_file(col_name, method_name, t1, t2):
        f_benchmark.write(
            "{0}: {1}, t1: {2} t2:{3} time:{4} \n".format(
                col_name, method_name, t1, t2, (t2-t1),))
        f_benchmark.flush()

class CompositeRegistry(Registry):
    '''Composite registry.'''

    def _register_resource(self, col, uid, spec, root, rpath, rkwargs,
                              path_semantics):

        run_start=None
        ignore_duplicate_error=False
        duplicate_exc=None

        if root is None:
            root = ''

        resource_kwargs = dict(rkwargs)
        if spec in self.known_spec:
            js_validate(resource_kwargs, self.known_spec[spec]['resource'])

        resource_object = dict(spec=str(spec),
                               resource_path=str(rpath),
                               root=str(root),
                               resource_kwargs=resource_kwargs,
                               path_semantics=path_semantics,
                               uid=uid)

        try:
            col.insert_one(resource_object)
        except duplicate_exc:
            if ignore_duplicate_error:
                warnings.warn("Ignoring attempt to insert Datum with duplicate "
                          "datum_id, assuming that both ophyd and bluesky "
                          "attempted to insert this document. Remove the "
                          "Registry (`reg` parameter) from your ophyd "
                          "instance to remove this warning.")
            else:
                raise

        resource_object['id'] = resource_object['uid']
        resource_object.pop('_id', None)
        ret = resource_object['uid']

        return ret


    def register_resource(self, spec, root, rpath, rkwargs,
                              path_semantics='posix'):

        uid = str(uuid.uuid4())

        method_name = "register_resource"

        # db2 database

        col_db2 = fs_db2['resource']

        t1 = datetime.now();
        ret_db2 = self._register_resource(col_db2, uid, spec, root, rpath,
                                            rkwargs, path_semantics='posix')
        t2 = datetime.now()

        _write_to_file(db2_name, method_name, t1, t2);

        # db1 database

        col = self._resource_col

        t1 = datetime.now();
        ret = self._register_resource(col, uid, spec, root, rpath,
                                      rkwargs, path_semantics='posix')
        t2 = datetime.now()

        _write_to_file(db1_name, method_name, t1, t2);

        return ret

    def register_datum(self, resource_uid, datum_kwargs, validate=False):

        if validate:
            raise RuntimeError('validate not implemented yet')

        datum_uid = str(uuid.uuid4())

        # db2 database

        col_db2 = fs_db2['datum']
        datum_db2 = self._api.insert_datum(col_db2, resource_uid, datum_uid, datum_kwargs, {}, None)
        ret_db2 = datum_db2['datum_id']

        # db1 database

        col = self._datum_col
        datum = self._api.insert_datum(col, resource_uid, datum_uid, datum_kwargs, {}, None)
        ret = datum['datum_id']

        return ret

    def _doc_or_uid_to_uid(self, doc_or_uid):

        if not isinstance(doc_or_uid, six.string_types):
            try:
                doc_or_uid = doc_or_uid['uid']
            except TypeError:
                pass

        return doc_or_uid

    def _bulk_insert_datum(self, col, resource, datum_ids,
                           datum_kwarg_list):

        resource_id = self._doc_or_uid_to_uid(resource)

        bulk = col.initialize_unordered_bulk_op()

        d_uids = deque()

        for d_id, d_kwargs in zip(datum_ids, datum_kwarg_list):
            dm = dict(resource=resource_id,
                      datum_id=str(d_id),
                      datum_kwargs=dict(d_kwargs))
            apply_to_dict_recursively(dm, sanitize_np)
            bulk.insert(dm)
            d_uids.append(dm['datum_id'])

        bulk_res = bulk.execute()

        # f_benchmark.write(" _bulk_insert_datum: bulk_res: {0}  \n".format(bulk_res))
        # f_benchmark.flush()

        return d_uids

    def bulk_register_datum_table(self, resource_uid, dkwargs_table,
                                  validate=False):

        if validate:
            raise RuntimeError('validate not implemented yet')

        d_ids = [str(uuid.uuid4()) for j in range(len(dkwargs_table))]

        dkwargs_table = pd.DataFrame(dkwargs_table)
        datum_kwarg_list = [ dict(r) for _, r in dkwargs_table.iterrows()]

        method_name = "bulk_register_datum_table"

        # db2 database

        col_db2 = fs_db2['datum']

        t1 = datetime.now();
        self._bulk_insert_datum(col_db2, resource_uid, d_ids, datum_kwarg_list)
        t2 = datetime.now()

        _write_to_file(db2_name, method_name, t1, t2);

        # db1 database

        t1 = datetime.now();
        self._bulk_insert_datum(self._datum_col, resource_uid, d_ids, datum_kwarg_list)
        t2 = datetime.now()

        _write_to_file(db1_name, method_name, t1, t2);

        ret = d_ids
        return ret

# Broker 1

mds_db1 = MDS(_mds_config_db1, auth=False)
db1 = Broker(mds_db1, CompositeRegistry(_fs_config_db1))

# Broker 2

mds_db2 = MDS(_mds_config_db2, auth=False)
db2 = Broker(mds_db2, CompositeRegistry(_fs_config_db2))


# wrapper for two databases
class Broker_New(Broker):

    def insert(self, name, doc):

        if name == "start":
            f_benchmark.write("\n scan_id: {} \n".format(doc['scan_id']))
            f_benchmark.flush()

        t1 = datetime.now();
        ret1 = db2.insert(name, doc)
        t2 = datetime.now()

        _write_to_file(db2_name, name, t1, t2);

        t3 = datetime.now();
        ret2 = db1.insert(name, doc)
        t4 = datetime.now()

        _write_to_file(db1_name, name, t3, t4);

        return ret2

db = Broker_New(mds_db1, CompositeRegistry(_fs_config_db1))

from hxntools.handlers import register as _hxn_register_handlers
_hxn_register_handlers(db)
#_hxn_register_handlers(db_old)
del _hxn_register_handlers
# do the rest of the standard configuration
from IPython import get_ipython
from nslsii import configure_base, configure_olog

# configure_base(get_ipython().user_ns, db_new, bec=False)
configure_base(get_ipython().user_ns, db, bec=False)
configure_olog(get_ipython().user_ns)

from bluesky.callbacks.best_effort import BestEffortCallback
bec = BestEffortCallback()

# un import *
ns = get_ipython().user_ns
for m in [bp, bps, bpp]:
    for n in dir(m):
        if (not n.startswith('_')
               and n in ns
               and getattr(ns[n], '__module__', '')  == m.__name__):
            del ns[n]
del ns
from bluesky.magics import BlueskyMagics

# set some default meta-data
RE.md['group'] = ''
RE.md['config'] = {}
RE.md['beamline_id'] = 'HXN'
RE.verbose = True

# set up some HXN specific callbacks
from ophyd.callbacks import UidPublish
from hxntools.scan_number import HxnScanNumberPrinter
from hxntools.scan_status import HxnScanStatus
from ophyd import EpicsSignal

uid_signal = EpicsSignal('XF:03IDC-ES{BS-Scan}UID-I', name='uid_signal')
uid_broadcaster = UidPublish(uid_signal)
scan_number_printer = HxnScanNumberPrinter()
hxn_scan_status = HxnScanStatus('XF:03IDC-ES{Status}ScanRunning-I')

# Pass on only start/stop documents to a few subscriptions
for _event in ('start', 'stop'):
    RE.subscribe(scan_number_printer, _event)
    RE.subscribe(uid_broadcaster, _event)
    RE.subscribe(hxn_scan_status, _event)


def ensure_proposal_id(md):
    if 'proposal_id' not in md:
        raise ValueError("You forgot the proposal id.")
# RE.md_validator = ensure_proposal_id


# be nice on segfaults
import faulthandler
faulthandler.enable()

# set up logging framework
import logging
import sys

handler = logging.StreamHandler(sys.stderr)
fmt = logging.Formatter("%(asctime)-15s [%(name)5s:%(levelname)s] %(message)s")
handler.setFormatter(fmt)
handler.setLevel(logging.INFO)

logging.getLogger('hxntools').addHandler(handler)
logging.getLogger('hxnfly').addHandler(handler)
logging.getLogger('ppmac').addHandler(handler)

logging.getLogger('hxnfly').setLevel(logging.DEBUG)
logging.getLogger('hxntools').setLevel(logging.DEBUG)
logging.getLogger('ppmac').setLevel(logging.INFO)


# Flyscan results are shown using pandas. Maximum rows/columns to use when
# printing the table:
pd.options.display.width = 180
pd.options.display.max_rows = None
pd.options.display.max_columns = 10

# enable < shortcut to replace RE(
from bluesky.plan_stubs import  mov
from bluesky.utils import register_transform
register_transform('RE', prefix='<')
