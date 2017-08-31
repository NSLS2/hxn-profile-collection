from functools import partial

from metadatastore.mds import MDS
from databroker import Broker
from databroker.core import register_builtin_handlers
from filestore.fs import FileStore

_mds_config = {'host': 'xf03id-ca1',
               'port': 27017,
               'database': 'datastore-new',
               'timezone': 'US/Eastern'}
mds = MDS(_mds_config, auth=False)

_fs_config = {'host': 'xf03id-ca1',
              'port': 27017,
              'database': 'filestore-new'}
db_new = Broker(mds, FileStore(_fs_config))

_mds_config_old = {'host': 'xf03id-ca1',
               'port': 27017,
               'database': 'datastore',
               'timezone': 'US/Eastern'}
mds_old = MDS(_mds_config_old, auth=False)

_fs_config_old = {'host': 'xf03id-ca1',
              'port': 27017,
              'database': 'filestore'}
db_old = Broker(mds_old, FileStore(_fs_config_old))


from hxntools.handlers.xspress3 import Xspress3HDF5Handler
from hxntools.handlers.timepix import TimepixHDF5Handler

register_builtin_handlers(db_new.fs)

db_new.fs.register_handler(Xspress3HDF5Handler.HANDLER_NAME,
                           Xspress3HDF5Handler)
db_new.fs.register_handler(TimepixHDF5Handler._handler_name,
                           TimepixHDF5Handler, overwrite=True)


register_builtin_handlers(db_old.fs)
db_old.fs.register_handler(Xspress3HDF5Handler.HANDLER_NAME,
                           Xspress3HDF5Handler)
db_old.fs.register_handler(TimepixHDF5Handler._handler_name,
                           TimepixHDF5Handler, overwrite=True)


# wrapper for two databases
class Broker_New(Broker):

    def __getitem__(self, key):
        try:
            return db_new[key]
        except ValueError:
            return db_old[key]

    def get_table(self, *args, **kwargs):
        try:
            return db_new.get_table(*args, **kwargs)
        except:
            return db_old.get_table(*args, **kwargs)

    def get_images(self, *args, **kwargs):
        try:
            return db_new.get_images(*args, **kwargs)
        except:
            return db_old.get_images(*args, **kwargs)



db = Broker_New(mds, FileStore(_fs_config))






def ensure_proposal_id(md):
    if 'proposal_id' not in md:
        raise ValueError("You forgot the proposal id.")

import ophyd
from ophyd import EpicsSignal
from ophyd.callbacks import UidPublish

from hxntools.scan_number import HxnScanNumberPrinter
from hxntools.scan_status import HxnScanStatus

from bluesky.global_state import get_gs
from bluesky.callbacks.olog import logbook_cb_factory

from pyOlog import SimpleOlogClient

# Set up the logbook. This configures bluesky's summaries of
# data acquisition (scan type, ID, etc.).

LOGBOOKS = ['Data Acquisition']

logbook = SimpleOlogClient()

uid_signal = EpicsSignal('XF:03IDC-ES{BS-Scan}UID-I')
uid_broadcaster = UidPublish(uid_signal)
scan_number_printer = HxnScanNumberPrinter()
hxn_scan_status = HxnScanStatus('XF:03IDC-ES{Status}ScanRunning-I')
logger_callback = logbook_cb_factory(partial(logbook.log, logbooks=LOGBOOKS))

RE = get_gs().RE

# Save all scan data to metadatastore:

#from bluesky.register_mds import register_mds
# RE.md_validator = ensure_proposal_id
RE.subscribe('all', mds.insert)


# Pass on only start/stop documents to a few subscriptions
for _event in ('start', 'stop'):
    RE.subscribe(_event, scan_number_printer)
    RE.subscribe(_event, uid_broadcaster)
    RE.subscribe(_event, hxn_scan_status)

RE.subscribe('start', logger_callback)
