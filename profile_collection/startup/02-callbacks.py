from functools import partial

from metadatastore.mds import MDS
from databroker import Broker
from databroker.core import register_builtin_handlers
from filestore.fs import FileStore

_mds_config = {'host': 'xf03id-ca1',
               'port': 27017,
               'database': 'datastore',
               'timezone': 'US/Eastern'}
mds = MDS(_mds_config, auth=False)

_fs_config = {'host': 'xf03id-ca1',
               'port': 27017,
               'database': 'filestore'}
db = Broker(mds, FileStore(_fs_config))
register_builtin_handlers(db.fs)


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
RE.subscribe('all', mds.insert)


# Pass on only start/stop documents to a few subscriptions
for _event in ('start', 'stop'):
    RE.subscribe(_event, scan_number_printer)
    RE.subscribe(_event, uid_broadcaster)
    RE.subscribe(_event, hxn_scan_status)

RE.subscribe('start', logger_callback)
