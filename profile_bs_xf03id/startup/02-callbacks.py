from functools import partial

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

LOGBOOKS = ['bs-testing']

logbook = SimpleOlogClient()

uid_signal = EpicsSignal('XF:03IDC-ES{BS-Scan}UID-I')
uid_broadcaster = UidPublish(uid_signal)
hxn_scan_status = HxnScanStatus('XF:03IDC-ES{Status}ScanRunning-I')
cb = logbook_cb_factory(partial(logbook.log, logbooks=LOGBOOKS))

RE = get_gs().RE
RE.subscribe('all', HxnScanNumberPrinter())
RE.subscribe('all', uid_broadcaster)
RE.subscribe('all', hxn_scan_status)
RE.subscribe('start', cb)
