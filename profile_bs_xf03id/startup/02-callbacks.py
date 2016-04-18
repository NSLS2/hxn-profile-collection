import ophyd
from ophyd import EpicsSignal
from ophyd.callbacks import UidPublish

from hxntools.scan_number import HxnScanNumberPrinter
from hxntools.scan_status import HxnScanStatus

from bluesky.global_state import get_gs


uid_signal = EpicsSignal('XF:03IDC-ES{BS-Scan}UID-I')
uid_broadcaster = UidPublish(uid_signal)
hxn_scan_status = HxnScanStatus('XF:03IDC-ES{Status}ScanRunning-I')

gs = get_gs()
gs.RE.subscribe('all', HxnScanNumberPrinter())
gs.RE.subscribe('all', uid_broadcaster)
gs.RE.subscribe('all', hxn_scan_status)
