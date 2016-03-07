import asyncio
from functools import partial
from bluesky.standard_config import *
from bluesky.scans import *
from bluesky.callbacks import *
from bluesky.broker_callbacks import *
from bluesky.hardware_checklist import *
from bluesky.qt_kicker import install_qt_kicker

from hxntools.uid_broadcast import HxnUidBroadcast
from hxntools.scan_status import HxnScanStatus
from hxntools.scan_number import HxnScanNumberPrinter


# The following line allows bluesky and pyqt4 GUIs to play nicely together:
install_qt_kicker()


RE = gs.RE
abort = RE.abort
resume = RE.resume
stop = RE.stop

RE.md['group'] = ''
RE.md['config'] = {}
RE.md['beamline_id'] = 'HXN'
# RE.ignore_callback_exceptions = False

uid_broadcaster = HxnUidBroadcast('XF:03IDC-ES{BS-Scan}UID-I')
hxn_scan_status = HxnScanStatus('XF:03IDC-ES{Status}ScanRunning-I')

RE.subscribe('all', HxnScanNumberPrinter())
RE.subscribe('all', uid_broadcaster)
RE.subscribe('all', hxn_scan_status)

loop = asyncio.get_event_loop()
loop.set_debug(False)
# RE.verbose = True

# sr_shutter_status = EpicsSignal('SR-EPS{PLC:1}Sts:MstrSh-Sts', rw=False,
#                                 name='sr_shutter_status')
# sr_beam_current = EpicsSignal('SR:C03-BI{DCCT:1}I:Real-I', rw=False,
#                               name='sr_beam_current')

checklist = partial(basic_checklist,
                    ca_url='http://xf03id-ca.cs.nsls2.local:4800',
                    disk_storage=[('/', 10e9), ('/data', 10e9), ('/xspress3_data', 10e9)],
                    # pv_names=['XF:23ID1-ES{Dif-Ax:SY}Pos-SP'],
                    # pv_conditions=[('XF:23ID-PPS{Sh:FE}Pos-Sts', 'front-end shutter is open', assert_pv_equal, 0),
                    # 		   ('XF:23IDA-PPS:1{PSh}Pos-Sts', 'upstream shutter is open', assert_pv_equal, 0),
                    #                ('XF:23ID1-PPS{PSh}Pos-Sts', 'downstream shutter is open', assert_pv_equal, 0)],
		    )


# # TODO: figure out olog issues
# olog_client = None
# # TODO
RE.logbook = olog_wrapper(olog_client, 'Experiments')
