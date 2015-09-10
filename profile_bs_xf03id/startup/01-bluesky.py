import asyncio
from functools import partial
from bluesky.standard_config import *
from bluesky.scans import *
from bluesky.callbacks import *
from bluesky.broker_callbacks import *
from bluesky.hardware_checklist import *


gs.RE.md['group'] = ''
gs.RE.md['config'] = {}
gs.RE.md['beamline_id'] = 'HXN'

class HxnScanNumberPrinter:
    def __init__(self):
        self._last_start = None

    def __call__(self, name, doc):
        if name == 'start':
            self._last_start = doc
        if self._last_start is None:
            return
        if name in ('start', 'stop'):
            print('Scan ID: {scan_id} [{uid}]'.format(**self._last_start))


gs.RE.subscribe('all', HxnScanNumberPrinter())

loop = asyncio.get_event_loop()
loop.set_debug(False)
# gs.RE.verbose = True

# sr_shutter_status = EpicsSignal('SR-EPS{PLC:1}Sts:MstrSh-Sts', rw=False,
#                                 name='sr_shutter_status')
# sr_beam_current = EpicsSignal('SR:C03-BI{DCCT:1}I:Real-I', rw=False,
#                               name='sr_beam_current')

checklist = partial(basic_checklist,
                    ca_url='http://xf03id-ca.cs.nsls2.local:4800',
                    # disk_storage=[('/', 1e9)],
                    # pv_names=['XF:23ID1-ES{Dif-Ax:SY}Pos-SP'],
                    # pv_conditions=[('XF:23ID-PPS{Sh:FE}Pos-Sts', 'front-end shutter is open', assert_pv_equal, 0),
                    # 		   ('XF:23IDA-PPS:1{PSh}Pos-Sts', 'upstream shutter is open', assert_pv_equal, 0),
                    #                ('XF:23ID1-PPS{PSh}Pos-Sts', 'downstream shutter is open', assert_pv_equal, 0)],
		    )
