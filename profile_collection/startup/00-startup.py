import pandas as pd

# Make ophyd listen to pyepics.
from ophyd import setup_ophyd
setup_ophyd()

# Set up a Broker.
# TODO clean this up
from databroker import Broker
from databroker.headersource.mongo import MDS
from databroker.assets.mongo import Registry

_mds_config = {'host': 'xf03id-ca1',
               'port': 27017,
               'database': 'datastore-new',
               'timezone': 'US/Eastern'}
mds = MDS(_mds_config, auth=False)

_fs_config = {'host': 'xf03id-ca1',
              'port': 27017,
              'database': 'filestore-new'}
db_new = Broker(mds, Registry(_fs_config))

_mds_config_old = {'host': 'xf03id-ca1',
                   'port': 27017,
                   'database': 'datastore',
                   'timezone': 'US/Eastern'}
mds_old = MDS(_mds_config_old, auth=False)

_fs_config_old = {'host': 'xf03id-ca1',
                  'port': 27017,
                  'database': 'filestore'}
db_old = Broker(mds_old, Registry(_fs_config_old))


# wrapper for two databases
class Broker_New(Broker):

    def __getitem__(self, key):
        try:
            return db_new[key]
        except ValueError:
            return db_old[key]

    def get_table(self, *args, **kwargs):
        result_old = db_old.get_table(*args, **kwargs)
        result_new = db_new.get_table(*args, **kwargs)
        result = [result_old, result_new]
        return pd.concat(result)

    def get_images(self, *args, **kwargs):
        try:
            result = db_new.get_images(*args, **kwargs)
        except IndexError:
            result = db_old.get_images(*args, **kwargs)
        return result


db = Broker_New(mds, Registry(_fs_config))

from hxntools.handlers import register as _hxn_register_handlers
_hxn_register_handlers(db_new)
_hxn_register_handlers(db_old)
del _hxn_register_handlers
# do the rest of the standard configuration
from IPython import get_ipython
from nslsii import configure_base, configure_olog

configure_base(get_ipython().user_ns, db_new, bec=False)
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
