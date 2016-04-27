from ophyd.commands import setup_ophyd
setup_ophyd()

from ophyd.commands import *
from bluesky.callbacks import *
from bluesky.plans import *
# from bluesky.spec_api import *
from bluesky.utils import install_qt_kicker
from bluesky.global_state import (get_gs, abort, stop, resume, panic,
                                  all_is_well)

from databroker import (DataBroker as db, get_events, get_images, get_table,
                        get_fields, restream, process)


# The following line allows bluesky and pyqt4 GUIs to play nicely together:
install_qt_kicker()

gs = get_gs()
RE = gs.RE
RE.md['group'] = ''
RE.md['config'] = {}
RE.md['beamline_id'] = 'HXN'

RE.verbose = True
# RE.ignore_callback_exceptions = False
