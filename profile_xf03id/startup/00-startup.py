import logging

session_mgr._logger.setLevel(logging.INFO)
# session_mgr._logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stderr)
fmt = logging.Formatter("%(asctime)-15s [%(name)5s:%(levelname)s] %(message)s")
handler.setFormatter(fmt)
logging.getLogger('hxntools').addHandler(handler)
logging.getLogger('hxnfly').addHandler(handler)
logging.getLogger('hxnfly').setLevel(logging.INFO)
logging.getLogger('ppmac').setLevel(logging.INFO)

try:
    from ophyd.commands import *
except ImportError:
    from ophyd.userapi import *

# NOTE: may need to have enaml.qt.prepare_pyqt() in IPython configuration before
# anything else is loaded...
from dataportal import (DataBroker as db,
                        StepScan as ss,
                        DataBroker, StepScan,
                        DataMuxer)
