import logging

session_mgr._logger.setLevel(logging.INFO)
# To help with debugging scanning-related problems, uncomment the following:
# session_mgr._logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stderr)
fmt = logging.Formatter("%(asctime)-15s [%(name)5s:%(levelname)s] %(message)s")
handler.setFormatter(fmt)
logging.getLogger('hxntools').addHandler(handler)
logging.getLogger('hxnfly').addHandler(handler)
logging.getLogger('ppmac').addHandler(handler)

logging.getLogger('hxnfly').setLevel(logging.DEBUG)
logging.getLogger('hxntools').setLevel(logging.DEBUG)
logging.getLogger('ppmac').setLevel(logging.INFO)

from ophyd.commands import *
import databroker
from databroker import DataBroker as db
