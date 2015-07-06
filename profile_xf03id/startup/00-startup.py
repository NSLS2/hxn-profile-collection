import logging

session_mgr._logger.setLevel(logging.INFO)
# session_mgr._logger.setLevel(logging.DEBUG)

handler = logging.StreamHandler(sys.stderr)
fmt = logging.Formatter("%(asctime)-15s [%(name)5s:%(levelname)s] %(message)s")
handler.setFormatter(fmt)
logging.getLogger('hxntools').addHandler(handler)
logging.getLogger('hxnfly').addHandler(handler)

logging.getLogger('hxnfly').setLevel(logging.INFO)
logging.getLogger('hxntools').setLevel(logging.INFO)
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

def run_engine_monkeypatch():
    from ophyd.runengine import RunEngine
    if hasattr(RunEngine, '_demunge_patch'):
        return
    elif not hasattr(RunEngine, '_demunge_names'):
        return
    elif not hasattr(RunEngine, '_demunge_values'):
        return

    RunEngine._demunge_patch = True
    orig_values = RunEngine._demunge_values
    orig_names = RunEngine._demunge_names

    skip_names = ['xspress']
    
    def remove_names(names):
        names = list(names)
        to_remove = set()
        for name in names:
            for skip in skip_names:
                if skip in name:
                    to_remove.add(name)
                    break
        for name in to_remove:
            names.remove(name)
        return names
    
    def names(self, names):
        return orig_names(self, remove_names(names))
        
    def values(self, vals, keys):
        return orig_values(self, vals, remove_names(keys))

    RunEngine._demunge_values = values
    RunEngine._demunge_names = names


run_engine_monkeypatch()
