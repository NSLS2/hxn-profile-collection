from bluesky.suspenders import (SuspendFloor, SuspendBoolHigh, SuspendBoolLow)
from bluesky.global_state import get_gs

gs = get_gs()
RE = gs.RE


# Here are some conditions that will cause scans to pause automatically:
# - when the beam current goes below a certain threshold
susp_current = SuspendFloor(beamline_status.beam_current,
                            suspend_thresh=100.0,
                            resume_thresh=105.0,
                            # message='beam current too low',
                            )
# RE.install_suspender(susp_current)

# - when the shutter is closed
susp_shutter = SuspendBoolLow(beamline_status.shutter_status,
                              # message='shutter not open',
                              )
# RE.install_suspender(susp_shutter)

# - if the beamline isn't enabled
susp_enabled = SuspendBoolLow(beamline_status.beamline_enabled,
                              # message='beamline is not enabled',
                              )
# RE.install_suspender(susp_enabled)
