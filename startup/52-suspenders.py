print(f"Loading {__file__!r} ...")

from bluesky.suspenders import (SuspendFloor, SuspendBoolHigh, SuspendBoolLow)


# Here are some conditions that will cause scans to pause automatically:
# - when the beam current goes below a certain threshold
susp_current = SuspendFloor(beamline_status.beam_current,
                            suspend_thresh=150.0,
                            resume_thresh=240.0,
                            tripped_message='beam current too low',
                            )

# - when the shutter is closed
susp_shutter = SuspendBoolLow(beamline_status.shutter_status,
                              tripped_message='shutter not open',
                              )

# - if the beamline isn't enabled
susp_enabled = SuspendBoolLow(beamline_status.beamline_enabled,
                              tripped_message='beamline is not enabled',
                              )

susp_cryo_cooler = SuspendBoolHigh(beamline_status.cryo_filling,
                                   sleep=60*10,
                                   tripped_message='Cyrocooler is refilling')

# NOTE: to enable or disable the suspenders, (un)comment the following:
# Install all suspenders:
#RE.install_suspender(susp_current)
#RE.install_suspender(susp_shutter)
#RE.install_suspender(susp_enabled)
#RE.install_suspender(susp_cryo_cooler)
