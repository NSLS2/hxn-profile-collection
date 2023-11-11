print(f"Loading {__file__!r} ...")

# vim: sw=4 ts=4 sts expandtab smarttab
# HXN step-scan configuration

import hxntools.scans
import bluesky.plans as bp
import bluesky.utils as bu
import ophyd

import functools
from bluesky import plans, Msg
from bluesky import plan_patterns

def _pre_scan(dets, total_points, count_time):
    yield Msg('hxn_next_scan_id')
    yield Msg('hxn_scan_setup', detectors=dets, total_points=total_points,
              count_time=count_time)


@functools.wraps(plans.relative_spiral)
def relative_spiral(dets, x_motor, y_motor, x_range, y_range, dr, nth,
                    time=None, *, per_step=None, md=None, tilt=0.0):

    cyc = plan_patterns.spiral(x_motor, y_motor, x_motor.position,
                                      y_motor.position, x_range, y_range, dr,
                                      nth, tilt=tilt)
    total_points = len(cyc)

    yield from _pre_scan(dets, total_points=total_points, count_time=time)
    return (yield from plans.relative_spiral(
        dets, x_motor, y_motor, x_range,
        y_range, dr, nth, per_step=per_step,
        md=md, tilt=tilt))

hxntools.scans.setup(RE=RE)
ct = bpp.subs_decorator(bec)(hxntools.scans.count)
ascan = bpp.subs_decorator(bec)(hxntools.scans.absolute_scan)
dscan = bpp.subs_decorator(bec)(hxntools.scans.relative_scan)

# fermat = bpp.subs_decorator(bec)(hxntools.scans.relative_fermat)
fermat = bpp.subs_decorator(bec_hxn)(hxntools.scans.relative_fermat)

# spiral = bpp.subs_decorator(bec)(hxntools.scans.relative_spiral)
spiral = bpp.subs_decorator(bec_hxn)(relative_spiral)

mesh = bpp.subs_decorator(bec)(hxntools.scans.absolute_mesh)
dmesh = bpp.subs_decorator(bec)(hxntools.scans.relative_mesh)
d2scan = bpp.subs_decorator(bec)(hxntools.scans.d2scan)
a2scan = bpp.subs_decorator(bec)(hxntools.scans.a2scan)
scan_steps = bpp.subs_decorator(bec)(hxntools.scans.scan_steps)

#Please dont comment out dets 1-4 and dets_fs; it affects the GUI
#dets1 = [zebra, sclr1, merlin1, xspress3]
dets1 = [fs,zebra,sclr1,merlin1, xspress3]
dets_fs = [fs,zebra, sclr1, xspress3]
dets6 = [zebra, sclr1, xspress3]
dets2 = [fs,zebra, sclr1, merlin2, xspress3]
#dets2 = [zebra, sclr1, xspress3, lakeshore2]
dets3 = [zebra, sclr1,xspress3]
#dets3 = [zebra, sclr1,eiger1,xspress3]
dets4 = [zebra, sclr1, merlin1, lakeshore2]
dets7 = [fs, zebra, sclr1, xspress3]
#dets5 = [zebra, sclr1, xspress3, dexela1]
#dets9 = [fs, zebra, sclr1, xspress3, eiger1m_single]
#dets10 = [zebra, sclr1, xspress3, dexela1]
dets8= [fs, zebra, sclr1]
#dets11= [fs, zebra, sclr1, xspress3, eiger1m_single]

# define all the position names and save them to baseline
# need to remove confict names
conflict_name = ['pmllf',
                 'zplab',
                 'pmllc']
#descs = {d.name: set(d.describe())
#         for d in bu.separate_devices
#         (ophyd.utils.instances_from_namespace(ophyd.PositionerBase))}

# sd.baseline = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, zp]
# sd.baseline = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, bpm1, bpm2, smlld]
sd.baseline = [e,
               dcm,
               m1,
               m2,
               smll,
               vmll,
               hmll,
               ssa2,
               s5,
               mllosa,
               zp,
               zps,
               zposa,
               zpbs,
               smlld,
               fdet1,
               diff,
               p,
               ps,
               pp,]
'''
sd.baseline = [ugap,
               e,
               dcm,
               m1,
               m2,
               beamline_status,
               smll,
               vmll,
               hmll,
               ssa2,
               s5,
               mllosa,
               zp,
               zps,
               zposa,
               zpbs,
               smlld,
               fdet1,
               diff,
               p,
               ps,
               pp,]
'''
# The following is a temporary solution, which replaces fixed list of devices passed to 'BlueskyMagics.positioners' (deprecated).
#   TODO: The proper way to label devices is to pass 'labels' to the device (Ophyd object) constructor. Use meaningful labels.
dev_list_motor = [d for  d in bu.separate_devices(ophyd.utils.instances_from_namespace((ophyd.EpicsMotor,)))]
dev_list_pseudo = [d for  d in bu.separate_devices(ophyd.utils.instances_from_namespace((ophyd.PseudoSingle,)))]
for d in dev_list_motor:
    d._ophyd_labels_ = set(["Epics Motor"])
for d in dev_list_pseudo:
    d._ophyd_labels_ = set(["Pseudo Single"])


# BlueskyMagics.positioners = [d for  d in
#                              bu.separate_devices(
#                                  ophyd.utils.instances_from_namespace(
#                                      (ophyd.EpicsMotor, ophyd.PseudoSingle)))]
