# vim: sw=4 ts=4 sts expandtab smarttab
# HXN step-scan configuration

import hxntools.scans
import bluesky.plans as bp
import ophyd

hxntools.scans.setup(RE)
ct = hxntools.scans.count
ascan = hxntools.scans.absolute_scan
dscan = hxntools.scans.relative_scan
fermat = hxntools.scans.relative_fermat
spiral = hxntools.scans.relative_spiral
mesh = hxntools.scans.absolute_mesh
dmesh = hxntools.scans.relative_mesh
d2scan = hxntools.scans.d2scan
a2scan = hxntools.scans.a2scan
scan_steps = hxntools.scans.scan_steps

dets1 = [zebra, sclr1, timepix1]
dets2 = [zebra, sclr1, xspress3, lakeshore2]
dets3 = [zebra, sclr1, merlin1, xspress3, lakeshore2,quad]
dets4 = [zebra, sclr1, merlin1, lakeshore2]
dets5 = [zebra, sclr1, merlin1, xspress3, lakeshore2,dexela1]


# define all the position names and save them to baseline
# need to remove confict names
conflict_name = ['pmllf', 'zplab', 'pmllc']
descs = {d.name: set(d.describe())
         for d in bp.separate_devices
         (ophyd.utils.instances_from_namespace(ophyd.PositionerBase))}

# sd.baseline = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, zp]
# sd.baseline = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, bpm1, bpm2, smlld]
sd.baseline = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, mllosa, zp, zps, smlld, fdet1]
