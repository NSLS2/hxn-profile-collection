# vim: sw=4 ts=4 sts expandtab smarttab
# HXN step-scan configuration

import hxntools.scans
from bluesky.global_state import get_gs
import bluesky.plans as bp

gs = get_gs()

hxntools.scans.setup()
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

#gs.DETS = [zebra, sclr1, timepix1, mercury1]#,dexela1]
#gs.DETS = [zebra, sclr1, xspress3, lakeshore2]
gs.DETS = [zebra, sclr1, merlin1, xspress3, lakeshore2]
#gs.DETS = [zebra, sclr1, merlin1, lakeshore2]

#gs.DETS = [zebra, sclr1, merlin1, xspress3, lakeshore2,dexela1]


#gs.DETS = [zebra, sclr1, timepix1]
#gs.TABLE_COLS = ['p_ssx', 'p_ssy', 'p_ssz','sclr1_ch2','sclr1_ch3', 'sclr1_ch4']
gs.TABLE_COLS = ['sclr1_ch2','sclr1_ch3', 'sclr1_ch4', 'sclr1_ch5', 'sclr1_ch5_calc', 'zpssx', 'zpssy', 'zpssz',
                 't_base', 't_sample', 't_vlens', 't_hlens']

# define all the position names and save them to baseline
# need to remove confict names
conflict_name = ['pmllf', 'zplab', 'pmllc']
descs = {d.name: set(d.describe()) for d in bp.separate_devices(ophyd.commands.instances_from_namespace(ophyd.PositionerBase))}

# Plot this by default versus motor position:
gs.PLOT_Y = 'sclr1_ch5_calc' #'dexela_roi1_tot' #'sclr1_ch5_calc' #'Det2_Cr'
gs.OVERPLOT = False
#gs.BASELINE_DEVICES = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, zp]
gs.BASELINE_DEVICES = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, mllosa, zp, zps,smlld]
#gs.BASELINE_DEVICES = [dcm, m1, m2, beamline_status, smll, vmll, hmll, ssa2, bpm1, bpm2, smlld]
