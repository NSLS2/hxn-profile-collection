# vim: sw=4 ts=4 sts expandtab smarttab
# HXN step-scan configuration

import hxntools.scans
from bluesky.global_state import get_gs

gs = get_gs()

hxntools.scans.setup()
ct = hxntools.scans.count
ascan = hxntools.scans.absolute_scan
dscan = hxntools.scans.relative_scan
fermat = hxntools.scans.relative_fermat
spiral = hxntools.scans.relative_spiral
mesh = hxntools.scans.absolute_mesh
dmesh = hxntools.scans.relative_mesh

gs.DETS = [zebra, sclr1, merlin1, xspress3, smll, lakeshore2, xbpm, dcm, s1]
gs.TABLE_COLS = ['sclr2_ch2','sclr2_ch3', 'sclr2_ch4', 'ssx', 'ssy', 'ssz',
                 't_base', 't_sample', 't_vlens', 't_hlens']


# Plot this by default versus motor position:
gs.PLOT_Y = 'Det2_V'
gs.BASELINE_DEVICES = [smll]
