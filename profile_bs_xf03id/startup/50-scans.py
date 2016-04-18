# vim: sw=4 ts=4 sts expandtab smarttab
# HXN step-scan configuration

import hxntools.scans
from hxntools.spiral_scans import HxnFermatScan

from bluesky.global_state import get_gs

gs = get_gs()


# Set up regular ascans/dscans to work with HXN detector triggering methods:
hxntools.scans.setup()

# Define how spiral scans should work
fermat = HxnFermatScan()


# default_detectors = [sclr1, det_beamstatus, sclr1_ch2, sclr1_ch3, sclr1_ch4,
#                      ssx_rbv, ssy_rbv, ssz_rbv, t_base, t_sample, t_vlens,
#                      t_hlens]

# gs.DETS = [zebra, master_sclr1, sclr1_ch2, sclr1_ch3, sclr1_ch4, ssx, ssy,
#   ssz, t_base, t_sample, t_vlens, t_hlens, xbpm_x, xbpm_y, quad_x, quad_y,
#   dcm_th, dcm_p, angle_x, angle_y, slit1_top, slit1_bottom, slit1_right,
#   slit1_left, slit1_xpos, slit1_ypos]

gs.TABLE_COLS = ['sclr2_ch2','sclr2_ch3', 'sclr2_ch4', 'ssx', 'ssy', 'ssz',
                 't_base', 't_sample', 't_vlens', 't_hlens']


# Plot this by default versus motor position:
gs.PLOT_Y = 'Det2_V'
