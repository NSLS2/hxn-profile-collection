# vim: sw=4 ts=4 sts expandtab smarttab
# HXN step-scan configuration

import hxntools.scans
from hxntools.detectors.master_detector import MasterDetector
from hxntools.spiral_scans import HxnFermatScan


# Set up regular ascans/dscans to work with HXN detector triggering methods:
hxntools.scans.setup()

# Define how spiral scans should work
fermat = HxnFermatScan()


# default_detectors = [sclr1, det_beamstatus,
#                      sclr1_ch2, sclr1_ch3, sclr1_ch4,
#                      ssx_rbv, ssy_rbv, ssz_rbv, t_base, t_sample, t_vlens,
#                      t_hlens]

# When the scaler is triggered, the xspress3 is externally triggered with this
# MasterDetector:
master_sclr1 = MasterDetector(sclr1, slaves=[xspress3.filestore, merlin1.filestore])
#master_sclr1 = MasterDetector(sclr1)
#master_sclr1 = MasterDetector(sclr1, slaves=[merlin1.filestore])


# NOTE: master_sclr1 has SUB-detectors which are not in this list (see above)
#gs.DETS = [zebra, master_sclr1, sclr1_ch2, sclr1_ch3, sclr1_ch4, ssx, ssy, ssz, t_base, t_sample,
#           t_vlens, t_hlens, xbpm_x, xbpm_y, quad_x, quad_y, dcm_th, dcm_p, angle_x, angle_y, slit1_top, slit1_bottom,
#           slit1_right, slit1_left, slit1_xpos, slit1_ypos]
gs.DETS = [zebra,master_sclr1, sclr1_ch2, sclr1_ch3, sclr1_ch4, sclr1_ch4_calc, ssx, ssy, ssz, t_base, t_sample,
           t_vlens, t_hlens, xbpm_x, xbpm_y, dcm_th, dcm_p, angle_x, angle_y, slit1_top, slit1_bottom,
           slit1_right, slit1_left, slit1_xpos, slit1_ypos, int_zpssx, int_zpssy, int_zpssz]


gs.TABLE_COLS = ['sclr1_ch2','sclr1_ch3', 'sclr1_ch4', sclr1_ch4_calc, 'ssx', 'ssy', 'ssz',
                 't_base', 't_sample', 't_vlens', 't_hlens', 'int_zpssx', 'int_zpssy', 'int_zpssz']

for roi in xspress3.rois.get_epics_rois(channels=[1, 2, 3]):
    # Add the ROI to the list of detectors to be recorded
    gs.DETS.append(roi)

    if roi.channel == 1:
        # Add the ROI to the list of columns to be printed in the table:
        gs.TABLE_COLS.append(roi.name)


# Plot this by default versus motor position:
gs.PLOT_Y = 'Det2_V'
