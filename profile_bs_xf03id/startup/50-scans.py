# vim: sw=4 ts=4 sts expandtab smarttab
# HXN step-scan configuration

import hxntools.scans
from hxntools.detectors.master_detector import MasterDetector
from hxntools.spiral_scans import HxnFermatScan


# TODO: figure out olog issues
olog_client = None
# TODO

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

# NOTE: master_sclr1 has SUB-detectors which are not in this list (see above)
gs.DETS = [zebra, master_sclr1, tpx1_roi, sclr2_ch2,sclr2_ch3, sclr2_ch4, ssx, ssy, ssz, t_base, t_sample,
           t_vlens, t_hlens, timepix1.filestore]

gs.TABLE_COLS = ['tpx1_roi','sclr2_ch2','sclr2_ch3', 'sclr2_ch4', 'ssx', 'ssy', 'ssz',
                 't_base', 't_sample', 't_vlens', 't_hlens']

for roi in xspress3.rois.get_epics_rois(channels=[1, 2, 3]):
    # Add the ROI to the list of detectors to be recorded
    gs.DETS.append(roi)

    if roi.channel == 1:
        # Add the ROI to the list of columns to be printed in the table:
        gs.TABLE_COLS.append(roi.name)


# Plot this by default versus motor position:
gs.PLOT_Y = 'Det2_Pt'
