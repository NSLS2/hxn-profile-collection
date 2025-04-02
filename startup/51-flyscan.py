print(f"Loading {__file__!r} ...")

# vim: sw=4 ts=4 sts=4 expandtab smarttab
# HXN fly-scan configuration
from hxnfly.bs import (FlyPlan1D, FlyPlan2D, FlyStep1D, maybe_a_table)
from hxnfly.hxn_fly import (Fly1D_MLL, Fly1D_ZP, Fly2D_MLL, Fly2D_ZP,
        Fly1D_Diffraction, Fly2D_Diffraction)

# These define which scans can be done by which class:
# 1D scans:
FlyPlan1D.scans = {frozenset({smll.ssx}): Fly1D_MLL,
                   frozenset({smll.ssy}): Fly1D_MLL,
                   frozenset({smll.ssz}): Fly1D_MLL,
                   frozenset({zps.zpssx}): Fly1D_ZP,
                   frozenset({zps.zpssy}): Fly1D_ZP,
                   frozenset({zps.zpssz}): Fly1D_ZP,
                   frozenset({smlld.dssx}): Fly1D_Diffraction,
                   frozenset({smlld.dssy}): Fly1D_Diffraction,
                   frozenset({smlld.dssz}): Fly1D_Diffraction,
                   }

# Mixed 1D fly/1D step scan, same as above:
FlyStep1D.scans = FlyPlan1D.scans

# 2D scans:
FlyPlan2D.scans = {frozenset({smll.ssx, smll.ssy}): Fly2D_MLL,
                   frozenset({smll.ssx, smll.ssz}): Fly2D_MLL,
                   frozenset({smll.ssy, smll.ssz}): Fly2D_MLL,

                   frozenset({zps.zpssx, zps.zpssy}): Fly2D_ZP,
                   frozenset({zps.zpssx, zps.zpssz}): Fly2D_ZP,
                   frozenset({zps.zpssy, zps.zpssz}): Fly2D_ZP,

                   frozenset({smlld.dssx, smlld.dssy}): Fly2D_Diffraction,
                   frozenset({smlld.dssx, smlld.dssz}): Fly2D_Diffraction,
                   frozenset({smlld.dssy, smlld.dssz}): Fly2D_Diffraction,
                   }


# from hxnfly.callbacks import FlyLivePlot
from hxnfly.callbacks import (FlyRoiPlot, FlyLiveImage)
from hxnfly.callbacks import FlyLiveCrossSection

# def _sum_func(*values):
#     return np.sum(values, axis=0)
#
# A live plot using scaler 1's first three MCA channels and sums them
# together:
# flyplot = FlyLivePlot(sclr1_mca[:3], data_func=_sum_func,
#                       labels=['test'])
#
# A 1D ROI plot of aluminum from channels 1 to 3:
# flyplot = FlyRoiPlot('Al', channels=[1, 2, 3])

#live_im_plot = FlyLiveImage(['Ca','W_L','Fe','Pt_L'], channels=[1, 2, 3])

#live_im_plot = FlyLiveImage(['Fe','Ca','Cl','S'], channels=[1, 2, 3])
live_im_plot = FlyLiveImage(live_plot_elems, channels=[1, 2, 3])
# fly2dplot1 = FlyLiveCrossSection(['V'], channels=[1, 2, 3)

#pt_plot = FlyRoiPlot(['Cr'],
#                     channels=[1, 2, 3],
#                     )

pt_plot = FlyRoiPlot(line_plot_elem,
                     channels=[1, 2, 3],
                    )

'''
pt_plot = FlyRoiPlot(['Fe'],
                     channels=[1, 2, 3],
                    )
'''
# NOTE: indicate which detectors can be used in fly scans.
# fly_scannable_detectors = [xspress3, zebra, sclr1, dexela1]
fly_scannable_detectors = [xspress3, zebra, sclr1]
_fly1d = FlyPlan1D(usable_detectors=fly_scannable_detectors,
                  scaler_channels=range(1,17))

_fly1d.sub_factories = [maybe_a_table]
_fly1d.subs = [pt_plot, ]
fly1d = _fly1d.__call__

_fly2d = FlyPlan2D(usable_detectors=fly_scannable_detectors,
                  scaler_channels=range(1,17))
_fly2d.sub_factories = [maybe_a_table]
_fly2d.subs = [pt_plot, live_im_plot, ]

flystep = FlyStep1D(usable_detectors=fly_scannable_detectors,
                    scaler_channels=[1, 2, 3, 4, 5, 6, 7, 8])
flystep.sub_factories = [maybe_a_table]
flystep.subs = [pt_plot, live_im_plot, ]

def fly2d(dets,
          motor1, scan_start1, scan_end1, num1,
          motor2, scan_start2, scan_end2, num2,
          exposure_time,
          *,
          dead_time=0.004, fly_type=None,
          return_speed=None, max_points=None, md=None, pulse_dur = 1):

      for det in dets:
            if det.name in ['eiger1']:
                  det.hdf5.frame_per_point = num1*num2
                  yield from bps.abs_set(det.cam.num_triggers, num1*num2, wait=True)
            if det.name == 'xspress3':
                det.hdf5.filestore_spec = 'XSP3'
      #zebra.stage_sigs.update([(zebra.pulse[1].width,exposure_time*pulse_dur)])
      yield from _fly2d(
        dets=dets,
        motor1=motor1, scan_start1=scan_start1, scan_end1=scan_end1, num1=num1,
        motor2=motor2, scan_start2=scan_start2, scan_end2=scan_end2, num2=num2,
        exposure_time=exposure_time,
        dead_time=dead_time, fly_type=fly_type,
        return_speed=return_speed, max_points=max_points, md=md
  )
