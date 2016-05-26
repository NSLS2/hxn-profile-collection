# vim: sw=4 ts=4 sts=4 expandtab smarttab
# HXN fly-scan configuration
from hxnfly.bs import (FlyPlan1D, FlyPlan2D, FlyStep1D, maybe_a_table)
from hxnfly.hxn_fly import (Fly1D_MLL, Fly1D_ZP, Fly2D_MLL, Fly2D_ZP)

# These define which scans can be done by which class:
# 1D scans:
FlyPlan1D.scans = {frozenset({ssx}): Fly1D_MLL,
                   frozenset({ssy}): Fly1D_MLL,
                   frozenset({ssz}): Fly1D_MLL,
                   frozenset({zpssx}): Fly1D_ZP,
                   frozenset({zpssy}): Fly1D_ZP,
                   frozenset({zpssz}): Fly1D_ZP,
                   }

# Mixed 1D fly/1D step scan, same as above:
FlyStep1D.scans = FlyPlan1D.scans

# 2D scans:
FlyPlan2D.scans = {frozenset({ssx, ssy}): Fly2D_MLL,
                   frozenset({ssx, ssz}): Fly2D_MLL,
                   frozenset({ssy, ssz}): Fly2D_MLL,
                   frozenset({zpssx, zpssy}): Fly2D_ZP,
                   frozenset({zpssx, zpssz}): Fly2D_ZP,
                   frozenset({zpssy, zpssz}): Fly2D_ZP,
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
#
live_im_plot = FlyLiveImage(['V','P','Ag','Pt','Si'], channels=[1, 2, 3])
# fly2dplot1 = FlyLiveCrossSection(['V'], channels=[1, 2, 3)

pt_plot = FlyRoiPlot(['Pt', 'V', 'Ag'], channels=[1, 2, 3])


# NOTE: indicate which detectors can be used in fly scans. When a
#       fly scan is run, all of those matching in gs.DETS will be
#       used.
fly_scannable_detectors = [xspress3, merlin1, zebra, sclr1]
fly1d = FlyPlan1D(usable_detectors=fly_scannable_detectors,
                  scaler_channels=[1, 2, 3, 4],
                  )

fly1d.sub_factories = [maybe_a_table]
fly1d.subs = [pt_plot, ]

fly2d = FlyPlan2D(usable_detectors=fly_scannable_detectors,
                  scaler_channels=[1, 2, 3, 4])
fly2d.sub_factories = [maybe_a_table]
fly2d.subs = [pt_plot, live_im_plot, ]

flystep = FlyStep1D(usable_detectors=fly_scannable_detectors,
                    scaler_channels=[1, 2, 3, 4])
flystep.sub_factories = [maybe_a_table]
flystep.subs = [pt_plot, live_im_plot, ]
