# vim: sw=4 ts=4 sts expandtab smarttab
# HXN fly-scan configuration
from hxnfly.bs import (FlyPlan1D, FlyPlan2D)
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

# 2D scans:
FlyPlan2D.scans = {frozenset({ssx, ssy}): Fly2D_MLL,
                   frozenset({ssx, ssz}): Fly2D_MLL,
                   frozenset({ssy, ssz}): Fly2D_MLL,
                   frozenset({zpssx, zpssy}): Fly2D_ZP,
                   frozenset({zpssx, zpssz}): Fly2D_ZP,
                   frozenset({zpssy, zpssz}): Fly2D_ZP,
                   }


# # from hxnfly.callbacks import FlyLivePlot
# from hxnfly.callbacks import (FlyRoiPlot, FlyLiveImage)
# from hxnfly.callbacks import FlyLiveCrossSection
#
# # def _sum_func(*values):
# #     return np.sum(values, axis=1)
# #
# # A live plot using scaler 1's first three MCA channels and sums them
# # together:
# # flyplot = FlyLivePlot(sclr1_mca[:3], data_func=_sum_func,
# #                       labels=['test'])
# #
# # A 1D ROI plot of aluminum from channels 1 to 3:
# # flyplot = FlyRoiPlot('Al', channels=[1, 2, 3])
# # (optionally sum them together with use_sum=True or adding a function to
# #  calculate with data_func=sum_func)
# #
# # flyplot = FlyRoiPlot(['Pt'], channels=[1, 2, 3], use_sum=True)
# # fly2dplot = FlyLiveImage(['V','P','Ag','Pt','Si'], channels=[1, 2, 3], use_sum=True)
# # fly2dplot1 = FlyLiveCrossSection(['V'], channels=[1, 2, 3], use_sum=True)
#
fly1d = FlyPlan1D(detectors=[xspress3],
                  scaler_channels=[1, 2, 3, 4])

fly2d = FlyPlan2D(detectors=[xspress3],
                  scaler_channels=[1, 2, 3, 4])
