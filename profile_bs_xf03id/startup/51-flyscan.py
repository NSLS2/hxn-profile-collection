# vim: sw=4 ts=4 sts expandtab smarttab
# HXN fly-scan configuration

from hxnfly.bs import (FlyScan1D, FlyScan2D)
# from hxnfly.callbacks import FlyLivePlot
from hxnfly.callbacks import (FlyRoiPlot, FlyLiveImage)
# from hxnfly.callbacks import FlyLiveCrossSection

# def _sum_func(*values):
#     return np.sum(values, axis=1)
#
# A live plot using scaler 1's first three MCA channels and sums them
# together:
# flyplot = FlyLivePlot(sclr1_mca[:3], data_func=_sum_func,
#                       labels=['test'])
#
# A 1D ROI plot of aluminum from channels 1 to 3:
# flyplot = FlyRoiPlot('Al', channels=[1, 2, 3])
# (optionally sum them together with use_sum=True or adding a function to
#  calculate with data_func=sum_func)
#
flyplot = FlyRoiPlot(['BadROI', 'Al'], channels=[1, 2, 3], use_sum=True)
fly2dplot = FlyLiveImage(['BadROI', 'Al'], channels=[1, 2, 3], use_sum=True)
# fly2dplot = FlyLiveCrossSection(['BadROI'], channels=[1, 2, 3], use_sum=True)


fly1d = FlyScan1D([flyplot], detectors=[xspress3])
fly2d = FlyScan2D([fly2dplot], detectors=[xspress3])
