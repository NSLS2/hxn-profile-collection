print(f"Loading {__file__!r} ...")
from collections import OrderedDict

from epics import caput, caget
import os
import threading
import h5py
from ophyd.sim import NullStatus
import numpy as np
import time as ttime
from ophyd import Device, EpicsSignal, EpicsSignalRO
from ophyd import Component as Cpt
from bluesky.plan_stubs import (kickoff, collect,
                                complete, abs_set, mv, checkpoint)
from hxntools.detectors.zebra import Zebra, EpicsSignalWithRBV
from databroker.assets.handlers import HandlerBase
from ophyd.areadetector.filestore_mixins import resource_factory
from bluesky.preprocessors import (stage_decorator,
                                   run_decorator, subs_decorator,
                                   monitor_during_decorator, finalize_wrapper)

from event_model import compose_resource

class ROIHDF5Handler(HandlerBase):
    HANDLER_NAME = "ROI_HDF5_FLY"

    def __init__(self, resource_fn):
        self._handle = h5py.File(resource_fn, "r", libver='latest', swmr=True)

    def __call__(self, *, det_elem):
        ds = self._handle[det_elem]
        ds.id.refresh()
        return ds[:]

    # def close(self):
    #     self._handle.close()
    #     self._handle = None
    #     super().close()


db.reg.register_handler(ROIHDF5Handler.HANDLER_NAME, ROIHDF5Handler, overwrite=True)

class ExportXpsROI:
    def __init__(self):
        self._fp = None
        self._filepath = None

    def open(self, filepath, xspress3):
        self.close()
        self._filepath = filepath
        self._fp = h5py.File(filepath, "w", libver="latest")

        self._fp.swmr_mode = True

        self._xsp = xspress3

        def create_ds(det_name):
            if not det_name in self._fp:
                self._fp.create_dataset(det_name, data=np.array([], dtype="f"), maxshape=(None,), dtype="f")

        for det_name in [roi.name for roi in self._xsp.enabled_rois]:
            create_ds(det_name)

        self._fp.flush()

    def close(self):
        if self._fp:
            self._fp.close()
            self._fp = None

    def __del__(self):
        self.close()

    def export(self, npoints):
        def add_data(det_name, data):
            ds = self._fp[det_name]
            ds.resize((npoints,))
            ds[:len(data)] = np.array(data)
            ds[len(data):] = data[-1]

        for roi in self._xsp.enabled_rois:
            if hasattr(roi,'settings'):
                add_data(roi.name, roi.settings.array_data.get())
            else:
                add_data(roi.name, roi.ts_total.get())
        self._fp.flush()

class PandaLivePlot():
    def __init__(self):
        self.do_plot = False
        self.fig = plt.figure(123)
        self.scan_id = 0

    def setup_plot(self,scan_input,det):
        self.fig.clear()
        self.ax = self.fig.add_subplot()
        self.xsp = det
        self.scan_input = scan_input.copy()
        self.total_points = scan_input[2] * scan_input[5]
        self.do_plot = True

    def update_plot(self, finished = False):
        if not self.do_plot:
            return
        #Live plot
        n = 0
        for roi in self.xsp.enabled_rois:
            if hasattr(roi,'settings'):
                fluo_data = roi.settings.array_data.get()
            else:
                fluo_data = roi.ts_total.get()
            if not finished:
                fluo_data = np.pad(fluo_data,(0,int(self.total_points-len(fluo_data))))
            else:
                fluo_data = np.pad(fluo_data,(0,int(self.total_points-len(fluo_data))),'edge')
                self.do_plot = False

            self.ax.set_title(roi.name)
            self.ax.imshow(fluo_data.reshape((int(self.scan_input[5]),int(self.scan_input[2]))),extent = [self.scan_input[0],self.scan_input[1],self.scan_input[4],self.scan_input[3]])
            self.ax.set_aspect('equal','box')
            n +=1
            if n>0:
                break
        self.fig.canvas.manager.set_window_title('Scan %d'%(self.scan_id))
        self.fig.canvas.manager.show()
        self.fig.canvas.draw()
        self.fig.canvas.flush_events()

panda_live_plot = PandaLivePlot()

plt.close(panda_live_plot.fig)
