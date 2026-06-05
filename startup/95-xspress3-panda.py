print(f"Loading {__file__!r} ...")

import h5py
import numpy as np


#def xspress3_det2_roi_setup(det):
#    from hxntools.elem_fluo_lines import elem_fluo_lines
#    for channel in det.iterate_channels():
#        i = 0
#        for roi in channel.iterate_mcarois():
#            elem = roi_elems[i]
#            try:
#                energy = elem_fluo_lines[elem]
#                roi.configure_roi((energy-150)//10,(energy+150)//10)
#                roi.name = 'Det'+str(channel.channel_number) + '_'+ elem
#            except:
#                pass
#            i += 1
##try:
##    xspress3_det2_roi_setup(xspress3_det2)
##except:
##    pass


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
            if ds.size == 0:
                ds.resize((npoints,))
                ds[:len(data)] = np.array(data)
                ds[len(data):] = ds[len(data)-1]

        for roi in self._xsp.enabled_rois:
            if hasattr(roi,'settings'):
                add_data(roi.name, roi.settings.array_data.get())
            else:
                add_data(roi.name, roi.ts_total.get())
        self._fp.flush()
