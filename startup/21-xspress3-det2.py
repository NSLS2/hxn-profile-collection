print(f"Loading {__file__!r} ...")

from ophyd.device import (Component as Cpt)

from nslsii.detectors.xspress3 import (Xspress3FileStore,
                                       Xspress3Channel)
from hxntools.detectors.hxn_xspress3 import HxnXspress3DetectorBase
import threading
from ophyd import DeviceStatus
from ophyd.areadetector import Xspress3Detector
from nslsii.areadetector.xspress3 import (
    build_xspress3_class,
    Xspress3HDF5Plugin,
    Xspress3Trigger,
    Xspress3FileStore,
)

# this is the community IOC package
from nslsii.areadetector.xspress3 import (
    build_xspress3_class
)

xspress3_channel_numbers=(1, 2, 3)

num_elem = np.size(roi_elems)
if num_elem > 48:
    num_elem = 48

class Xspress3HDF5PluginHXN(Xspress3HDF5Plugin):
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)
        self.kind = 1
    def describe(self):
        desc = super().describe()

        spec = {
            "external": "FileStore:",
            "dtype": "array",
            # TODO do not hard code
            "shape": (self.parent.cam.num_images.get(), 4, 4096),
            "source": self.prefix,
        }
        return {'xspress3': spec}


CommunityXspress3_8Channel = build_xspress3_class(
    channel_numbers=xspress3_channel_numbers,
    mcaroi_numbers=tuple(i for i in range(1,num_elem+1)),
    image_data_key="fluor",
    xspress3_parent_classes=(Xspress3Detector, Xspress3Trigger),
    extra_class_members={
        "hdf5": Cpt(
            Xspress3HDF5PluginHXN,
            "HDF1:",
            name="hdf5",
            resource_kwargs={},
            # These are overriden by properties.
            path_template='/nsls2/data/hxn/legacy/%Y/%m/%d/',
            root_path='/nsls2/data/hxn/legacy',
        )
    }

)
class CommunitySrxXspress3Detector(CommunityXspress3_8Channel):
    # replace HDF5:FileCreateDir with HDF1:FileCreateDir
    create_dir = Cpt(EpicsSignal, "HDF1:FileCreateDir")
    erase = Cpt(EpicsSignal, "det1:ERASE")

    # this is used as a latch to put the xspress3 into 'bulk' mode
    # for fly scanning.  Do this is a signal (rather than as a local variable
    # or as a method so we can modify this as part of a plan
    fly_next = Cpt(Signal, value=False)

    def __init__(
        self,
        prefix,
        *,
        f_key="fluor",
        configuration_attrs=None,
        read_attrs=None,
        **kwargs,
    ):
        self._f_key = f_key
        if configuration_attrs is None:
            configuration_attrs = [
                "external_trig",
                "total_points",
                "spectra_per_point",
                "cam",  # replaced settings with cam
                "rewindable",
            ]
        super().__init__(
            prefix,
            configuration_attrs=configuration_attrs,
            read_attrs=read_attrs,
            **kwargs,
        )

        if read_attrs is None:
            pass
            # JL removed channels from read_attrs
            #read_attrs = ["channel1", "channel2", "channel3", "channel4", "hdf5"]
            # JL read all rois on all channels
            #read_attrs = [
            #    xs_mcaroi.total_rbv.name
            #    for xs_channel
            #    in self.iterate_channels()
            #    for xs_mcaroi
            #    in xs_channel.iterate_mcarois()
            #]
            #read_attrs.append("hdf5")
            #print(f"read_attrs: {read_attrs}")
        # this is possiblely one too many places to store this
        # in the parent class it looks at if the extrenal_trig signal is high
        self._mode = 1
        self.channel_num = self.channel_numbers


        # 2020-01-24
        # Commented out by AMK for using the xs3-server-IOC from TES
        # self.create_dir.put(-3)
    def read(self):
        res = super().read()
        res.update({'xspress3':{'value':xspress3_det2.hdf5._datum_factory(datum_kwargs={})['resource'] + '/0'}})
        return res

    def stop(self, *, success=False):
        ret = super().stop()
        # todo move this into the stop method of the settings object?
        self.cam.acquire.put(0)
        self.hdf5.stop(success=success)
        return ret

    def _compute_total_capture(self):
        total_points = self.total_points.get()
        if total_points < 1:
            raise RuntimeError("You must set the total points")
        spec_per_point = self.spectra_per_point.get()
        total_capture = total_points * spec_per_point
        return total_points, spec_per_point, total_capture

    def stage(self):
        # if should external trigger
        ext_trig = True

        # really force it to stop acquiring
        self.cam.acquire.put(0, wait=True)
        self.erase.put(1)

        total_points, spec_per_point, total_capture = self._compute_total_capture()

        for channel in self.iterate_channels():
            channel.mcaroi.ts_control.put(0)
            channel.mcaroi.ts_num_points.put(total_points)
        # # stop previous acquisition
        # self.stage_sigs[self.parent.cam.acquire] = 0

        # # re-order the stage signals and disable the calc record which is
        # # interfering with the capture count
        # self.stage_sigs.pop(self.num_capture, None)
        # self.stage_sigs.pop(self.parent.cam.num_images, None)
        # self.stage_sigs[self.num_capture_calc_disable] = 1

        if ext_trig:
            self.stage_sigs[self.cam.trigger_mode] = 'TTL Veto Only'
            self.stage_sigs[self.cam.num_images] = total_capture
        else:
            # self.settings.trigger_mode.put('Internal')
            # self.settings.num_images.put(1)
            self.stage_sigs[self.cam.trigger_mode] = 'Internal'
            self.stage_sigs[self.cam.num_images] = spec_per_point
            # Failed attempt to fix expected shape in tiled

        self.stage_sigs[self.hdf5.auto_save] = 'Yes'

        # print("stage!")
        # Erase what is currently in the system
        # This prevents a single hot pixel in the upper-left corner of a map
        # JL replaced xs.erase.put(0) with self.cam.erase.put(0)
        #    why was xs.erase.put(0) not self.erase.put(0) ?
        #xs.erase.put(0)
        # JL commented out the next line because it caused a significant delay in starting acqusitions
        #self.cam.erase.put(0)
        # JL added the next line, it is not pretty

        # file_write_mode = self.hdf5.file_write_mode.get(as_string=True)
        # print(f"{file_write_mode = }")

        # self.previous_file_write_mode_value = self.hdf5.file_write_mode.get()
        # JL added the next 2 lines
        #   should use stage_sigs for file_write_mode?
        # self.hdf5.file_write_mode.put(1)
        #self.hdf5.auto_save.put(1)  # using stage_sigs for this
        # do the latching
        if self.fly_next.get():
            self.fly_next.put(False)
            self._mode = 2
        return super().stage()

    def unstage(self):
        # print("unstage!")
        # JL added the next two lines
        #self.hdf5.auto_save.put(0)
        # self.hdf5.file_write_mode.put(self.previous_file_write_mode_value)
        # JL removed the next line
        #self.hdf5.capture.put(0)  # this PV un-sets itself
        try:
            ret = super().unstage()
        finally:
            self._mode = 1
        return ret
    @property
    def enabled_rois(self):
        for ch in self.iterate_channels():
            yield from ch.iterate_mcarois()

xspress3_det2 = CommunitySrxXspress3Detector("XF:03IDC-ES{Xsp:3}:", name="xspress3_det2", f_key="fluor")

def xspress3_det2_roi_setup(det):
    from hxntools.elem_fluo_lines import elem_fluo_lines
    for channel in det.iterate_channels():
        i = 0
        for roi in channel.iterate_mcarois():
            elem = roi_elems[i]
            try:
                energy = elem_fluo_lines[elem]
                roi.configure_roi((energy-150)//10,(energy+150)//10)
                roi.name = 'Det'+str(channel.channel_number) + '_'+ elem
            except:
                pass
            i += 1
try:
    xspress3_det2_roi_setup(xspress3_det2)
except:
    pass
