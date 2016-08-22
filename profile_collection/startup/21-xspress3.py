from ophyd.device import (Component as Cpt)

from hxntools.detectors.xspress3 import (Xspress3FileStore,
                                         Xspress3Channel)
from hxntools.detectors.hxn_xspress3 import HxnXspress3DetectorBase


class HxnXspress3Detector(HxnXspress3DetectorBase):
    channel1 = Cpt(Xspress3Channel, 'C1_', channel_num=1)
    channel2 = Cpt(Xspress3Channel, 'C2_', channel_num=2)
    channel3 = Cpt(Xspress3Channel, 'C3_', channel_num=3)
    # Currently only using three channels. Uncomment these to enable more
    # channels:
    # channel4 = C(Xspress3Channel, 'C4_', channel_num=4)
    # channel5 = C(Xspress3Channel, 'C5_', channel_num=5)
    # channel6 = C(Xspress3Channel, 'C6_', channel_num=6)
    # channel7 = C(Xspress3Channel, 'C7_', channel_num=7)
    # channel8 = C(Xspress3Channel, 'C8_', channel_num=8)

    hdf5 = Cpt(Xspress3FileStore, 'HDF5:',
               write_path_template='/data/%Y/%m/%d/',
               mds_key_format='xspress3_ch{chan}',
               )

    def __init__(self, prefix, *, configuration_attrs=None, read_attrs=None,
                 **kwargs):
        if configuration_attrs is None:
            configuration_attrs = ['external_trig', 'total_points',
                                   'spectra_per_point']
        if read_attrs is None:
            read_attrs = ['channel1', 'channel2', 'channel3', 'hdf5']
        super().__init__(prefix, configuration_attrs=configuration_attrs,
                         read_attrs=read_attrs, **kwargs)


xspress3 = HxnXspress3Detector('XF:03IDC-ES{Xsp:1}:', name='xspress3')


# Create directories on the xspress3 server, otherwise scans can fail:
xspress3.make_directories.put(True)


def xspress3_roi_setup():
    for channel in [xspress3.channel1, xspress3.channel2, xspress3.channel3]:
        #channel.set_roi(1, 9300, 9600, name='Pt')
        channel.set_roi(1, 1590, 1890, name='Si')
        channel.set_roi(2, 2150, 2450, name='S')
        #channel.set_roi(2, 4690, 4990, name='Ce')
        #channel.set_roi(3, 5750, 6050, name='Mn')
        channel.set_roi(3, 4360, 4660, name='Ti')
        channel.set_roi(4, 3600, 4200, name='Te')
        channel.set_roi(5, 2810, 3110, name='Ag')
        channel.set_roi(6, 6780, 7080, name='Co')
        #channel.set_roi(7, 3542, 3842, name='Ca')
        channel.set_roi(7, 4800, 5100, name='V')
        channel.set_roi(8, 1850, 2150, name='P')
        channel.set_roi(9, 5270, 5570, name='Cr')
        channel.set_roi(10, 3160, 3460, name='K')
        channel.set_roi(11, 6250, 6550, name='Fe')
        channel.set_roi(12, 6530, 6940, name='Gd')
        channel.set_roi(13, 8487, 8787, name='Zn')
        channel.set_roi(14, 7330, 7630, name='Ni')
        channel.set_roi(15, 7900, 8200, name='Cu')
        channel.set_roi(16, 9150, 9350, name='Ga')
        #channel.set_roi(16, 9736, 10036, name='Ge')
        # channel.set_roi(17, 8250, 8550, 'W')
        # channel.set_roi(18, 9600, 9750, 'Au')
        # channel.set_roi(19, 11500, 12500, 'EL')
        # channel.set_roi(20, 1900, 2000, 'Y')
        # channel.set_roi(21, 1340, 1640, 'Al')
        # channel.set_roi(22, 4360, 4660, 'Ti')
        # channel.set_roi(23, 4550, 4750, 'La')
        # channel.set_roi(24, 9150, 9350, 'Ga')


try:
    print('Configuring Xspress3 ROIs...')
    xspress3_roi_setup()
    print('Done')
except KeyboardInterrupt:
    print('Xspress3 ROI configuration cancelled.')
