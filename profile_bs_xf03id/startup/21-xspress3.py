from hxntools.detectors import Xspress3Detector


xspress3 = Xspress3Detector('XF:03IDC-ES{Xsp:1}:', cam='', files=['HDF5:'],
                            name='xspress3',
                            file_path='/xspress3_data/',
                            ioc_file_path='/xspress3_data/',
                            default_channels=[1, 2, 3],
                            num_roi=16, channel_prefix='Det')


def xspress3_roi_setup():
    # ROIs (added to all default channels if unspecified)
    xspress3.clear_all_rois()
    xspress3.add_roi(0, 1000, 'BadROI')
    xspress3.add_roi(1340, 1640, 'Al')
    xspress3.add_roi(1590, 1890, 'Si')
    xspress3.add_roi(2150, 2450, 'S')
    xspress3.add_roi(2810, 3110, 'Ar')
    xspress3.add_roi(3540, 3840, 'Ca')
    xspress3.add_roi(4360, 4660, 'Ti')
    xspress3.add_roi(4800, 5100, 'V')
    xspress3.add_roi(5270, 5570, 'Cr')
    xspress3.add_roi(5750, 6050, 'Mn')
    xspress3.add_roi(6250, 6550, 'Fe')
    xspress3.add_roi(6780, 7080, 'Co')
    xspress3.add_roi(7330, 7630, 'Ni')
    xspress3.add_roi(7900, 8200, 'Cu')
    xspress3.add_roi(8490, 8790, 'Zn')
    xspress3.add_roi(1970, 2270, 'Au')
    xspress3.add_roi(9300, 9600, 'Pt')
    xspress3.add_roi(4610, 5070, 'Ce')
    xspress3.add_roi(2000, 2300, 'Zr')
    xspress3.add_roi(1900, 2000, 'Y')
    xspress3.add_roi(6530, 6940, 'Gd')


try:
    print('Configuring Xspress3 ROIs...')
    xspress3_roi_setup()
    print('Done')
except KeyboardInterrupt:
    print('Xspress3 ROI configuration cancelled.')
