from hxntools.detectors import Xspress3Detector


xspress3 = Xspress3Detector('XF:03IDC-ES{Xsp:1}:', cam='', files=['HDF5:'],
                            name='xspress3',
                            file_path='/xspress3_data/',
                            ioc_file_path='/xspress3_data/',
                            default_channels=[1, 2, 3],
                            num_roi=16, channel_prefix='Det')


def xspress3_roi_setup():
    # ROIs (added to all default channels if unspecified)
    rois = xspress3.rois
    # rois.clear_all()
    rois.add(9300, 9600, 'Pt')
    rois.add(1590, 1890, 'Si')
    rois.add(2150, 2450, 'S')
    rois.add(2000, 2300, 'Zr')
    rois.add(2810, 3110, 'Ar')
    rois.add(3540, 3840, 'Ca')
    rois.add(4610, 5070, 'Ce')
    rois.add(4000, 4200, 'Xe')
    rois.add(5270, 5570, 'Cr')
    rois.add(5750, 6050, 'Mn')
    rois.add(6250, 6550, 'Fe')
    rois.add(6530, 6940, 'Gd')
    rois.add(6780, 7080, 'Co')
    rois.add(7330, 7630, 'Ni')
    rois.add(7900, 8200, 'Cu')
    rois.add(8250, 8550, 'W')
    rois.add(8490, 8790, 'Zn')
    rois.add(9600, 9750, 'Au')
    rois.add(11500, 12500,'EL' )
    rois.add(1900, 2000, 'Y')
    rois.add(1340, 1640, 'Al')
    rois.add(4360, 4660, 'Ti')
    rois.add(4550, 4750, 'La')

try:
    print('Configuring Xspress3 ROIs...')
    xspress3_roi_setup()
    print('Done')
except KeyboardInterrupt:
    print('Xspress3 ROI configuration cancelled.')
