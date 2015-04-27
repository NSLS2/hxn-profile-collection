# vim: sw=4 ts=4 sts expandtab smarttab

import os
import time

from ophyd.userapi.scan_api import Scan, AScan, DScan, Count
from hxntools.ophyd_tools import HXNDScan, HXNAScan


# Use ct as a count which is a single scan.
ct = Count()


def setup_xrfscan(scan):
    # TODO: wrap this up in a subclassed SignalDetector
    # Stop acquiring
    xrf_acquire.put(0)
    xrf_erase.put(1)
    
    nfs_path = '/data/%s/scan%.5d/' % (time.strftime('%Y%m%d'), scan.scan_id)

    if False:  # when software trigger decides to work...
        xrf_acquire.put(1)
        xrf_save.put(1)
        xrf_num_images.put(scan.npts + 1)
    else:
        # internal triggering
        xrf_num_images.put(1)
        xrf_trig_mode.put('Internal')
        xrf_path = os.path.join(nfs_path, 'xspress')
        try:
            os.makedirs(xrf_path)
        except OSError:
            pass

        xrf_filepath.put(xrf_path)
        xrf_filename.put('scan_%.5d_' % scan.scan_id)
        xrf_filenumber.put(1)
   
    tpx1_path = os.path.join(nfs_path, 'timepix1')

    try:
        os.makedirs(tpx1_path)
    except OSError:
        pass

    tpx1_filepath.put(tpx1_path)
    tpx1_filename.put('scan_%.5d' % scan.scan_id)
    tpx1_filenumber.put(0)
    tpx1_autosave.put(1)


def teardown_xrfscan(scan):
    tpx1_autosave.put(0)


class HXN_XRFDScan(HXNDScan):
    def pre_scan(self):
        super(HXNDScan, self).pre_scan()
        setup_xrfscan(self)
    
    def post_scan(self):
        super(HXNDScan, self).post_scan()
        teardown_xrfscan(self)



class HXN_XRFAScan(HXNAScan):
    def pre_scan(self):
        super(HXNAScan, self).pre_scan()
        setup_xrfscan(self)
    
    def post_scan(self):
        super(HXNAScan, self).post_scan()
        teardown_xrfscan(self)


# scan = Scan()
ascan = HXN_XRFAScan()
ascan.default_detectors = [det_xrf_erase, det_xrf_save, det_xrf_acquire, det_sclr2, det_tpx1_acquire, 
                           ion0, ion1, ion3, Pt_ch1, Pt_ch2, Pt_ch3, ssx_rbv, ssy_rbv, ssz_rbv, tpx1_filenumber]

# example with saving:
# ascan.default_detectors = [det_xrf_save, det_xrf_acquire, det_sclr2, 
#                            ion0, ion1, ionN, ion3, Pt_ch1, Pt_ch2, Pt_ch3]
# example without saving:
# ascan.default_detectors = [det_xrf_acquire, det_sclr2, 
#                            ion0, ion1, ionN, ion3, Pt_ch1, Pt_ch2, Pt_ch3]

dscan = HXN_XRFDScan()
dscan.default_detectors = ascan.default_detectors


def synchronize(detectors, integration_time):
    # EpicsScaler calls it preset_time;
    # AreaDetector calls it exposure_time.
    attrs = ['preset_time', 'exposure_time']
    for det in detectors:
        for attr in attrs:
            try:
                setattr(det, attr, integration_time)
            except AttributeError:
                pass
            # TODO Raise better.


def sync_dscan(positioners, start, stop, step, acquisition_time):
    """Perform a normal dscan, but first sync acquisition time.

    Parameters
    ----------
    positioners : same as dscan
    start : same as dscan
    stop : same as dscan
    step : same as dscan
    acquisition_time : sensor integration time in seconds
    """
    synchronize(dscan.detectors, acquisition_time)
    return dscan(positioners, start, stop, step)


def d2_scan(*args):
    """
    Parameters
    ----------
    args : (positioner, start, stop)
    steps : number of steps in the scan

    Example
    -------
    # In five steps, move sli from 0 to 5 and slo from 0 to 10.
    >>> d2_scan([sli, 0, 5], [slo, 0, 10], 5)
    """

    args = list(args)
    steps = args.pop()  # The last argument is the number of steps.
    args = [[list(lst)] for lst in zip(*args)]
    args.append(steps)
    return dscan(*args)


def mesh_scan(*args):
    '''
    Parameters
    ----------
    args: (positioner, start, stop)
    steps: [size1, size2, ...]

    Example
    -------
    >>> mesh([sli, 0, 5], [slo, 0, 10], [5, 5])
    >>> mesh([sli, 0, 5], [slo, 0, 10], [slt, 0, 10], [5, 5, 5])
    '''
    args = list(args)
    steps = args.pop()  # The last argument is the number of steps.
   
    args = [list(lst) for lst in zip(*args)]
    args = _elements_to_singlets(args)
    args.append(steps)
    return dscan(*args)


def _elements_to_singlets(a):
    """ [[a, b], c] -> [[[a], [b]], [c]] """
    return [_elements_to_singlets(element) if isinstance(element, (list, tuple))
            else [element]
            for element in a]


def movr_hth(angle):
    radian = angle*pi/180.0
    correction = -1.*tan(radian)*34376.6
    movr(hth, angle)
    movr(hx,correction)


def movr_ssx(d):
    dx_n = d * cos(15.*pi/180.)
    dz_n = d * sin(15.*pi/180.)
    movr(ssx, dx_n)
    movr(ssz, dz_n)

def movr_sx(d):
    dx_n = d * cos(15.*pi/180.)
    dz_n = d * sin(15.*pi/180.)
    movr(sx, dx_n)
    movr(sz, dz_n)


def movr_sz(d):
    dx_n = -d * sin(15.*pi/180.)
    dz_n = d * cos(15.*pi/180.)
    movr(sx, dx_n)
    movr(sz, dz_n)


def movr_ssz(d):
    dx_n = -d * sin(15.*pi/180.)
    dz_n = d * cos(15.*pi/180.)
    movr(ssx, dx_n)
    movr(ssz, dz_n)


