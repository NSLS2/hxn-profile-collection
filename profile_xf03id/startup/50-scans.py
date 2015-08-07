# vim: sw=4 ts=4 sts expandtab smarttab

from hxntools.scans import (HXNDScan, HXNAScan, HXNCount)

from hxnfly import fly1d, fly2d

ct = HXNCount()

ascan = HXNAScan()
ascan.default_detectors = [det_sclr1, det_beamstatus,
                           sclr1_ch2, sclr1_ch3, sclr1_ch4,
                           det1_Pt, det2_Pt, det3_Pt, det1_Al, det2_Al, det3_Al, det1_Si, det2_Si,
                           det3_Si, det1_Gd, det2_Gd, det3_Gd, det1_YSZ, det2_YSZ, det3_YSZ, det1_Ca, det2_Ca,
                           det3_Ca, det1_Ti, det2_Ti, det3_Ti, det1_V, det2_V, det3_V, det1_Cr, det2_Cr,
                           det3_Cr, det1_Mn, det2_Mn, det3_Mn, det1_Fe, det2_Fe, det3_Fe, det1_Co,
                           det2_Co, det3_Co, det1_Ni, det2_Ni, det3_Ni, det1_Cu, det2_Cu, det3_Cu,
                           det1_Ce, det2_Ce, det3_Ce, det1_Au, det2_Au, det3_Au,
                           ssx_rbv, ssy_rbv, ssz_rbv, t_base, t_sample, t_vlens, t_hlens]


ascan.user_detectors = [xspress3.filestore, merlin1.filestore]

fly1d.detectors = [xspress3, merlin1, det_sclr1]
fly2d.detectors = [xspress3, merlin1, det_sclr1]


dscan = HXNDScan()
# Detectors are shared among the scans


def synchronize(detectors, integration_time):
    # EpicsScaler calls it preset_time;
    # AreaDetector calls it exposure_time.
    integration_time = float(integration_time)

    attrs = ['preset_time', 'exposure_time']
    for det in detectors:
        for attr in attrs:
            try:
                setattr(det, attr, integration_time)
            except AttributeError:
                pass
            # TODO Raise better.

    # For a step scan using the zebra or scaler 1, the scaler is the master
    # Scaler LNE -> Zebra IN1_TTL -> Pulse1 -> fanout
    if zebra in detectors or det_sclr1 in detectors:
        sclr1._preset_time.put(integration_time)
        zebra.preset_time = integration_time

    for tpx in [timepix1, timepix2]:
        if tpx.filestore in detectors:
            tpx.acquire_time.put(integration_time)
            tpx.acquire_period.put(integration_time)

    for merlin in [merlin1, ]:
        if merlin.filestore in detectors:
            merlin.acquire_time.put(integration_time)
            merlin.acquire_period.put(integration_time + 0.001)


def sync_dscan(positioners, start, stop, step, exposure_time):
    """Perform a normal dscan, but first sync exposure_time time.

    Parameters
    ----------
    positioners : same as dscan
    start : same as dscan
    stop : same as dscan
    step : same as dscan
    exposure_time : sensor integration time in seconds
    """
    synchronize(dscan.detectors, exposure_time)
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


def mesh(*args):
    '''
    Parameters
    ----------
    args: (positioner, start, stop)
    steps: [size1, size2, ...]
    exposure_time: float

    Example
    -------
    >>> mesh([sli, 0, 5], [slo, 0, 10], [5, 5], 1.0)
    >>> mesh([sli, 0, 5], [slo, 0, 10], [slt, 0, 10], [5, 5, 5], 0.5)
    '''
    args = list(args)
    exposure_time = args.pop()
    steps = args.pop()

    args = [list(lst) for lst in zip(*args)]
    args = _elements_to_singlets(args)
    args.append(steps)
    args.append(exposure_time)
    return sync_dscan(*args)


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
