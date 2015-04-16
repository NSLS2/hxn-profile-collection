# vim: sw=4 ts=4 sts expandtab smartab

from ophyd.userapi.scan_api import Scan, AScan, DScan, Count

scan = Scan()
ascan = AScan()
ascan.default_triggers = [sclr_trig]
ascan.default_detectors = [ion0, ion1, ionN, ion3]
dscan = DScan()

dscan.default_detectors = ascan.default_detectors

# Use ct as a count which is a single scan.

ct = Count()

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
    args = _elements_become_singlets(args)
    args.append(steps)
    return dscan(*args)


def _elements_become_singlets(a):
    """ [[a, b], c] -> [[[a], [b]], [c]] """
    return [_elements_become_singlets(element) if isinstance(element, (list, tuple))
            else [element]
            for element in a]

def movr_hth(angle):
    radian=angle*pi/180.0
    correction=-1.*tan(radian)*34376.6
    movr(hth, angle)
    movr(hx,correction)
#    return movr_hth(angle)

    
