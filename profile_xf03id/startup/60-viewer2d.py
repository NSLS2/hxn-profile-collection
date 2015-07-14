from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from xray_vision.qt_widgets import CrossSectionMainWindow
from xray_vision.backend.mpl.cross_section_2d import CrossSection
from datetime import datetime
import os

def plot2d(scan, name,row,col):
    num_det = 4;
    if name == "all":
        for i in range(num_det):
            det_name = dscan.detectors[i].name
            det = getattr(dscan.data[scan], det_name)
            plt.figure(i)
            data = np.reshape(det, (row, col))
            plt.imshow(data)
    else:
        plt.figure()
        det = getattr(dscan.data[scan], name)
        data = np.reshape(det, (row, col))
        plt.imshow(data,interpolation='None')
    # return data


def dev(scan,namex,namey):
    dety  = getattr(dscan.data[scan], namey)
    d = 3.13559

    if namex == "energy":
        detx = getattr(dscan.data[scan],"dcm_th")
        num_points = len(detx)
        data = np.zeros((num_points-1,2))
        for i in range(num_points-1):
            data[i,1] = (dety[i+1] - dety[i])/(detx[i+1]-detx[i])
            tmp = (detx[i+1] + detx[i])/2 
            data[i,0] = 12.398/(2*d*np.sin(np.pi*(tmp + (-0.0135))/180))  
            
    else:
        detx = getattr(dscan.data[scan],namex)
        num_points = len(detx)
        data = np.zeros((num_points-1,2))
        for i in range(num_points-1):
            data[i,1] = (dety[i+1] - dety[i])/(detx[i+1]-detx[i])
            data[i,0] = (detx[i+1] + detx[i])/2

 
    plt.figure(20)
    plt.plot(data[:,0],data[:,1])
    
    #return data

def plot(namex,namey,norm):
    plt.figure()
    plt.clf()
    if namey == "Pt":
        data1 = dscan.data[-1].Pt_ch1
        data2 = dscan.data[-1].Pt_ch2
        data3 = dscan.data[-1].Pt_ch3
        data = data1 + data2 + data3
    else:
        data = getattr(dscan.data[-1], namey)
    x =  getattr(dscan.data[-1], namex)
    
    if norm != 'None':
        norm_v = getattr(dscan.data[-1], norm)
        plt.plot(x,data/(norm_v+1.e-8))
        plt.plot(x,data/(norm_v+1.e-8),'bo')
    else:
        plt.plot(x,data)
        plt.plot(x,data,'bo')

    plt.show()

def plotfly(namex):
    plt.figure()
    x = fly_data[namex]
    data1 = fly_data['Ch1 [9300:9600]']
    data2 = fly_data['Ch2 [9300:9600]']
    data3 = fly_data['Ch3 [9300:9600]']
    Pt = data1 + data2 + data3
    Pt[0] = Pt[1]
    plt.plot(x,Pt)
    plt.plot(x,Pt,'bo')
    plt.show()


if 'data_cache' not in globals():
    # Don't erase the cache when reloading this module via %run -i
    data_cache = {}


def _scan_starts(df):
    t = df['time']
    # TODO using simple time threshold to determine where 2d scans start/stop
    scan_starts, = np.where(np.diff(t) > 5)
    return [0] + list(scan_starts) + [len(t)]


def _shift_zeros(df, spectrum):
    zeros, = np.where(spectrum == 0)
    print('number of zeros=%d' % len(zeros))
    scan_starts = _scan_starts(df)
    last_start = 0

    spectrum = list(spectrum)
    for zero_i in reversed(zeros):
        # get end of scan index
        next_scan = np.where(zero_i >= scan_starts)[0][-1] + 1
        scan_end = scan_starts[next_scan]
        # shift all spectrum data down
        print('inserting zero at scan_end [idx %d], removing zero '
              'at index %d' % (scan_end, zero_i))
        spectrum.insert(scan_end, 5)
        del spectrum[zero_i]

    return np.array(spectrum, dtype=float)

def _load_scan(scan_id, fill_events=False):
    '''Load scan from databroker by scan id'''

    if scan_id > 0 and scan_id in data_cache:
        df = data_cache[scan_id]
    else:
        hdr = DataBroker[scan_id]
        scan_id = hdr.scan_id
        if scan_id in data_cache:
            df = data_cache[scan_id]
        else:
            data = DataBroker.fetch_events(hdr, fill=fill_events)
            dm = DataMuxer.from_events(data)
            df = dm.to_sparse_dataframe()
            data_cache[scan_id] = df

    return df

#todo change l, h to clim which defaults to 'auto'
def plot2dfly(scan_id, x='ssx[um]', y='ssy[um]', clim=None, 
              fill_events=False, cmap='Oranges', shift_zeros=False, cols=None):
    """Plot the results of a 2d fly scan

    Parameters
    ----------
    scan_id : int
        Any valid input to DataBroker[] or StepScan
    x : str
        The data key that corresponds to the x axis
        Defaults to 'ssx[um]'
    y : str
        The data key that corresponds to the y axis
        Defaults to 'ssy[um]'
    clim : tuple, optional
        formtted as (min, max)
        If None, defaults to min/max of the data
    fill_events : bool, optional
        Fill the events with data from filestore
        Defaults to False (and is much much faster)
    cmap : str, optional
        Defaults to "Oranges"
        The colormap to use. See the pyplot.cm module for valid color maps
    shift_zeros : bool, optional
        Shift all zeros in the scan [NOTE: this is for the detector triggering-
        related bug only]
    cols : int, optional
        The number of columns in the scan. Automatically detected when possible
    """
    
    df = _load_scan(scan_id, fill_events=fill_events)

    title = 'Scan id %s. ' % scan_id
    rois = [ch for ch in list(df) if ch.startswith('Ch')]
    spectrum = np.sum([getattr(df, roi) for roi in rois], axis=0)
    x_data = df[x]
    y_data = df[y]

    if cols is None:
        prev_x = x_data[0]
        for idx, pos in enumerate(x_data[1:]):
            if pos > prev_x:
                prev_x = pos
            else:
                break
        
        cols = idx + 1
        print('Detected scan columns: ', cols)

    if shift_zeros:
        spectrum = _shift_zeros(df, spectrum)

    if clim is None:
        clim = (np.nanmin(spectrum), np.nanmax(spectrum))
    extent = (np.min(x_data), np.max(x_data), np.max(y_data), np.min(y_data))

    # save it to disk
    dt = datetime.utcnow()
    folder = '/data/{}{:0>2}{:0>2}/'.format(dt.year, dt.month, dt.day)
    if not os.path.exists(folder):
        os.makedirs(folder)
    
    try:
        spectrum2 = spectrum.reshape((-1, cols))
    except Exception as ex:
        print('Unable to reshape data to width: %d (%s: %s)' % (cols, ex.__class__.__name__, ex))
        fig = plt.figure()
        ax2 = plt.subplot(111)
    else:
        print('Reshaped to %s' % (spectrum2.shape, ))
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,5))
        fig.set_tight_layout(True)

        im = ax1.imshow(spectrum2, extent=extent, interpolation='None', cmap=cmap, 
                        vmin=clim[0], vmax=clim[1])
        np.savetxt(folder + 'data_scan_'+np.str(scan_id), spectrum2)

        ax1.set_title('IMSHOW. ' + title)

    # create the scatter plot version
    scatter = ax2.scatter(x_data, y_data, c=spectrum, marker='s', s=250, 
                          cmap=getattr(mpl.cm, cmap), linewidths=0, alpha=.8,
                          vmin=clim[0], vmax=clim[1])
    ax2.set_xlim(np.min(x_data), np.max(x_data))
    ax2.set_ylim(np.min(y_data), np.max(y_data))
    ax2.set_title('SCATTER. ' + title)
    ax2.set_aspect('equal')
    fig.colorbar(scatter)
    fig.savefig(folder + 'data_scan_'+np.str(scan_id)+'.png')
    np.savetxt(folder + 'data_x_y_ch_'+np.str(scan_id),
               np.vstack((x_data, y_data, spectrum)).T)
    
    '''
    plt.figure()
    data1 = fly_data['Ch1 [930:960]'].reshape(row, col)
    data2 = fly_data['Ch2 [930:960]'].reshape(row, col)
    data3 = fly_data['Ch3 [930:960]'].reshape(row, col)
    data = data1 + data2 + data3
    data[0,0] = data[0,1]
    plt.imshow(data, interpolation='none')
    plt.show()
    '''


def export(sid):
    data = StepScan[sid]
    np.savetxt('/data/txt/scan_'+np.str(sid)+'.txt',data,fmt='%1.5e')


def _fly2d_points_from_olog(scan_id):
    '''Messy way to get fly2d tile scan points from olog

    Need to improve olog properties in the future to simplify this.
    '''

    cmd = None
    for log in olog_client.find(property='FlyScan'):
        prop = log['properties']['FlyScan']
        if int(prop['id']) == scan_id:
            cmd = prop['command']
            print('Scan command was: %s' % cmd)
            break
    
    if cmd is None:
        raise ValueError('Fly scan not found in olog')

    # ugh, put these in the properties list...
    def _get_args(x_motor, x1, x2, total_pointsx,
                  y_motor, y1, y2, total_pointsy,
                  *args, **kwargs):
        return (x1, x2, total_pointsx,
                y1, y2, total_pointsy)

    point_args = eval(cmd, {'fly2d': _get_args,
                            'x': None, 'y': None, 'z': None})
    from hxnfly.fly import _get_scan_points
    return list(_get_scan_points(*point_args)), log


def flyimshow(scan_id, x='ssx[um]', y='ssy[um]', clim=None, 
              fill_events=False, cmap='Oranges', shift_zeros=False):
    """Use imshow to plot all of the tiles of a 2d flyscan

    Parameters
    ----------
    scan : DataFrame
        Scan dataframe
    x : str
        The data key that corresponds to the x axis
        Defaults to 'ssx[um]'
    y : str
        The data key that corresponds to the y axis
        Defaults to 'ssy[um]'
    clim : tuple, optional
        formtted as (min, max)
        If None, defaults to min/max of the data
    fill_events : bool, optional
        Fill the events with data from filestore
        Defaults to False (and is much much faster)
    cmap : str, optional
        Defaults to "Oranges"
        The colormap to use. See the pyplot.cm module for valid color maps
    shift_zeros : bool, optional
        Shift all zeros in the scan [NOTE: this is for the detector triggering-
        related bug only]
    """
    df = _load_scan(scan_id, fill_events=False)

    rois = [ch for ch in df.keys() if ch.startswith('Ch')]
    spectrum = np.sum([getattr(df, roi) for roi in rois], axis=0)
    x_data = df[x]
    y_data = df[y]
   
    if shift_zeros:
        spectrum = _shift_zeros(df, spectrum)
    
    fig = plt.figure()
    ax = plt.subplot(111)

    if clim is None:
        clim = (np.nanmin(spectrum), np.nanmax(spectrum))

    npts_total = 0

    points, log = _fly2d_points_from_olog(scan_id)
    for point_info in points:
        xi, xj, pointsx, yi, yj, pointsy = point_info
        pointsx += 1
        pointsy += 1

        npts = pointsx * pointsy
        slice_ = slice(npts_total, npts_total + npts)
        x = list(x_data[slice_])
        y = list(y_data[slice_])
        spec = spectrum[slice_]
        
        x0 = min((x[0], x[-1]))
        x1 = max((x[0], x[-1]))
        y0 = min((y[0], y[-1]))
        y1 = max((y[0], y[-1]))
        extent = (x0, x1, y1, y0)

        im = ax.imshow(spec.reshape((pointsy, pointsx)),
                       extent=extent, interpolation='None', cmap=cmap, 
                       vmin=clim[0], vmax=clim[1])

        npts_total += npts
    
    ax.set_xlim(np.min(x_data), np.max(x_data))
    ax.set_ylim(np.min(y_data), np.max(y_data))
    fig.colorbar(im)
