from __future__ import print_function
import matplotlib.pyplot as plt
import numpy as np
from xray_vision.qt_widgets import CrossSectionMainWindow
from xray_vision.backend.mpl.cross_section_2d import CrossSection
from datetime import datetime
import os


def plot2d(scan, name,row,col):

    det = getattr(dscan.data[scan], 'Det1_'+name)+getattr(dscan.data[scan], 'Det2_'+name)+getattr(dscan.data[scan], 'Det3_'+name)

    plt.figure()
    data = np.reshape(det, (row, col))
    plt.imshow(data,interpolation='None')
    plt.colorbar()
    '''
    num_det = 1;
    if name == "all":
        for i in range(num_det):
            #det_name = dscan.detectors[i].name
            det = getattr(dscan.data[scan], 'Fe_ch1')+getattr(dscan.data[scan], 'Fe_ch2')+getattr(dscan.data[scan], 'Fe_ch3')
            plt.figure(i)
            data = np.reshape(det, (row, col))
            plt.imshow(data,interpolation='None')
    else:
        plt.figure()
        det = getattr(dscan.data[scan], name)
        data = np.reshape(det, (row, col))
        plt.imshow(data,interpolation='None')
    # return data
    '''

def dev(scan, namex, namey):
    dety = getattr(dscan.data[scan], namey)
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


def plot(scan, namex, elem='Pt', channels=None, norm='None'):
    plt.figure()
    plt.clf()
    if channels is None:
        channels = [1, 2, 3]

    data = np.sum(getattr(dscan.data[scan],'Det%d_%s' % (chan, elem))
                for chan in channels)

    # data = getattr(dscan.data[scan], 'Det1_'+namey)+getattr(dscan.data[scan], 'Det2_'+namey)+getattr(dscan.data[scan], 'Det3_'+namey)

    x = getattr(dscan.data[scan], namex)

    if norm != 'None':
        norm_v = getattr(dscan.data[-1], norm)
        plt.plot(x,data/(norm_v+1.e-8))
        plt.plot(x,data/(norm_v+1.e-8),'bo')
    else:
        plt.plot(x,data)
        plt.plot(x,data,'bo')
        plt.figure()
        plt.plot(x[:-1],data[1:]-data[:-1])
        plt.plot(x[:-1],data[1:]-data[:-1],'bo')
        plt.title('derivative')
    plt.show()


def plotfly(namex, elem='Pt', channels=None):
    if channels is None:
        channels = [1, 2, 3]

    plt.figure()
    x = fly_data[namex]
    Pt = np.sum(fly_data['Det%d_%s' % (chan, elem)]
                for chan in channels)
    Pt[0] = Pt[1]
    Pt = np.array(Pt)

    plt.subplot(121)
    plt.plot(x,Pt)
    plt.plot(x,Pt,'bo')
    plt.title(elem)

    plt.subplot(122)
    plt.plot(x[1:],Pt[1:]-Pt[:-1])
    plt.plot(x[1:],Pt[1:]-Pt[:-1],'bo')
    plt.title('%s (deriv)' % elem)

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

    return scan_id, df


# TODO: change l, h to clim which defaults to 'auto'
def plot2dfly(scan_id, x='ssx[um]', y='ssy[um]', elem='Pt', clim=None,
              fill_events=False, cmap='Oranges', shift_zeros=False, cols=None,
              channels=None):
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
    elem : str
        The element to display
        Defaults to 'Pt'
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
    channels : list, optional
        The channels to use (defaults to 1 to 3)
    """

    if channels is None:
        channels = range(1, 4)

    scan_id, df = _load_scan(scan_id, fill_events=fill_events)

    title = 'Scan id %s. ' % scan_id + elem
    roi_keys = ['Det%d_%s' % (chan, elem) for chan in channels]

    for key in roi_keys:
        if key not in df:
            raise KeyError('ROI %s not found' % (key, ))

    spectrum = np.sum([getattr(df, roi) for roi in roi_keys], axis=0)
    spectrum[0] = spectrum[1]
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
        spectrum2 = spectrum.copy()
        spectrum2 = spectrum2.reshape((-1, cols))
    except Exception as ex:
        print('Unable to reshape data to width: %d (%s: %s)' % (cols, ex.__class__.__name__, ex))
        fig = plt.figure()
        ax2 = plt.subplot(111)
    else:
        print('Reshaped to %s' % (spectrum2.shape, ))
        print(np.shape(spectrum2))
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,5))
        fig.set_tight_layout(True)
  #      rows = np.size(spectrum2)/cols
  #      for row in range(1, rows, 2):
  #          spectrum2[row, :] = spectrum2[row, ::-1]

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
    globals()['S_%d_%s' % (scan_id, elem)] = spectrum2

def export(sid):
    data = StepScan[sid]
    np.savetxt('/data/txt/scan_'+np.str(sid)+'.txt',data,fmt='%1.5e')
