import os
import sys
import numpy as np
import filestore
import filestore.api
from datetime import datetime

import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

import databroker
from databroker import (get_table, db)
# from xray_vision.qt_widgets import CrossSectionMainWindow
# from xray_vision.backend.mpl.cross_section_2d import CrossSection
from scipy.interpolate import interp1d, interp2d


def plot2d(scan_id, elem, norm='sclr1_ch4', det_type='elem'):
    scan_id, df = _load_scan(scan_id, fill_events=False)
    scan_info=db[scan_id]
    tmp = scan_info['start']
    #tmp = tmp.split()
    x_motor = tmp['motors'][0]
    y_motor = tmp['motors'][1]

    x_start = tmp['plan_args']['args'][1]
    x_end = tmp['plan_args']['args'][2]
    col = tmp['plan_args']['args'][3]

    y_start = tmp['plan_args']['args'][5]
    y_end = tmp['plan_args']['args'][6]
    row = tmp['plan_args']['args'][7]

    if det_type == 'elem':
        det = (df['Det1_{}'.format(elem)] +
               df['Det2_{}'.format(elem)] +
               df['Det3_{}'.format(elem)])
    elif det_type == 'scalar':
        det = df[elem]
    else:
        det = df[elem]

    if norm is not None:
        mon = np.reshape(df[norm], (row,col))
        plt.figure()
        data = np.reshape(det, (row, col))
        plt.title('Scan %d: %s (normalized to %s)' % (scan_id, elem, norm))
        plt.imshow(data/mon, interpolation='None',extent=[x_start,x_end,y_end,y_start],cmap='bone')
        plt.xlabel(x_motor)
        plt.ylabel(y_motor)
        plt.colorbar()
    else:
        plt.figure()
        data = np.reshape(det, (row, col))
        plt.title('Scan %d: %s' % (scan_id, elem))
        plt.imshow(data, interpolation='None',extent=[x_start,x_end,y_end,y_start])
        plt.xlabel(x_motor)
        plt.ylabel(y_motor)
        plt.colorbar()


def dev(scan_id, namex, namey):
    d = 3.13559

    scan_id, df = _load_scan(scan_id, fill_events=False)
    dety = df[namey]

    if namex == "energy":
        detx = df["dcm_th"]
        num_points = len(detx)
        data = np.zeros((num_points - 1, 2))
        for i in range(num_points - 1):
            data[i, 1] = (dety[i + 1] - dety[i]) / (detx[i + 1] - detx[i])
            tmp = (detx[i + 1] + detx[i]) / 2
            s = np.sin(np.pi * (tmp + (-0.0135)) / 180)
            data[i, 0] = 12.398 / (2 * d * s)
    else:
        detx = df[namex]
        num_points = len(detx)
        data = np.zeros((num_points - 1, 2))
        for i in range(num_points - 1):
            data[i, 1] = (dety[i + 1] - dety[i]) / (detx[i + 1] - detx[i])
            data[i, 0] = (detx[i + 1] + detx[i]) / 2

    plt.figure(20)
    plt.plot(data[:, 0], data[:, 1])
    # return data


def scatter_plot(scan_id, namex,namey, elem='Pt', channels=None, norm=None):
    plt.figure()
    plt.title(elem)
    if channels is None:
        channels = [1, 2, 3]
    scan_id, df = _load_scan(scan_id, fill_events=False)
    x = df[namex]
    y = df[namey]
    data = np.sum(df['Det%d_%s' % (chan, elem)]
                  for chan in channels)
    #data = df[elem]
    x = np.asarray(x)
    y = np.asarray(y)
    data = np.asarray(data)
    if norm is not None:
        norm_v = df[norm]
        plt.scatter(x,y, c=data / (norm_v + 1.e-8),s=200)
        plt.gca().invert_yaxis()
        plt.axes().set_aspect('equal','datalim')
        plt.xlabel(namex)
        plt.ylabel(namey)
    else:
        plt.scatter(x,y, c=data,s=200)
        plt.gca().invert_yaxis()
        plt.axes().set_aspect('equal','datalim')
        plt.xlabel(namex)
        plt.ylabel(namey)
    plt.show()

def plot(scan_id, elem='Pt', norm=None,center_method='com',log=0,e_flag=0):
    plt.figure()
    scan_id, df = _load_scan(scan_id, fill_events=False)
    hdr = db[scan_id]['start']
    scan_start_time = datetime.isoformat(datetime.fromtimestamp(hdr['time']))

    if elem in df:
        data = np.asarray(df[elem])
    else:
        channels = [1, 2, 3]
        roi_keys = ['Det%d_%s' % (chan, elem) for chan in channels]
        for key in roi_keys:
            if key not in df:
                raise KeyError('ROI %s not found' % (key, ))
        data = np.sum([getattr(df, roi) for roi in roi_keys], axis=0)

    scanned_axis = hdr['motors'][0]

    if scanned_axis == 'ugap':
        scanned_axis = 'ugap_readback'
    x = df[scanned_axis]

    if e_flag:
        x = 12.39842 / (2.*3.1355893*np.sin(np.deg2rad(x)))
    '''
    if channels is 'sum':
        channels = [1, 2, 3]
        data = np.sum(df['Det%d_%s' % (chan, elem)]
                      for chan in channels)
    else:
        data = df[elem]
    '''
    x = np.asarray(x)
    data = np.asarray(data)

    if norm is not None:
        norm_v = df[norm]
        if log:
            plt.plot(x,np.log10(data / (norm_v+1.e-8)))
            plt.plot(x,np.log10(data / (norm_v + 1.e-8)),'bo')
        else:
            plt.plot(x, data / (norm_v + 1.e-8))
            plt.plot(x, data / (norm_v + 1.e-8), 'bo')
        if e_flag:
            plt.xlabel('Energy (keV)')
        else:
            plt.xlabel(scanned_axis)
        plt.ylabel(elem)
        plt.title('Scan %d' % (scan_id))
    else:
        if log:
            plt.plot(x,np.log10(data+1.e-8))
            plt.plot(x,np.log10(data+1.e-8),'bo')
        else:
            plt.plot(x, data)
            plt.plot(x, data, 'bo')
        if e_flag:
            plt.xlabel('Energy (keV)')
        else:
            plt.xlabel(scanned_axis)
        plt.ylabel(elem)
        plt.title('Scan %d' % (scan_id))
        try:
            diff = np.diff(data)
            plt.figure()
            plt.plot(x[:-1], diff)
            plt.plot(x[:-1], diff, 'bo')
            mc = find_mass_center(data)
            #plt.title('center of mass: %d',x[mc])
        except Exception as ex:
            print('Failed to plot derivative: ({}) {}'
                  ''.format(ex.__class__.__name__, ex))
            raise

    plt.title('Scan %d: %s    Start time: %s' % (scan_id, elem, scan_start_time))
    plt.show()


def plot_all(scan_id, namex=None, diff=False, channels=None,
             same_axis=False):
    plt.figure()

    if channels is None:
        channels = [1, 2, 3]

    scan_id, df = _load_scan(scan_id, fill_events=False)
    plt.title('Scan id: {}'.format(scan_id))

    x = df[namex]
    elems = set(key.split('_', 1)[1] for key in df
                if key.startswith('Det'))

    if same_axis:
        ax = plt.subplot(111)
    else:
        n_elem = len(elems)
        cols = rows = int(np.ceil(np.sqrt(n_elem)))
        gs = gridspec.GridSpec(rows, rows)

    print('All elements:', list(elems))
    for i, elem in sorted(enumerate(elems)):
        if not same_axis:
            ax = plt.subplot(gs[i])
            ax.set_title(elem)

            # share the x-axes in columns
            if i < (n_elem - cols):
                plt.setp(ax.get_xticklabels(), visible=False)

        data = np.sum(df['Det%d_%s' % (chan, elem)]
                      for chan in channels)

        ax.plot(x, data, label=elem)
        ax.plot(x, data, 'bo')

    if same_axis:
        plt.legend(loc='best')

    plt.show()


def find_mass_center(array):
    n = np.size(array)
    tmp = 0
    for i in range(n):
        tmp += i * array[i]
    mc = np.round(tmp / np.sum(array))
    return mc

def plotfly(scan_id, elem='Pt', norm=None,center_method='com'):
    plt.figure()
    scan_id, df = _load_scan(scan_id, fill_events=False)
    hdr = db[scan_id]['start']
    scan_start_time = datetime.isoformat(datetime.fromtimestamp(hdr['time']))
    if elem in df:
        roi_data = np.asarray(df[elem])
    else:
        channels = [1, 2, 3]
        roi_keys = ['Det%d_%s' % (chan, elem) for chan in channels]
        for key in roi_keys:
            if key not in df:
                raise KeyError('ROI %s not found' % (key, ))
        roi_data = np.sum([getattr(df, roi) for roi in roi_keys], axis=0)

    scanned_axis = hdr['motor']
    x = df[scanned_axis]

    if norm is not None:
        norm_tot = df[norm]
        roi_data = roi_data/(norm_tot + 1e-8)

    try:
        diff = np.diff(roi_data)
        plt.subplot(122)
        plt.plot(x[1:], diff)
        plt.plot(x[1:], diff, 'bo')
        #if center_method == 'com':
        #    i_center = find_mass_center(roi_data)
        #else:
        i_max = np.where(diff == np.max(diff))
        i_min = np.where(diff == np.min(diff))
        i_center = np.round((i_max[0][0]+i_min[0][0])/2)+1
        plt.title('Scan %d: %s (deriv)' % (scan_id, elem)+' Center: '+np.str(x[i_center]))
    except Exception as ex:
        print('Failed to plot derivative: ({}) {}'
              ''.format(ex.__class__.__name__, ex))
        plt.clf()
        plt.subplot(111)
    else:
        plt.subplot(121)

    plt.plot(x, roi_data)
    plt.plot(x, roi_data, 'bo')
    plt.xlabel(scanned_axis)
    plt.ylabel(elem)
    plt.title('Scan %d: %s    Start time: %s' % (scan_id, elem, scan_start_time))
    plt.show()


if 'data_cache' not in globals():
    # Don't erase the cache when reloading this module via %run -i
    data_cache = {}


def _load_scan(scan_id, fill_events=False):
    '''Load scan from databroker by scan id'''

    if scan_id > 0 and scan_id in data_cache:
        df = data_cache[scan_id]
    else:
        hdr = db[scan_id]
        scan_id = hdr['start'].scan_id
        if scan_id not in data_cache:
            data_cache[scan_id] = db.get_table(hdr, fill=fill_events)

        df = data_cache[scan_id]

    return scan_id, df


def get_flyscan_dimensions(hdr):
    if 'dimensions' in hdr:
        return hdr['dimensions']
    else:
        return hdr['shape']


def fly2d_grid(hdr, x_data=None, y_data=None, plot=False):
    '''Get ideal gridded points for a 2D flyscan'''
    try:
        nx, ny = get_flyscan_dimensions(hdr)
    except ValueError:
        raise ValueError('Not a 2D flyscan')

    rangex, rangey = hdr['scan_range']
    width = rangex[1] - rangex[0]
    height = rangey[1] - rangey[0]

    if 'scan_starts' in hdr:
        start_x, start_y = hdr['scan_starts'][0]
    else:
        macros = eval(hdr['subscan_0']['macros'], dict(array=np.array))
        start_x, start_y = macros['scan_starts']

    dx = width / nx
    dy = height / ny
    grid_x = np.linspace(start_x, start_x + width + dx / 2, nx)
    grid_y = np.linspace(start_y, start_y + height + dy / 2, ny)

    if plot:
        mesh_x, mesh_y = np.meshgrid(grid_x, grid_y)
        plt.figure()
        if x_data is not None and y_data is not None:
            plt.scatter(x_data, y_data, c='blue', label='actual')
        plt.scatter(mesh_x, mesh_y, c='red', label='gridded',
                    alpha=0.5)
        plt.legend()
        plt.show()

    return grid_x, grid_y


def interp2d_scan(hdr, x_data, y_data, spectrum, *, kind='linear',
                  plot_points=False, **kwargs):
    '''Interpolate a 2D flyscan over a grid'''

    new_x, new_y = fly2d_grid(hdr, x_data, y_data, plot=plot_points)

    f = interp2d(x_data, y_data, spectrum, kind=kind, **kwargs)
    return f(new_x, new_y)


def interp1d_scan(hdr, x_data, y_data, spectrum, kind='linear',
                  plot_points=False, **kwargs):
    '''Interpolate a 2D flyscan only over the fast-scanning direction'''
    grid_x, grid_y = fly2d_grid(hdr, x_data, y_data, plot=plot_points)
    x_data = fly2d_reshape(hdr, x_data, verbose=False)

    spectrum2 = np.zeros_like(spectrum)
    for row in range(len(grid_y)):
        spectrum2[row, :] = interp1d(x_data[row, :], spectrum[row, :],
                                     kind=kind, bounds_error=False,
                                     **kwargs)(grid_x)

    return spectrum2


def fly2d_reshape(hdr, spectrum, verbose=True):
    '''Reshape a 1D array to match the shape of a 2D flyscan'''
    try:
        nx, ny = get_flyscan_dimensions(hdr)
    except ValueError:
        raise ValueError('Not a 2D flyscan')

    try:
        spectrum2 = spectrum.copy().reshape((ny, nx))
    except Exception as ex:
        if verbose:
            print('\tUnable to reshape data to (%d, %d) (%s: %s)'
                  '' % (nx, ny, ex.__class__.__name__, ex))
    else:
        fly_type = hdr['fly_type']
        if fly_type in ('pyramid', ):
            # Pyramid scans' odd rows are flipped:
            if verbose:
                print('\tPyramid scan. Flipping odd rows.')
            spectrum2[1::2, :] = spectrum2[1::2, ::-1]

        return spectrum2


# TODO: change l, h to clim which defaults to 'auto'
def plot2dfly(scan_id, elem='Pt', norm=None, *, x=None, y=None, clim=None,
              fill_events=False, cmap='jet', cols=None,
              channels=None, interp=None, interp2d=None):
    """Plot the results of a 2d fly scan

    Parameters
    ----------
    scan_id : int
        Any valid input to databroker[] or StepScan
    elem : str
        The element to display
        Defaults to 'Pt'
    norm : str, optional
        scaler for intensity normalization
    x : str, optional
        The data key that corresponds to the x axis
    y : str, optional
        The data key that corresponds to the y axis
    clim : tuple, optional
        formatted as (min, max)
        If None, defaults to min/max of the data
    fill_events : bool, optional
        Fill the events with data from filestore
        Defaults to False (and is much much faster)
    cmap : str, optional
        Defaults to "Oranges"
        The colormap to use. See the pyplot.cm module for valid color maps
    channels : list, optional
        The channels to use (defaults to 1 to 3)
    interp : {'linear', 'cubic', 'quintic'}, optional
        Interpolate the data on the 2D mesh defined by positioners x and y,
        only in the x direction
    interp2d : {'linear', 'cubic', 'quintic'}, optional
        Interpolate the data on the 2D mesh defined by positioners x and y,
        in both the x and y directions (NOTE: _extremely_ slow)
    """

    if channels is None:
        channels = [1, 2, 3]

    scan_id, df = _load_scan(scan_id, fill_events=fill_events)

    title = 'Scan id %s. ' % scan_id + elem
    if elem in df:
        spectrum = np.asarray(df[elem], dtype=np.float32)
    else:
        roi_keys = ['Det%d_%s' % (chan, elem) for chan in channels]

        for key in roi_keys:
            if key not in df:
                raise KeyError('ROI %s not found' % (key, ))

        spectrum = np.sum([getattr(df, roi) for roi in roi_keys], axis=0)

    hdr = db[scan_id]['start']
    if x is None:
        x = hdr['motor1']
    x_data = np.asarray(df[x])

    if y is None:
        y = hdr['motor2']
    y_data = np.asarray(df[y])

    if norm is not None:
        monitor = np.asarray(df[norm],dtype=np.float32)
        spectrum = spectrum/(monitor + 1e-8)


    nx, ny = get_flyscan_dimensions(hdr)
    total_points = nx * ny

    if clim is None:
        clim = (np.nanmin(spectrum), np.nanmax(spectrum))
    extent = (np.nanmin(x_data), np.nanmax(x_data),
              np.nanmax(y_data), np.nanmin(y_data))

    # these values are also used to set the limits on the value
    if ((abs(extent[0] - extent[1]) <= 0.001) or
            (abs(extent[2] - extent[3]) <= 0.001)):
        extent = None

    dt = datetime.utcnow()
    folder = os.path.join('/data/output/',
                          '{}{:0>2}{:0>2}/'.format(dt.year, dt.month, dt.day))

    if not os.path.exists(folder):
        os.makedirs(folder)

    print('Scan {}. Saving to: {}'.format(scan_id, folder))

    if len(spectrum) != total_points:
        print('Padding data (points=%d expected=%d)' % (len(spectrum),
                                                        total_points))

        _spectrum = np.zeros(total_points, dtype=spectrum.dtype)
        _spectrum[:len(spectrum)] = spectrum
        spectrum = _spectrum

    if interp2d is not None:
        print('\tUsing 2D %s interpolation...' % (interp2d, ), end=' ')
        sys.stdout.flush()
        spectrum = interp2d_scan(hdr, x_data, y_data, spectrum,
                                 kind=interp2d)
        print('done')

    spectrum2 = fly2d_reshape(hdr, spectrum)

    if interp is not None:
        print('\tUsing 1D %s interpolation...' % (interp, ), end=' ')
        sys.stdout.flush()
        spectrum2 = interp1d_scan(hdr, x_data, y_data, spectrum2, kind=interp)
        print('done')

    fig = None
    ax1 = None
    ax2 = None

    if spectrum2 is None:
        fig = plt.figure()
        ax2 = plt.subplot(111)
    else:
        fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10, 5))
        fig.set_tight_layout(True)
        ax1.imshow(spectrum2, extent=extent, interpolation='None', cmap=cmap,
                   vmin=clim[0], vmax=clim[1])
        np.savetxt(os.path.join(folder, 'data_scan_{}'.format(scan_id)),
                   spectrum2)

        ax1.set_title('IMSHOW. ' + title)
        ax1.set_xlabel(x)
        ax1.set_ylabel(y)

    if extent is not None:
        # create the scatter plot version
        scatter = ax2.scatter(x_data, y_data, c=spectrum, marker='s', s=250,
                              cmap=getattr(mpl.cm, cmap), linewidths=0,
                              alpha=.8, vmin=clim[0], vmax=clim[1])
        ax2.set_xlabel(x)
        ax2.set_xlim(np.min(x_data), np.max(x_data))
        ax2.set_ylabel(y)
        ax2.set_ylim(np.min(y_data), np.max(y_data))
        ax2.set_title('SCATTER. ' + title)
        ax2.set_aspect('equal')
        ax2.invert_yaxis()
        fig.colorbar(scatter)


    fig_path = os.path.join(folder,'data_scan_{}.png'.format(scan_id))
    print('\tSaving figure to: {}'.format(fig_path))
    fig.savefig(fig_path)

    text_path = os.path.join(folder,'data_x_y_ch_{}'.format(scan_id))
    print('\tSaving text positions to: {}'.format(text_path))
    np.savetxt(text_path, np.vstack((x_data, y_data, spectrum)).T)

    var_name = 'S_%d_%s' % (scan_id, elem)
    globals()[var_name] = spectrum2
    print('\tScan data available in variable: {}'.format(var_name))
    return fig, ax1, ax2


def export(sid,num=1, export_folder='/home/xf03id/data_analysis/',
           fields_excluded=['xspress3_ch1', 'xspress3_ch2', 'xspress3_ch3', 'merlin1']):
    for i in range(num):
        sid, df = _load_scan(sid, fill_events=False)
        path = os.path.join(export_folder, 'scan_{}.txt'.format(sid))
        print('Scan {}. Saving to {}'.format(sid, path))
        #non_objects = [name for name, col in df.iteritems()
        #               if col.dtype.name not in ('object', )]
        non_objects = [name for name in df.keys() if name not in fields_excluded]
        print('fields inclued: {}'.format(sorted(non_objects)))
        #dump all data
        #non_objects = [name for name, col in df.iteritems()]
        df.to_csv(path, float_format='%1.5e', sep='\t',
                  columns=sorted(non_objects))

        #path = os.path.join('/home/xf03id/data_analysis/Amy_Aug2016/', 'scan_{}_raw.txt'.format(sid))
        #np.savetxt(path, (df['sclr1_ch4'], df['zpssx'], df['zpssy']), fmt='%1.5e')

        sid = sid + 1
    # return df


def get_all_filenames(scan_id, key='merlin1'):
    scan_id, df = _load_scan(scan_id, fill_events=False)
    from filestore.path_only_handlers import (AreaDetectorTiffPathOnlyHandler,
                                              RawHandler)
    handlers = {'AD_TIFF': AreaDetectorTiffPathOnlyHandler,
                'XSP3': RawHandler,
                'AD_HDF5': RawHandler,
                'TPX_HDF5': RawHandler,
                }
    filenames = [filestore.api.retrieve(uid, handlers)[0]
                 for uid in list(df[key])]

    if len(set(filenames)) != len(filenames):
        return set(filenames)
    return filenames
