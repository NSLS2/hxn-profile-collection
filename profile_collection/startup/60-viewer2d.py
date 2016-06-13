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
    #data = np.sum(df['Det%d_%s' % (chan, elem)]
    #              for chan in channels)
    data = df[elem]
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

def plot(scan_id, namex, elem='Pt', channels=None, norm=None):
    plt.figure()
    plt.title(elem)

    scan_id, df = _load_scan(scan_id, fill_events=False)
    x = df[namex]

    if channels is 'sum':
        channels = [1, 2, 3]

    #scan_id, df = _load_scan(scan_id, fill_events=False)
    #x = df[namex]
        data = np.sum(df['Det%d_%s' % (chan, elem)]
                      for chan in channels)
    else:
        data = df[elem]

    x = np.asarray(x)
    data = np.asarray(data)

    if norm is not None:
        norm_v = df[norm]
        plt.plot(x, data / (norm_v + 1.e-8))
        plt.plot(x, data / (norm_v + 1.e-8), 'bo')
    else:
        plt.plot(x, data)
        plt.plot(x, data, 'bo')
        try:
            diff = np.diff(data)
            plt.figure()
            plt.plot(x[:-1], diff)
            plt.plot(x[:-1], diff, 'bo')
            plt.title('derivative')
        except Exception as ex:
            print('Failed to plot derivative: ({}) {}'
                  ''.format(ex.__class__.__name__, ex))
            raise

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


def plotfly(scan_id, elem='Pt', channels=None):
    if channels is None:
        channels = [1, 2, 3]

    plt.figure()

    scan_id, df = _load_scan(scan_id, fill_events=False)
    hdr = db[scan_id]['start']
    namex = hdr['fast_axis']

    x = df[namex]
    roi_data = np.sum(df['Det%d_%s' % (chan, elem)]
                      for chan in channels)
    try:
        diff = np.diff(roi_data)
        plt.subplot(122)
        plt.plot(x[1:], diff)
        plt.plot(x[1:], diff, 'bo')
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
    plt.title('Scan %d: %s' % (scan_id, elem))
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
def plot2dfly(scan_id, elem='Pt', *, x='ssx[um]', y='ssy[um]', clim=None,
              fill_events=False, cmap='Oranges', cols=None,
              channels=None, interp=None, interp2d=None):
    """Plot the results of a 2d fly scan

    Parameters
    ----------
    scan_id : int
        Any valid input to databroker[] or StepScan
    elem : str
        The element to display
        Defaults to 'Pt'
    x : str, optional
        The data key that corresponds to the x axis
        Defaults to 'ssx[um]'
    y : str, optional
        The data key that corresponds to the y axis
        Defaults to 'ssy[um]'
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
        channels = range(1, 4)

    scan_id, df = _load_scan(scan_id, fill_events=fill_events)

    title = 'Scan id %s. ' % scan_id + elem
    if elem in df:
        spectrum = np.asarray(df[elem])
    else:
        roi_keys = ['Det%d_%s' % (chan, elem) for chan in channels]

        for key in roi_keys:
            if key not in df:
                raise KeyError('ROI %s not found' % (key, ))

        spectrum = np.sum([getattr(df, roi) for roi in roi_keys], axis=0)
    x_data = np.asarray(df[x])
    y_data = np.asarray(df[y])

    hdr = db[scan_id]['start']

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

    if extent is not None:
        # create the scatter plot version
        scatter = ax2.scatter(x_data, y_data, c=spectrum, marker='s', s=250,
                              cmap=getattr(mpl.cm, cmap), linewidths=0,
                              alpha=.8, vmin=clim[0], vmax=clim[1])
        ax2.set_xlim(np.min(x_data), np.max(x_data))
        ax2.set_ylim(np.min(y_data), np.max(y_data))
        ax2.set_title('SCATTER. ' + title)
        ax2.set_aspect('equal')
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


def export(sid,num=1):
    for i in range(num):
        sid, df = _load_scan(sid, fill_events=False)
        path = os.path.join('/data/output/txt/', 'scan_{}.txt'.format(sid))
        print('Scan {}. Saving to {}'.format(sid, path))
        non_objects = [name for name, col in df.iteritems()
                       if col.dtype.name not in ('object', )]
        df.to_csv(path, float_format='%1.5e', sep='\t',
                  columns=sorted(non_objects))
        sid = sid + 1
    # return df


def get_all_filenames(scan_id, key='merlin1'):
    scan_id, df = _load_scan(scan_id, fill_events=False)
    from filestore.path_only_handlers import (AreaDetectorTiffPathOnlyHandler,
                                              RawHandler)
    handlers = {'AD_TIFF': AreaDetectorTiffPathOnlyHandler,
                'XSP3': RawHandler}
    filenames = [filestore.api.retrieve(uid, handlers)[0]
                 for uid in list(df[key])]

    if len(set(filenames)) != len(filenames):
        return set(filenames)
    return filenames
