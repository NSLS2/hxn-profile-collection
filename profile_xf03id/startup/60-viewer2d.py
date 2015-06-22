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

def plot(namex,namey):
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

data_cache = {}

#todo change l, h to clim which defaults to 'auto'
def plot2dfly(scan_id, x='ssx[um]', y='ssy[um]', clim=None, fill_events=False, as_image=True, cmap='CMRmap'):
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
    """
    
    title = 'Scan id %s. ' % scan_id
    if scan_id > 0 and scan_id in data_cache:
        df = data_cache[scan_id]
    else:
        hdr = DataBroker[scan_id]
        if hdr.scan_id in data_cache:
            df = data_cache[scan_id]
        else:
            data = DataBroker.fetch_events(hdr, fill=fill_events)
            dm = DataMuxer.from_events(data)
            df = dm.to_sparse_dataframe()
            data_cache[scan_id] = df
    rois = [ch for ch in list(df) if ch.startswith('Ch')]
    spectrum = np.sum([getattr(df, roi) for roi in rois], axis=0)
    x_data = df[x]
    y_data = df[y]
    prev_x = x_data[0]
    for idx, pos in enumerate(x_data[1:]):
        if pos > prev_x:
            prev_x = pos
        else:
            break
    spectrum2 = spectrum.reshape((-1, idx+1))
    if clim is None:
        clim = (np.nanmin(spectrum2), np.nanmax(spectrum2))
    extent = (np.min(x_data), np.max(x_data), np.max(y_data), np.min(y_data))
    fig, (ax1, ax2) = plt.subplots(ncols=2, figsize=(10,5))
    fig.set_tight_layout(True)
    im = ax1.imshow(spectrum2, extent=extent, interpolation='None', cmap=cmap, 
                   vmin=clim[0], vmax=clim[1])
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
    # save it to disk
    dt = datetime.utcnow()
    folder = '/data/{}{:0>2}{:0>2}/'.format(dt.year, dt.month, dt.day)
    if not os.path.exists(folder):
        os.makedirs(folder)
    fig.savefig(folder + 'data_scan_'+np.str(scan_id)+'.png')
    np.savetxt(folder + 'data_scan_'+np.str(scan_id), spectrum2)
    np.savetxt(folder + 'data_x_y_ch_'+np.str(scan_id),
               np.vstack((x_data, y_data, spectrum)))
    
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
