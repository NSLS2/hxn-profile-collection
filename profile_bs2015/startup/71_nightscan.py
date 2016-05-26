import numpy as np
import datetime
from dataportal import DataBroker, DataMuxer


def _load_scan(scan_id, fill_events=False):
    '''Load scan from databroker by scan id'''

    if scan_id > 0 and scan_id in data_cache:
        df = data_cache[scan_id]
    else:
        hdr = DataBroker[scan_id]
        scan_id = hdr['start'].scan_id
        if scan_id in data_cache:
            df = data_cache[scan_id]
        else:
            data = DataBroker.fetch_events(hdr, fill=fill_events)
            dm = DataMuxer.from_events(data)
            df = dm.to_sparse_dataframe()
            data_cache[scan_id] = df
    return scan_id, df

def find_mass_center(array):
    n = np.size(array)
    tmp = 0
    sumtmp = 0
    for i in range(n):
        if array[i] > 50:
            tmp += i * array[i]
            sumtmp += array[i]
    if sumtmp > 0:
        mc = np.round(tmp/sumtmp)
    else:
        mc = np.round(n/2)
    return mc


def mov_to_center(scan_id, elem='Ga', channels=None, norm=None):
    if channels is None:
        channels = [1, 2, 3]

    scan_id, df = _load_scan(scan_id, fill_events=False)
    hdr = DataBroker[scan_id]['start']
    namex = hdr['fast_axis']

    x = df[namex]
    roi_data = np.sum(df['Det%d_%s' % (chan, elem)]
                      for chan in channels)
    x = np.asarray(x)
    roi_data = np.asarray(roi_data)
    mc = np.int(find_mass_center(roi_data))
    mov(zpssx,x[mc])


def night_scan(angle_list,x_range,y_range,step_size,exposure_time):
    angle_list = np.array(angle_list)
    num_angle = np.size(angle_list)
    now=datetime.datetime.now()
    hh = open('itiff_prefixes_nightscan_'+str(now.isoformat())+'.txt','w')

    for i in range(num_angle):
        angle = angle_list[i]
        mov(zpsth,angle)
        hh.write(str(angle))
        hh.write('\n')
        sleep(2)

        fly1d(zpssx, -2, 2, 100, 0.5)
        mov_to_center(-1,'Ga')
        sleep(1)

        fly2d(zpssx, -.75*y_range/2, .75*y_range/2, 30, zpssy, -y_range/2, y_range/2, 30, 0.5, return_speed=40)
        export(-1)
        hh.write(caget("XF:03IDC-ES{Merlin:1}TIFF1:FileName_RBV",as_string=True))
        hh.write('\n')
        fermat(zpssx,zpssy, x_range, y_range, step_size, 1, exposure_time)
        hh.write(caget("XF:03IDC-ES{Merlin:1}TIFF1:FileName_RBV",as_string=True))
        hh.write('\n')
        export(-1)
        scatter_plot(-1, 'zpssx', 'zpssy', 'Det2_Ga')

    hh.close()

