import matplotlib.pyplot as plt
import numpy as np


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
        plt.imshow(data)
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
    data1 = fly_data['Ch1 [930:960]']
    data2 = fly_data['Ch2 [930:960]']
    data3 = fly_data['Ch3 [930:960]']
    Pt = data1 + data2 + data3
    plt.plot(x,Pt)
    plt.plot(x,Pt,'bo')
    plt.show()

def plot2dfly(row, col,l,h,sid):

    import epics
    roi_pvs = ['XF:03IDC-ES{Xsp:1}:C1_ROI1:ArrayData_RBV',
               'XF:03IDC-ES{Xsp:1}:C2_ROI1:ArrayData_RBV',
               'XF:03IDC-ES{Xsp:1}:C3_ROI1:ArrayData_RBV']
    #roi_pvs = ['XF:03IDC-ES{Xsp:1}:C1_ROI2:ArrayData_RBV']
    counter_pv = 'XF:03IDC-ES{Xsp:1}:ArrayCounter_RBV'
    spectrum = np.sum([epics.caget(pv).astype(float)
                               for pv in roi_pvs], axis=0)
    data = spectrum[:row*col]
    data.resize(row,col)
    #sid = session.get_next_scan_id() - 1
    plt.figure()
    plt.imshow(data, interpolation='none')
    plt.colorbar()
    plt.clim([l,h])
    plt.title('Scan '+np.str(sid))
    plt.savefig('/data/20150425/data_scan_'+np.str(sid)+'.png')
    np.savetxt('/data/20150425/data_scan_'+np.str(sid),data)

    '''
    plt.figure()
    data1 = fly_data['Ch1 [930:960]'].reshape(row, col)
    data2 = fly_data['Ch2 [930:960]'].reshape(row, col)
    data3 = fly_data['Ch3 [930:960]'].reshape(row, col)
    data = data1 + data2 + data3
    plt.imshow(data, interpolation='none')
    plt.show()
    '''
