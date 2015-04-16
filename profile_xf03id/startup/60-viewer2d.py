import matplotlib.pyplot as plt
import numpy as np
def show2d(scan, name,row,col):
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
    #reture data


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
    plt.figure(1)
    plt.clf()
    data = getattr(dscan.data[-1], namey)
    x =  getattr(dscan.data[-1], namex)
    plt.plot(x,data)
    plt.show()

    
