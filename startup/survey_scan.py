import numpy as np
import matplotlib.pyplot as plt
import sys
from scipy import ndimage
#from databroker import get_table, db
#from skimage.filters.rank import median
from skimage.morphology import disk
from skimage import io

'''
from hxntools.handlers import register
import filestore
register()

import yaml
with open("./hxn.yml","r") as read_file:
    data=yaml.load(read_file)

from metadatastore.mds import MDS
from databroker import Broker
from databroker.core import register_builtin_handlers
from filestore.fs import FileStore

_mds_config = {'host': data['metadatastore']['config']['host'],
               'port': 27017,
               'database': data['metadatastore']['config']['database'],
               'timezone': 'US/Eastern'}
mds = MDS(_mds_config)
_fs_config = {'host': data['assets']['config']['host'],
              'port': 27017,
              'database': data['assets']['config']['database']}
db = Broker(mds, FileStore(_fs_config))

from hxntools.handlers.timepix import TimepixHDF5Handler
from hxntools.handlers.xspress3 import Xspress3HDF5Handler
db.fs.register_handler(TimepixHDF5Handler._handler_name, TimepixHDF5Handler, overwrite=True)
'''

def rm_pixel(data,ix,iy):
    data[ix,iy] = np.median(data[ix-1:ix+1,iy-1:iy+1])
    return data

def display_frame(index):
    #t = np.flipud(np.squeeze(np.asarray(filestore.api.get_data(df['merlin1'][np.ceil(index)]))))
    #print(index)
    t = diff_array[:,:,index]
    fig3 = plt.figure(3)
    plt.clf()
    cl = np.percentile(diff_array,15)
    ch = np.percentile(diff_array,99.95)
    #im2 = plt.imshow(np.flipud(np.log10(t+.001)), cmap = 'spectral', interpolation = 'none')
    im2 = plt.imshow(np.flipud(t.T), cmap = 'hot', interpolation = 'none',clim=[cl,ch])
    #im2 = plt.imshow(np.flipud(np.log10(t+0.001).T), cmap = 'hot', interpolation = 'none')
    iy = np.floor(index/num_x)
    ix = np.mod(index, num_x)
    index_new = (num_y-iy) * num_x + (num_x - ix)
    #plt.title('x: '+np.str(x_data[index])+'um, y: '+np.str(y_data[index])+'um')
    plt.title('x: %.3f' % x_data[index]+' um, y: %.3f' %y_data[index]+' um')
    plt.colorbar()
    plt.draw()

    '''
    t2 = diff_array_2[:,:,index]
    fig4 = plt.figure(4)
    plt.clf()
    #im2 = plt.imshow(np.flipud(np.log10(t+.001)), cmap = 'spectral', interpolation = 'none')
    im2 = plt.imshow(np.flipud(np.log10(t2+0.001).T), cmap = 'hot', interpolation = 'none')#,clim=[cl2,ch2])
    iy = np.floor(index/num_x)
    ix = np.mod(index, num_x)
    index_new = (num_y-iy) * num_x + (num_x - ix)
    #plt.title('x: '+np.str(x_data[index])+'um, y: '+np.str(y_data[index])+'um')
    plt.title('x: %.3f' % x_data[index]+' um, y: %.3f' %y_data[index]+' um')
    plt.colorbar()
    plt.draw()
    '''

def onclick(event):
    #print(w,l,event.xdata,event.ydata)
    index = w*(np.int(np.round(event.ydata)))  + np.int(np.round(event.xdata))
    fig = plt.figure(0)
    plt.clf()
    plt.imshow(xrf,interpolation = 'none',aspect='auto')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' '+elem)
    plt.draw()

    fig1 = plt.figure(1)
    plt.clf()
    plt.imshow(roi,interpolation = 'none',aspect='auto')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' Diff ROI')
    plt.draw()

    '''
    fig2 = plt.figure(2)
    plt.clf()
    im1 = plt.imshow(roi2,interpolation = 'none')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' Trans ROI')
    plt.draw()
    '''
    #print(index)
    display_frame(index)


def onclick_roi(event):
    #print(w,l,event.xdata,event.ydata)
    index = w*(np.int(np.round(event.ydata)))  + np.int(np.round(event.xdata))
    fig1 = plt.figure(1)
    plt.clf()
    plt.imshow(roi,interpolation = 'none',aspect='auto')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' Diff ROI')
    plt.draw()


    fig1 = plt.figure(0)
    plt.clf()
    plt.imshow(xrf,interpolation = 'none',aspect='auto')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' '+elem)
    plt.draw()

    '''
    fig1 = plt.figure(2)
    plt.clf()
    im1 = plt.imshow(roi2,interpolation = 'none')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' Trans ROI')
    plt.draw()
    '''
    #print(index)
    display_frame(index)

def onclick_roi2(event):
    #print(w,l,event.xdata,event.ydata)
    index = w*(np.int(np.round(event.ydata)))  + np.int(np.round(event.xdata))
    fig1 = plt.figure(2)
    plt.clf()
    im1 = plt.imshow(roi2,interpolation = 'none')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' tran ROI')
    plt.draw()


    fig1 = plt.figure(0)
    plt.clf()
    im1 = plt.imshow(xrf,interpolation = 'none')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' '+elem)
    plt.draw()


    fig1 = plt.figure(1)
    plt.clf()
    im1 = plt.imshow(roi,interpolation = 'none')
    plt.scatter(np.int(np.round(event.xdata)),np.int(np.round(event.ydata)),zorder=1)
    plt.title('#'+scan_num +' Diff ROI')
    plt.draw()

    #print(index)
    display_frame(index)


def onclick_fermat(event):
    index = (np.abs(x-event.xdata)**2+np.abs(y-event.ydata)**2).argmin()+1
    fig1 = plt.figure(1)
    plt.clf()
    props = dict(alpha=0.8, edgecolors='none' )
    plt.scatter(x,y,c=xrf,s=50,marker='s',**props)
    plt.axes().set_aspect('equal')
    plt.gca().invert_yaxis()
    plt.scatter(event.xdata,event.ydata,zorder=1)
    plt.draw()
    display_frame(index)


def show_diff_data(sid,element,det_name='merlin1',fermat_flag=False, save_flag=False,zp_flag=False):

    #scan_num = sys.argv[1]
    #sid = np.int(scan_num)
    global elem
    elem = element
    global scan_num
    #scan_num = np.str(sid)
    #scan_id, df = _load_scan(sid, fill_events=False)
    scan_num = np.str(sid)
    df = db.get_table(db[sid],fill=False)
    num_frame, count = np.shape(df)
    hdr = db[sid]
    #fermat_flag = np.int(sys.argv[3])
    #elem = sys.argv[2]
    #det_name = sys.argv[4]
    ic = np.asfarray(df['sclr1_ch4'])
    #ic_0 = 153000

    #images = db.get_images(db[sid],name=det_name)
    images = np.array(np.squeeze(list(hdr.data(det_name))))
    print(np.shape(images))
    #mask = 1-io.imread('/data/users/2020Q3/Huang_2020Q3/NPO_mask.tif')
    #mm = np.load('/data/users/2021Q2/Huang_2021Q2/TMA_LCO_60C/mask.npy')
    #mm2 = np.load('/data/users/2021Q2/Huang_2021Q2/TMA_LCO_pristine/mask2.npy')
    #index = np.where(mask == 1)
    #mx = index[0]
    #my = index[1]
    #print(mx)
    #print(my)
    #m_num = np.shape(mx)
    #print('load mask 2')
    #mm = np.load('/data/users/2021Q2/Liu_2021Q2/NZ150_2/mask.npy')
    #mm2 = np.load('/data/users/2021Q2/Liu_2021Q2/NZ150_2/mask2.npy')
    #mm3 = np.load('/data/users/2021Q2/Liu_2021Q2/Z150_1/mask.npy')
    for i in range(num_frame):
        if np.mod(i,500) ==0:
            print('load frame ',i, '/', num_frame)
        #t = np.flipud(images.get_frame(i)[0]).T
        t = np.flipud(images[i,:,:]).T
        #t = t * mask
        t = t / ic[i]
        #t *= (1-mm)
        #t *= (1-mm2)
        ##t = np.flipud(t)
        #t[20,187] = 0
        #t *= (1-mm3)
        #t *= (1-mm3)
        #t *= (1-mm4)
        #t[mm2 == 1.] == 0
        #t[mm3 == 1.] == 0
        #t[mm4 == 1.] == 0
        #for jj in range(m_num[0]):
        #    t = rm_pixel(t,mx[jj],my[jj])
        '''
        plt.figure()
        plt.imshow(t)
        plt.show()
        ddd
        '''
        #if i == 0:
        #    index = np.where(t >= 5)
        #'''
        #t[96,222] = 0.
        #t[140,94] = 0.
        #t[57,138] = 0.
        #'''
        if i == 0:
            nx,ny = np.shape(t)
            global diff_array
            diff_array = np.zeros((nx,ny,num_frame))
        #t[index] = 0
        #t[mask == 1] = 0
        #t[164,107] = 0
        diff_array[:,:,i] = t #* mask


    diff_array[22,255,:] = 0
    #diff_array[296,281,:] = 0
    #diff_array[91,221,:] = 0
    #diff_array[440,381,:] = 0

    #diff_array[diff_array > 200000] = 0
    '''
    #diff_array[diff_array > 1e4] = 0
    diff_array[419,412,:] = 0
    diff_array[431,408,:] = 0
    diff_array[118,370,:] = 0
    diff_array[145,357,:] = 0
    diff_array[142,345,:] = 0
    diff_array[57,462,:] = 0
    diff_array[206,83,:] = 0
    diff_array[115,348,:] = 0
    diff_array[118,359,:] = 0
    diff_array[446,385,:] = 0
    '''
    '''
    diff_array[191,87,:] = 0
    diff_array[124,57,:] = 0
    diff_array[51,152,:] = 0
    diff_array[114,113,:] = 0
    diff_array[111,133,:] = 0
    diff_array[134,122,:] = 0
    '''
    #if elem == 'roi':
    for i in range(num_frame):
        if i == 0:
            global roi
            roi = np.zeros(num_frame)
        roi[i] = np.sum(diff_array[:,:,i])

    global xrf
    if elem in df:
        xrf = np.asarray(df[elem])
    else:
        xrf = np.asfarray(eval('df.Det1_'+elem)) + np.asarray(eval('df.Det2_'+elem)) + np.asarray(eval('df.Det3_'+elem))
    xrf = xrf * ic[0] / (ic + 1.e-9)


    if fermat_flag:
        x = np.asarray(df.dssx)
        y = np.asarray(df.dssy)
        fig = plt.figure(1)
        ax = fig.add_subplot(111)
        props = dict(alpha=0.8, edgecolors='none' )
        im = ax.scatter(x,y,c=xrf,s=50,marker='s',**props)
        plt.axes().set_aspect('equal')
        plt.gca().invert_yaxis()

        cid = fig.canvas.mpl_connect('button_press_event', onclick_fermat)
        fig2 = plt.figure(2)
        ax2 = fig2.add_subplot(111)
    else:
        global y_data
        global x_data
        global num_x
        global num_y
        try:
            hdr.start.plan_args['num']
            #print(hdr.start.plan_args['num'])
            xrf = np.reshape(xrf,(1,hdr.start.plan_args['num']))
            roi = np.reshape(roi,(1,hdr.start.plan_args['num']))
            #roi2 = np.reshape(roi2,(1,hdr.start.plan_args['num']))
        except:
            if hdr.start.plan_name == 'grid_scan':
                xrf = np.reshape(xrf,(hdr.start.shape[0],hdr.start.shape[1]))
                roi = np.reshape(roi,(hdr.start.shape[0],hdr.start.shape[1]))
                x_data = df[hdr.start.motors[1]]
                y_data = df[hdr.start.motors[0]]
                num_x = hdr.start.shape[0]
                num_y = hdr.start.shape[1]
                #roi2 = np.reshape(roi2,(hdr.start.shape[0],hdr.start.shape[1]))
                extent = (hdr.start.plan_args['args'][2], hdr.start.plan_args['args'][1],hdr.start.plan_args['args'][6],hdr.start.plan_args['args'][5])
            elif hdr.start.plan_name == 'FlyPlan2D':
                xrf = np.reshape(xrf,(hdr.start.shape[1],hdr.start.shape[0]))
                roi = np.reshape(roi,(hdr.start.shape[1],hdr.start.shape[0]))
                #roi2 = np.reshape(roi2,(hdr.start.shape[1],hdr.start.shape[0]))
                extent = (hdr.start.plan_args['scan_end1'], hdr.start.plan_args['scan_start1'],hdr.start.plan_args['scan_end2'],hdr.start.plan_args['scan_start2'])
                #x_motor = hdr['motor1']
                if zp_flag:
                    x_motor = hdr.start['motors'][0]
                    x_data = np.asarray(df[x_motor])
                    #y_motor = hdr['motor2']
                    y_motor = hdr.start['motors'][1]
                    y_data = np.asarray(df[y_motor])
                else:
                    x_motor = hdr.start['motors'][0]
                    x_data = np.asarray(df[x_motor])
                    y_motor = hdr.start['motors'][1]
                    y_data = np.asarray(df[y_motor])
                #global num_x
                #global num_y
                num_x = hdr.start.shape[0]
                num_y = hdr.start.shape[1]
            else:
                x_motor = hdr['motor1']
                x_data = np.asarray(df[x_motor])
                y_motor = hdr['motor2']
                y_data = np.asarray(df[y_motor])
                xrf = np.reshape(xrf,(hdr.start.shape[1],hdr.start.shape[0]))
                roi = np.reshape(roi,(hdr.start.shape[1],hdr.start.shape[0]))
                #roi2 = np.reshape(roi2,(hdr.start.shape[1],hdr.start.shape[0]))
                extent = (np.nanmin(x_data), np.nanmax(x_data),np.nanmax(y_data), np.nanmin(y_data))
                #print('no num')

    fig=plt.figure(0)
    ax = fig.add_subplot(111)
    im = ax.imshow(xrf,interpolation = 'none',aspect='auto')
    plt.title('#'+scan_num +' '+ elem)
    global w
    w = xrf.shape[1]
    cid = fig.canvas.mpl_connect('button_press_event', onclick)

    fig1 = plt.figure(1)
    ax1 = fig1.add_subplot(111)
    im = ax1.imshow(roi,interpolation = 'none',aspect='auto')
    plt.title('#'+scan_num+' Diff ROI')
    #w = roi.shape[1]
    cid = fig1.canvas.mpl_connect('button_press_event', onclick_roi)


    fig3 = plt.figure(3)
    ax3 = fig3.add_subplot(111)

    plt.show()

    if save_flag:
        io.imsave('/data/users/2021Q3/Singer_2021Q3/rock_'+scan_num+'_roi.tif',roi.astype(np.float32))
        io.imsave('/data/users/2021Q3/Singer_2021Q3/rock_'+scan_num+'_xrf.tif',xrf.astype(np.float32))
        io.imsave('/data/users/2021Q3/Singer_2021Q3/rock_'+scan_num+'_diff_data.tif',diff_array.astype(np.float32))


