import logging
import scipy
import pickle
import json
import numpy as np
import pandas as pd

from scipy.optimize import curve_fit
from epics import caget, caput

logging.basicConfig(level=logging.INFO)
logger = logging.getLogger()

def erfunc3(z,a,b,c,d,e):
    return d+e*z+c*(scipy.special.erf((z-a)/(b*np.sqrt(2.0)))+1.0)
def erfunc4(z,a,b,c,d,e):
    return d+e*z+c*(1.0-scipy.special.erf((z-a)/(b*np.sqrt(2.0))))

def erfunc1(z,a,b,c):
    return c*(scipy.special.erf((z-a)/(b*np.sqrt(2.0)))+1.0)
def erfunc2(z,a,b,c):
    return c*(1.0-scipy.special.erf((z-a)/(b*np.sqrt(2.0))))
def sicifunc(z,a,b,c):
    si, ci = scipy.special.sici(2.0*b*(z-a))
    #return c*((-1.0+np.cos(2.0*b*(a-z))+2.0*b*(a-z)*si)/(2.0*(b*b)*(z-a))+np.pi/(2.0*b))*b
    return c*(-scipy.sinc(b*(z-a)/np.pi)*scipy.sin(b*(z-a))+si+np.pi/2.0)/np.pi
def squarefunc(z,c,a1,b1,a2,b2):
    return c*(scipy.special.erf((z-a1)/(b1*np.sqrt(2.0)))-scipy.special.erf((z-a2)/(b2*np.sqrt(2.0))))
def erf_fit(sid,elem,mon='sclr1_ch4',linear_flag=True):
    h=db[sid]
    sid=h['start']['scan_id']
    df=h.table()
    mots=h.start['motors']
    xdata=df[mots[0]]
    xdata=np.array(xdata,dtype=float)
    ydata=(df['Det1_'+elem]+df['Det2_'+elem]+df['Det3_'+elem])/df[mon]
    ydata=np.array(ydata,dtype=float)
    ydata[ydata==np.nan] = np.nanmean(ydata)#patch for ic3 returns zero
    ydata[ydata==np.inf] = np.nanmean(ydata)#patch for ic3 returns zero
    ydata[0] = ydata[1] #patch for drop point issue
    ydata[-1] = ydata[-2]#patch for drop point issue
    y_min=np.min(ydata)
    y_max=np.max(ydata)
    ydata=(ydata-y_min)/(y_max-y_min)
    plt.figure()
    plt.plot(xdata,ydata,'bo')
    y_mean = np.mean(ydata)
    half_size = int (len(ydata)/2)
    y_half_mean = np.mean(ydata[0:half_size])
    edge_pos=find_edge(xdata,ydata,10)
    if y_half_mean < y_mean:
        if linear_flag == False:
            popt,pcov=curve_fit(erfunc1,xdata,ydata, p0=[edge_pos,0.05,0.5])
            fit_data=erfunc1(xdata,popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(erfunc3,xdata,ydata, p0=[edge_pos,0.05,0.5,0,0])
            fit_data=erfunc3(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]);
    else:
        if linear_flag == False:
            popt,pcov=curve_fit(erfunc2,xdata,ydata,p0=[edge_pos,0.05,0.5])
            fit_data=erfunc2(xdata,popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(erfunc4,xdata,ydata,p0=[edge_pos,0.05,0.5,0,0])
            fit_data=erfunc4(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]);
    plt.plot(xdata,fit_data)
    plt.title(f'{sid = }, edge = {popt[0] :.3f}, FWHM = {popt[1]*2354.8 :.2f} nm \n {zp.zpz1.position = :.3f}')
    return (popt[0],popt[1]*2.3548*1000.0)


def sici_fit(sid,elem,mon='sclr1_ch4',linear_flag=True):
    h=db[sid]
    sid=h['start']['scan_id']
    df=h.table()
    mots=h.start['motors']
    xdata=df[mots[0]]
    xdata=np.array(xdata,dtype=float)
    ydata=(df['Det1_'+elem]+df['Det2_'+elem]+df['Det3_'+elem])/df[mon]
    ydata=np.array(ydata,dtype=float)
    y_min=np.min(ydata)
    y_max=np.max(ydata)
    ydata=(ydata-y_min)/(y_max-y_min)
    plt.figure()
    plt.plot(xdata,ydata,'bo')
    y_mean = np.mean(ydata)
    half_size = int (len(ydata)/2)
    y_half_mean = np.mean(ydata[0:half_size])
    edge_pos=find_edge(xdata,ydata,10)
    if y_half_mean < y_mean:
        if linear_flag == False:
            popt,pcov=curve_fit(sicifunc,xdata,ydata, p0=[edge_pos,0.05,1/np.pi])
            fit_data=sicifunc(xdata,popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(sicifunc,xdata,ydata, p0=[edge_pos,0.05,1/np.pi,0,0])
            fit_data=sicifunc(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]);
    else:
        if linear_flag == False:
            popt,pcov=curve_fit(sicifunc,-xdata,ydata,p0=[-edge_pos,0.05,1/np.pi])
            fit_data=sicifunc(-xdata,-popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(sicifunc,-xdata,ydata,p0=[-edge_pos,0.05,1/np.pi,0,0])
            fit_data=sicifunc(-xdata,-popt[0],popt[1],popt[2],popt[3],popt[4]);
    plt.plot(xdata,fit_data)
    plt.title('sid= %d edge = %.3f, FWHM = %.2f nm' % (sid,popt[0], 1.39156*2*1000.0/popt[1]))
    return (popt[0],1.39156*2*1000.0/popt[1])


def find_double_edge2(xdata,ydata):
    l = np.size(ydata)
    der = np.zeros(l-1)
    for i in range(l-1):
        der[i]=ydata[i+1]-ydata[i]
    ind1 = scipy.argmax(der)
    ind2 = scipy.argmin(der)
    return(xdata[ind1],xdata[ind2])

def square_fit(sid,elem,mon='sclr1_ch4',linear_flag=True):

    h=db[sid]
    sid=h['start']['scan_id']
    df=h.table()
    mots=h.start['motors']
    xdata=df[mots[0]]
    xdata=np.array(xdata,dtype=float)
    #df[mon][df[mon]==np.inf] = np.mean(df[mon])
    #ydata=(df['Det1_'+elem]+df['Det2_'+elem]+df['Det3_'+elem])/(df[mon]+1e-8)
    ydata=(df['Det1_'+elem]+df['Det2_'+elem]+df['Det3_'+elem])

    ydata=np.array(ydata,dtype=float)
    y_min=np.min(ydata)
    y_max=np.max(ydata)
    ydata=(ydata-y_min)/(y_max-y_min)
    plt.figure()
    plt.plot(xdata,ydata,'bo')
    edge_pos_1, edge_pos_2 = find_double_edge(xdata,ydata,10)
    #print('sid={}  e1={}  e2={}'.format(sid,edge_pos_1,edge_pos_2))
    popt,pcov=curve_fit(squarefunc,xdata,ydata,p0=[0.5,edge_pos_1,0.1,edge_pos_2,0.1])
    fit_data=squarefunc(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]);

    #print('a={} b={} c={}'.format(popt[0],popt[1],popt[2]))
    plt.plot(xdata,fit_data)
    plt.title('sid= %d cen = %.3f e1 = %.3f e2 = %.3f ' % (sid,(popt[1]+popt[3])*0.5, popt[1],popt[3]))
    plt.xlabel(mots[0])
    return (popt[1],popt[3],(popt[1]+popt[3])*0.5)
    #return(xdata, ydata, fit_data)



def data_erf_fit(sid,xdata,ydata,linear_flag=True):

    xdata=np.array(xdata,dtype=float)
    ydata=np.array(ydata,dtype=float)
    y_min=np.min(ydata)
    y_max=np.max(ydata)
    ydata=(ydata-y_min)/(y_max-y_min)
    plt.figure()
    plt.plot(xdata,ydata,'bo')
    y_mean = np.mean(ydata)
    half_size = int (len(ydata)/2)
    y_half_mean = np.mean(ydata[0:half_size])
    edge_pos=find_edge(xdata,ydata,10)
    if y_half_mean < y_mean:
        if linear_flag == False:
            popt,pcov=curve_fit(erfunc1,xdata,ydata, p0=[edge_pos,0.5,0.5])
            fit_data=erfunc1(xdata,popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(erfunc3,xdata,ydata, p0=[edge_pos,0.5,0.5,0,0])
            fit_data=erfunc3(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]);
    else:
        if linear_flag == False:
            popt,pcov=curve_fit(erfunc2,xdata,ydata,p0=[edge_pos,0.5,0.5])
            fit_data=erfunc2(xdata,popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(erfunc4,xdata,ydata,p0=[edge_pos,0.5,0.5,0,0])
            fit_data=erfunc4(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]);

    #print('a={} b={} c={}'.format(popt[0],popt[1],popt[2]))
    #plt.figure(1000)
    plt.plot(xdata,fit_data)
    plt.title('sid = %d, edge = %.3f, FWHM = %.2f nm' % (sid, popt[0], popt[1]*2.3548*1000.0))
    return (popt[0],popt[1]*2.3548*1000.0)
def data_sici_fit(xdata,ydata,linear_flag=True):

    xdata=np.array(xdata,dtype=float)
    ydata=np.array(ydata,dtype=float)
    y_min=np.min(ydata)
    y_max=np.max(ydata)
    ydata=(ydata-y_min)/(y_max-y_min)
    plt.figure(1000)
    plt.plot(xdata,ydata,'bo')
    y_mean = np.mean(ydata)
    half_size = int (len(ydata)/2)
    y_half_mean = np.mean(ydata[0:half_size])
    edge_pos=find_edge(xdata,ydata,10)
    if y_half_mean < y_mean:
        if linear_flag == False:
            popt,pcov=curve_fit(sicifunc,xdata,ydata, p0=[edge_pos,20,1])
            fit_data=sicifunc(xdata,popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(sicifunc,xdata,ydata, p0=[edge_pos,20,1])
            fit_data=sicifunc(xdata,popt[0],popt[1],popt[2]);
    else:
        if linear_flag == False:
            popt,pcov=curve_fit(sicifunc,-xdata,ydata,p0=[-edge_pos,20,1])
            fit_data=sicifunc(-xdata,-popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(sicifunc,-xdata,ydata,p0=[-edge_pos,20,1])
            fit_data=sicifunc(-xdata,-popt[0],popt[1],popt[2]);
    plt.plot(xdata,fit_data)
    plt.title('edge = %.3f, FWHM = %.2f nm' % (popt[0], 1.39156*2*1000.0/popt[1]))
    return (popt[0],1.39156*2*1000.0/popt[1])

def find_2D_edge(sid, axis, elem):
    df2 = db.get_table(db[sid],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    motors = db[sid].start['motors']
    x = np.array(df2[motors[0]])
    y = np.array(df2[motors[1]])
    #I0 = np.asfarray(df2.sclr1_ch4)
    I0 = np.asfarray(df2['sclr1_ch4'])
    scan_info=db[sid]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']
    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))
    #l = np.linspace(y[0],y[-1],ny)
    #s = xrf.sum(1)
    if axis == 'x':
        l = np.linspace(x[0],x[-1],nx)
        s = xrf.sum(0)
    else:
        l = np.linspace(y[0],y[-1],ny)
        s = xrf.sum(1)
    edge,fwhm = data_erf_fit(l, s)
    return edge


def mll_z_alignment(z_start, z_end, z_num, mot, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_sz = smlld.sbz.position
    yield from bps.movr(smlld.sbz, z_start)
    for i in range(z_num + 1):

        yield from fly1d(dets1, mot, start, end, num, acq_time)

        #plot(-1, elem, mon)
        #plt.title('sbz = %.3f' % smlld.sbz.position)
        '''
        h=db[-1]
        sid=h['start']['scan_id']
        df=db.get_table(h)
        xdata=df[mot]
        xdata=np.array(xdata,dtype=float)
        x_mean=np.mean(xdata)
        xdata=xdata-x_mean
        ydata=(df['Det1_'+elem]+df['Det2_'+elem]+df['Det3_'+elem])/df[mon]
        ydata=np.array(ydata,dtype=float)
        y_min=np.min(ydata)
        y_max=np.max(ydata)
        ydata=(ydata-y_min)/y_max
        y_mean = np.mean(ydata)
        half_size = int (len(ydata)/2)
        y_half_mean = np.mean(ydata[0:half_size])
        if y_half_mean < y_mean:
            popt,pcov=curve_fit(erfunc1,xdata,ydata)
            fit_data=erfunc1(xdata,popt[0],popt[1],popt[2]);
        else:
            popt,pcov=curve_fit(erfunc2,xdata,ydata)
            fit_data=erfunc2(xdata,popt[0],popt[1],popt[2]);
        plt.figure()
        plt.plot(xdata,ydata,'bo')
        plt.plot(xdata,fit_data)
        z_pos[i]=smlld.sbz.position
        fit_size[i]=popt[1]*2.3548*1000
        plt.title('sid = %d sbz = %.3f um FWHM = %.2f nm' %(sid,smlld.sbz.position,fit_size[i]))
        '''
        edge_pos,fwhm=erf_fit(-1,elem,mon,linear_flag=False)
        fit_size[i]= fwhm
        z_pos[i]=smlld.sbz.position
        yield from bps.movr(smlld.sbz, z_step)
    yield from bps.mov(smlld.sbz, init_sz)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')
    plt.xlabel('sbz')

def find_edge(xdata,ydata,size):
    set_point=0.5
    j=int (ceil(size/2.0))
    l=len(ydata)
    local_mean=np.zeros(l-size)
    for i in range(l-size):
        local_mean[i]=np.mean(ydata[i:i+size])
    zdata=abs(local_mean-np.array(set_point))
    index=np.argmin(zdata)
    index=index+j
    return xdata[index]

def find_double_edge(xdata, ydata, size):
    edge_1 = find_edge(xdata, ydata, size)
    l = np.size(ydata)
    index = np.argmax(ydata)
    cen = xdata[index]
    if cen > edge_1:
        edge_2 = find_edge(xdata[index:l],ydata[index:l],size)
        #edge_2 = (cen-edge_1) + cen
        return(edge_1,edge_2)
    else:
        #edge_2 = cen - (edge_1 - cen)
        edge_2 = find_edge(xdata[1:index],ydata[1:index],size)
        return(edge_2,edge_1)
def hmll_z_alignment(z_start, z_end, z_num, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_hz = hmll.hz.position
    yield from bps.movr(hmll.hz, z_start)
    for i in range(z_num + 1):
        yield from fly1d(dets1,dssx, start, end, num, acq_time)
        edge_pos,fwhm=erf_fit(-1,elem,mon)
        fit_size[i]=fwhm
        z_pos[i]=hmll.hz.position
        yield from bps.mov(dssx, edge_pos)
        yield from bps.movr(hmll.hz, z_step)
    yield from bps.mov(hmll.hz, init_hz)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')
    plt.xlabel('hz')

def mll_vchi_alignment(vchi_start, vchi_end, vchi_num, mot, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    vchi_pos = np.zeros(vchi_num+1)
    fit_size = np.zeros(vchi_num+1)
    vchi_step = (vchi_end - vchi_start)/vchi_num
    init_vchi = vmll.vchi.position
    yield from bps.movr(vmll.vchi, vchi_start)
    for i in range(vchi_num + 1):
        yield from fly1d(dets1, mot, start, end, num, acq_time,dead_time=0.002)
        edge_pos,fwhm=erf_fit(-1,elem,mon)
        fit_size[i]=fwhm
        vchi_pos[i]=vmll.vchi.position
        yield from bps.movr(vmll.vchi, vchi_step)
    yield from bps.mov(vmll.vchi, init_vchi)
    plt.figure()
    plt.plot(vchi_pos,fit_size,'bo')
    plt.xlabel('vchi')

def vmll_z_alignment(z_start, z_end, z_num, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_vz = vmll.vz.position
    yield from bps.movr(vmll.vz, z_start)
    for i in range(z_num + 1):
        yield from fly1d(dets1,dssy, start, end, num, acq_time)
        edge_pos,fwhm=erf_fit(-1,elem,mon)
        #plt.title('vz={}'.format(vmll.vz.position),loc='right')
        fit_size[i]=fwhm
        z_pos[i]=vmll.vz.position
        yield from bps.movr(vmll.vz, z_step)
    yield from bps.mov(vmll.vz, init_vz)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')
    plt.xlabel('vz')

def zp_z_alignment(z_start, z_end, z_num, mot, start, end, num, acq_time, elem=' ',linFlag = True,mon='sclr1_ch4'):
    
    print("moves the zone plate relatively and find the focus with a linescan at each position")
    
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_sz = zp.zpz1.position

    yield from movr_zpz1(z_start)
    for i in range(z_num + 1):

        yield from fly1d(dets_fs, mot, start, end, num, acq_time)
        edge_pos,fwhm=erf_fit(-1,elem,mon,linear_flag=linFlag)
        fit_size[i]= fwhm
        z_pos[i]=zp.zpz1.position
        yield from movr_zpz1(z_step)
        merlin1.unstage()
        xspress3.unstage()
    yield from movr_zpz1(-1*z_end)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')
    plt.xlabel('zpz1')


def pos2angle(col,row):

    # old Dexelar calibration
    '''
    pix = 74.8
    R = 2.315e5
    th1 = 0.7617
    phi1 = 3.0366
    th2 = 0.1796
    phi2 = 2.5335
    phi3 = -0.1246
    alpha = 8.5*np.pi/180
    '''
    # new Dexelar calibration at position 0
    pix = 74.8
    R = 2.6244e5
    th1 = 0.7685
    phi1 = 3.0238
    th2 = 0.1398
    phi2 = 2.9292
    phi3 = -0.1486
    alpha = 8.5*np.pi/180

    det_orig = R*np.array([np.sin(th1)*np.cos(phi1),np.sin(th1)*np.sin(phi1),np.cos(th1)])
    det_z = np.array([np.sin(th2)*np.cos(phi2), np.sin(th2)*np.sin(phi2),np.cos(th2)])
    th3 = np.arctan(-1.0/(np.cos(phi2-phi3)*np.tan(th2)))
    det_x = np.array([np.sin(th3)*np.cos(phi3),np.sin(th3)*np.sin(phi3),np.cos(th3)])
    det_y = np.cross(det_z,det_x)

    pos = det_orig + (col - 1)*pix*det_x + (row -1)*pix*det_y

    M = np.array([[np.cos(alpha),-np.sin(alpha),0],[np.sin(alpha), np.cos(alpha),0],[0,0,1]])

    pos = np.dot(M,pos)

    tth = np.arccos(pos[2]/np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2))*180.0/np.pi
    delta = np.arcsin(pos[1]/np.sqrt(pos[0]**2+pos[1]**2+pos[2]**2))*180.0/np.pi
    pos_xy = pos*np.array([1,0,1])
    gamma = np.arccos(pos[2]/np.sqrt(pos_xy[0]**2+pos_xy[1]**2+pos_xy[2]**2))*180.0/np.pi
    return (gamma,delta,tth)


def return_line_center(sid,elem='Cr',threshold=0.2, neg_flag=0):
    h = db[sid]

    df2 = h.table()
    xrf = np.array(df2['Det2_' + elem]+df2['Det1_' + elem] + df2['Det3_' + elem])

    xrf[xrf==np.nan] = np.nanmean(xrf)#patch for ic3 returns zero
    xrf[xrf==np.inf] = np.nanmean(xrf)#patch for ic3 returns zero

    #threshold = np.max(xrf)/10.0
    x_motor = h.start['motors']
    x = np.array(df2[x_motor[0]])
    if neg_flag == 1:
        xrf = xrf * -1
        xrf = xrf - np.min(xrf)

    #print(x)
    #print(xrf)
    xrf[xrf<(np.max(xrf)*threshold)] = 0.
    #index = np.where(xrf == 0.)
    #xrf[:index[0][0]] = 0.
    xrf[xrf>=(np.max(xrf)*threshold)] = 1.
    mc = find_mass_center_1d(xrf[:-2],x[:-2])
    return mc


def return_tip_pos(sid,elem='Cr'):
    h = db[sid]

    df2 = h.table()
    xrf = np.array(df2['Det2_' + elem]+df2['Det1_' + elem] + df2['Det3_' + elem])
    threshold = np.max(xrf)/10.0
    x_motor = h.start['motor']
    x = np.array(df2[x_motor])
    #print(x)
    #print(xrf)
    #xrf[xrf<(np.max(xrf)*0.5)] = 0.
    #xrf[xrf>=(np.max(xrf)*0.5)] = 1.
    #mc = find_mass_center_1d(xrf,x)
    xrf_d = np.diff(xrf)
    #peak_index = np.where(xrf_d == np.max(xrf_d))
    peak_index = np.where(xrf == np.max(xrf))
    #print(x[peak_index[0][0]+1])
    return x[peak_index[0][0]+1]

def zp_rot_alignment_edge(a_start, a_end, a_num, start, end, num, acq_time, elem='Pt_L', move_flag=0):
    a_step = (a_end - a_start)/a_num
    x = np.zeros(a_num+1)
    y = np.zeros(a_num+1)
    orig_th = zps.zpsth.position
    for i in range(a_num+1):
        x[i] = a_start + i*a_step
        yield from bps.mov(zps.zpsth, x[i])
        if np.abs(x[i]) > 45:
            yield from fly1d(dets1,zpssz,start,end,num,acq_time)
            #tmp = return_line_center(-1, elem=elem,threshold=0.5)
            edge,fwhm = erf_fit(-1,elem = elem,linear_flag=False)
            y[i] = edge*np.sin(x[i]*np.pi/180.0)
        else:
            yield from fly1d(dets1,zpssx,start,end,num,acq_time)
            #tmp = return_line_center(-1,elem=elem,threshold=0.5)
            edge,fwhm = erf_fit(-1,elem = elem,linear_flag=False)
            y[i] = edge*np.cos(x[i]*np.pi/180.0)
        print('y=',y[i])
    y = -1*np.array(y)
    x = np.array(x)
    r0, dr, offset = rot_fit_2(x,y)
    yield from bps.mov(zps.zpsth, orig_th)
    dx = -dr*np.sin(offset*np.pi/180)/1000.0
    dz = -dr*np.cos(offset*np.pi/180)/1000.0

    print('dx=',dx,'   ', 'dz=',dz)

    if move_flag:
        yield from bps.movr(zps.smarx, dx)
        yield from bps.movr(zps.smarz, dz)


    return x,y

def zp_rot_alignment(a_start, a_end, a_num, start, end, num, acq_time, elem='Pt_L', neg_flag = 0, move_flag=0):
    a_step = (a_end - a_start)/a_num
    x = np.zeros(a_num+1)
    y = np.zeros(a_num+1)
    orig_th = zps.zpsth.position
    for i in range(a_num+1):
        x[i] = a_start + i*a_step
        yield from bps.mov(zps.zpsth, x[i])
        if np.abs(x[i]) > 45.05:
            yield from fly1d(dets_fs,zpssz,start,end,num,acq_time)
            tmp = return_line_center(-1, elem=elem,threshold=0.8,neg_flag=neg_flag)
            #tmp = return_tip_pos(-1, elem=elem)
            #tmp,fwhm = erf_fit(-1,elem = elem,linear_flag=False)
            y[i] = tmp*np.sin(x[i]*np.pi/180.0)
        else:
            yield from fly1d(dets_fs,zpssx,start,end,num,acq_time)
            tmp = return_line_center(-1,elem=elem,threshold=0.5,neg_flag=neg_flag )
            #tmp = return_tip_pos(-1, elem=elem)
            #tmp,fwhm = erf_fit(-1,elem = elem,linear_flag=False)
            y[i] = tmp*np.cos(x[i]*np.pi/180.0)
        print('y=',y[i])
    y = -1*np.array(y)
    x = np.array(x)
    r0, dr, offset = rot_fit_2(x,y)
    yield from bps.mov(zps.zpsth, orig_th)
    dx = -dr*np.sin(offset*np.pi/180)/1000.0
    dz = -dr*np.cos(offset*np.pi/180)/1000.0

    print('dx=',dx,'   ', 'dz=',dz)

    if move_flag:
        yield from bps.movr(zps.smarx, dx)
        yield from bps.movr(zps.smarz, dz)


    return x,y




def mll_rot_alignment(a_start, a_end, a_num, start, end, num, acq_time, elem='Pt_L', move_flag=0):
    th_init = smlld.dsth.position
    y_init = dssy.position
    #y_init = -0.5 #remove this temp.
    a_step = (a_end - a_start)/a_num
    x = np.zeros(a_num+1)
    y = np.zeros(a_num+1)
    v = np.zeros(a_num+1)
    orig_th = smlld.dsth.position
    for i in range(a_num+1):
        yield from bps.mov(dssx,0)
        yield from bps.mov(dssz,0)
        x[i] = a_start + i*a_step
        yield from bps.mov(smlld.dsth, x[i])
        #angle = smlld.dsth.position
        #dy = -0.1+0.476*np.sin(np.pi*(angle*np.pi/180.0-1.26)/1.47)
        #ddy = (-0.0024*angle)-0.185
        #dy = dy+ddy
        #yield from bps.movr(dssy,dy)

        y_offset1 = sin_func(x[i], 0.110, -0.586, 7.85,1.96)
        y_offset2 = sin_func(th_init, 0.110, -0.586, 7.85,1.96)
        yield from bps.mov(dssy,y_init+y_offset1-y_offset2)


        if np.abs(x[i]) > 45.01:
            #yield from fly2d(dets1,dssz,start,end,num, dssy, -2,2,20,acq_time)
            #cx,cy = return_center_of_mass(-1,elem,0.3)
            #y[i] = cx*np.sin(x[i]*np.pi/180.0)
            yield from fly1d(dets1,dssz,start,end,num,acq_time)
            cen = return_line_center(-1, elem=elem,threshold = 0.3)
            #cen, edg1, edg2 = square_fit(-1,elem=elem)
            y[i] = cen*np.sin(x[i]*np.pi/180.0)
            # yield from bps.mov(dssz,cen)
        else:
            #yield from fly2d(dets1,dssx,start,end,num, dssy, -2,2,20,acq_time)
            #cx,cy = return_center_of_mass(-1,elem,0.3)
            #y[i] = cx*np.cos(x[i]*np.pi/180.0)
            yield from fly1d(dets1,dssx,start,end,num,acq_time)
            cen = return_line_center(-1,elem=elem,threshold = 0.3)
            #cen, edg1, edg2 = square_fit(-1,elem=elem)
            y[i] = cen*np.cos(x[i]*np.pi/180.0)
            #y[i] = tmp*np.cos(x[i]*np.pi/180.0)
            #y[i] = -tmp*np.cos(x[i]*np.pi/180.0)
            # yield from bps.mov(dssx,cen)

        ##v[i] = cy

        #yield from bps.mov(dssy,cy)
        #yield from fly1d(dets1,dssy,-2,2,100,acq_time)
        #tmp = return_line_center(-1, elem=elem)
        #yield from bps.mov(dssy,tmp)
        #v[i] = tmp
        #print('h_cen= ',y[i],'v_cen = ',v[i])
        #plot_data(-1,elem,'sclr1_ch4')
        #insertFig(note='dsth = {}'.format(check_baseline(-1,'dsth')))
        #plt.close()

    y = -1*np.array(y)
    x = np.array(x)
    r0, dr, offset = rot_fit_2(x,y)
    #insertFig(note='dsth: {} {}'.format(a_start,a_end))
    
    yield from bps.mov(smlld.dsth, th_init)
    
    dx = -dr*np.sin(offset*np.pi/180)
    dz = -dr*np.cos(offset*np.pi/180)
    
    #moving back to intial y position
    yield from bps.mov(dssy, y_init)
    print(f'{dx = :.2f}, {dz = :.2f}')

    #if move_flag:
        #yield from bps.movr(smlld.dsx, dx)
        #yield from bps.movr(smlld.dsz, dz)

    #plt.figure()
    #plt.plot(x,v)

    x = np.array(x)
    y = -np.array(y)

    print(x)
    print(y)
    print(v)
    #caliFIle = open('rotCali','wb')
    #pickle.dump(y,CaliFile)
    return v

def mll_rot_alignment_2D(th_start, th_end, th_num, x_start, x_end, x_num,
                         y_start, y_end, y_num, acq_time, elem='Pt_L', move_flag=0):
    
    th_list = np.linspace(th_start, th_end, th_num+1)
    
    x = th_list
    y = np.zeros(th_num+1)
    v = np.zeros(th_num+1)
    orig_th = smlld.dsth.position
    for i, th in enumerate(th_list):
        yield from bps.mov(dssx,0, dssz,0, smlld.dsth,th)

        if np.abs(x[i]) > 45.01:
            yield from fly2d(dets1,
                            dssz,
                            x_start,
                            x_end,
                            x_num, 
                            dssy, 
                            y_start,
                            y_end,
                            y_num,
                            acq_time
                            )

            cx,cy = return_center_of_mass(-1,elem,0.5)
            y[i] = cx*np.sin(x[i]*np.pi/180.0)

        else:
            yield from fly2d(dets1,dssx,start,end,num, dssy, -0.5,0.5,20,acq_time)
            cx,cy = return_center_of_mass(-1,elem,0.5)
            y[i] = cx*np.cos(x[i]*np.pi/180.0)

        yield from bps.mov(dssy,cy)

    y = -1*np.array(y)
    x = np.array(x)
    r0, dr, offset = rot_fit_2(x,y)
    yield from bps.mov(smlld.dsth, 0)
    dx = -dr*np.sin(offset*np.pi/180)
    dz = -dr*np.cos(offset*np.pi/180)

    print('dx=',dx,'   ', 'dz=',dz)

    x = np.array(x)
    y = -np.array(y)
    print(x)
    print(y)
    #caliFIle = open('rotCali','wb')
    #pickle.dump(y,CaliFile)


def mll_rot_v_alignment(a_start, a_end, a_num, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    a_step = (a_end - a_start)/a_num
    x = np.zeros(a_num+1)
    y = np.zeros(a_num+1)
    orig_th = smlld.dsth.position
    for i in range(a_num+1):
        x[i] = a_start + i*a_step
        yield from bps.mov(smlld.dsth, x[i])
        yield from fly1d(dets1,dssy,start,end,num,acq_time)
        edge_pos,fwhm=erf_fit(-1,elem=elem,mon=mon)
        y[i] = edge_pos
    y = -1*np.array(y)
    x = np.array(x)
    yield from bps.mov(smlld.dsth,0)
    #r0, dr, offset = rot_fit_2(x,y)
    #yield from bps.mov(smlld.dsth, 0)
    #dx = -dr*np.sin(offset*np.pi/180)
    #dz = -dr*np.cos(offset*np.pi/180)
    print(x,y)
    plt.figure()
    plt.plot(x,y)
    #print('dx=',dx,'   ', 'dz=',dz)

    return x,y

def check_baseline(sid,name):
    h = db[sid]
    bl = h.table('baseline')
    dsmll_list = ['dsx','dsy','dsz','dsth','sbx','sbz','dssx','dssy','dssz']
    vmll_list = ['vx','vy','vz','vchi','vth']
    hmll_list = ['hx','hy','hz','hth']
    mllosa_list = ['osax','osay','osaz']
    mllbs_list = ['mll_bsx','mll_bsy','mll_bsz','mll_bsth']

    if name =='dsmll':
        #print(bl[dsmll_list])
        return(bl[dsmll_list])
    elif name == 'vmll':
        #print(bl[vmll_list])
        return(bl[vmll_list])
    elif name == 'hmll':
        #print(bl[hmll_list])
        return(bl[hmll_list])
    elif name == 'mllosa':
        #print(bl[mllosa_list])
        return(bl[mllosa_list])
    elif name == 'mll':
        #print(bl[dsmll_list])
        #print(bl[vmll_list])
        #print(bl[hmll_list])
        #print(bl[mllosa_list])
        mot_pos = [bl[dsmll_list],bl[vmll_list],bl[hmll_list],bl[mllosa_list]]
        return(mot_pos)
    else:
        #print(name,bl[name])
        return(bl[name].values[0])


def check_info(sid):
    h = db[sid]
    sid = h.start['scan_id']
    scan_time = datetime.fromtimestamp(h.start['time'])
    end_time = datetime.fromtimestamp(h.stop['time'])
    scan_uid = h.start['uid']
    scan_type = h.start['plan_name']
    scan_motors = h.start['motors']
    num_motors = len(scan_motors)
    det_list = h.start['detectors']
    exp_time = h.start['exposure_time']
    print('sid = {}'.format(sid), 'uid = ', scan_uid)
    print('start time = ', scan_time, 'end time = ', end_time)
    if num_motors == 1:
        mot1 = scan_motors[0]
        s1 = h.start['scan_start1']
        e1 = h.start['scan_end1']
        n1 = h.start['num1']
        print(scan_type,mot1,s1,e1,n1,exp_time)
    elif num_motors == 2:
        mot1 = scan_motors[0]
        s1 = h.start['scan_start1']
        e1 = h.start['scan_end1']
        n1 = h.start['num1']
        mot2 = scan_motors[1]
        s2 = h.start['scan_start2']
        e2 = h.start['scan_end2']
        n2 = h.start['num2']
        print(scan_type, mot1,s1,e1,n1,mot2,s2,e2,n2,exp_time)

    print('detectors = ', det_list)

def scan_command(sid):
    h = db[sid]
    sid = h.start['scan_id']
    scan_type = h.start['plan_name']
    if scan_type == 'FlyPlan1D' or scan_type == 'FlyPlan2D':
        scan_motors = h.start['motors']
        num_motors = len(scan_motors)
        exp_time = h.start['exposure_time']
        if num_motors == 1:
            mot1 = scan_motors[0]
            s1 = h.start['scan_start']
            e1 = h.start['scan_end']
            n1 = h.start['num']
            return(mot1+' {:1.3f} {:1.3f} {:d} {:1.3f}'.format(s1,e1,n1,exp_time))
        elif num_motors == 2:
            mot1 = scan_motors[0]
            s1 = h.start['scan_start1']
            e1 = h.start['scan_end1']
            n1 = h.start['num1']
            mot2 = scan_motors[1]
            s2 = h.start['scan_start2']
            e2 = h.start['scan_end2']
            n2 = h.start['num2']
            return(mot1+' {:1.3f} {:1.3f} {:d}'.format(s1,e1,n1)+' '+mot2+' {:1.3f} {:1.3f} {:d} {:1.3f}'.format(s2,e2,n2,exp_time))

class ScanInfo:
    plan = ''
    time = ''
    command = ''
    status = ''
    det = ''
    sid = ''

def scan_info(sid):
    si = ScanInfo()
    h = db[sid]
    si.sid = '{:d}'.format(h.start['scan_id'])
    si.time = datetime.fromtimestamp(h.start['time']).isoformat()
    si.plan = h.start['plan_name']
    si.status = h.stop['exit_status']
    si.command = scan_command(sid)
    si.det = h.start['detectors']
    return(si)



def engage_mirror_feedback():
    
    """
    synchronizes necessary mirror motor positions and reengage the feedbacks
    TODO conditions for enganging and error handling
    
    """

    caput("XF:03IDC-CT{FbPid:01}PID:on",0)
    caput("XF:03IDC-CT{FbPid:02}PID:on",0)

    yield from bps.mov(m1.p,m1.p.position)
    print("HCM Pitch ; Done!")
    yield from bps.mov(m2.p,m2.p.position)
    print("HFM Pitch ; Done!")
    yield from bps.mov(dcm.p,dcm.p.position)
    print("DCM Pitch ; Done!")
    yield from bps.mov(dcm.r,dcm.r.position)
    print("DCM Roll ; Done!")

    print("Engaging feedbacks....")

    m1_p = m1.p.position
    m2_p = m2.p.position

    yield from bps.mov(m2.pf, 10)
    #yield from bps.mov(m1.pf, 10)
    caput("XF:03IDA-OP{HCM:1-Ax:PF}Mtr.VAL",10) #m1_pf
    yield from bps.mov(m1.p,m1_p)
    yield from bps.sleep(3)
    yield from bps.mov(m2.p,m2_p)
    yield from bps.sleep(3)

    caput("XF:03IDC-CT{FbPid:01}PID.I",0) #PID I value to zero
    caput("XF:03IDC-CT{FbPid:02}PID.I",0) #PID I value to zero

    caput("XF:03IDC-CT{FbPid:01}PID:on",1)
    caput("XF:03IDC-CT{FbPid:02}PID:on",1)

    print("Feedbacks Engaged....")


def mll_sample_out_docs():
    """
    vmll vz -8000 and confirm movement
    hmll hz -5000 and confirm movement

    sbz +5000
    sbx -4000
    """

def mll_optics_out_docs():
    """
    move the mll bsx 500 +x 
    move mll bsy 500 -y'

    vmll vy + 500
    hmll hx -500
    osax +2600 um


    """
    
    #ensure fluorescence detecor is out; otherwise move it out 

    #if fdet1.x.position > -30:
    #    yield from bps.mov(fdet1.x, -107)

    #yield from 


def set_motor_val_zero(pv):
    set_pv = pv+".SET"
    val_pv = pv+".VAL"
    caput(set_pv,1)
    caput(val_pv,0)
    caput(set_pv,0)


def zp_optics_to_zero():
    
    zpsx = "XF:03IDC-ES{ZpPI:1-zpsx}Mtr"
    zpsz = "XF:03IDC-ES{ZpPI:1-zpsz}Mtr"
    osax = "XF:03IDC-ES{ANC350:5-Ax:0}Mtr"
    osay = "XF:03IDC-ES{ANC350:5-Ax:1}Mtr"
    bsx = "XF:03IDC-ES{ANC350:8-Ax:0}Mtr"
    bsy = "XF:03IDC-ES{ANC350:8-Ax:1}Mtr"
    bsz = "XF:03IDC-ES{ANC350:8-Ax:2}Mtr"

    list_of_mtrs = [zpsx,zpsz,osax,osay,bsx,bsy,bsz]

    for mtr in list_of_mtrs:
        set_motor_val_zero(mtr)


def save_cam06_images(filename = "crl"):

    pv_filename = epics.PV("XF:03IDC-ES{CAM:06}TIFF1:FileName")
    
    e_ = np.round(e.position,2)
    th = caget("XF:03IDA-OP{Lens:CRL-Ax:P}Mtr.RBV")
    exp_time = caget("XF:03IDC-ES{CAM:06}cam1:AcquireTime_RBV")
    ic1 =sclr2_ch2.get()
    filename_ = f'{filename}_e_{e_}_th_{th :.2f}_exp_{exp_time}_ic1_{ic1}'
    pv_filename.put(filename_)

    for i in range(3):
        print(i)
        time.sleep(2)
        caput('XF:03IDC-ES{CAM:06}TIFF1:WriteFile',1)



def mov_diff(gamma, delta, r=500, calc=0, check_for_dexela = True):
    
    if check_for_dexela and caget("XF:03IDC-ES{Stg:FPDet-Ax:Y}Mtr.RBV")<380:

        raise ValueError("Dexela detector maybe IN, Please move it away and try again!")
        return

    else:

        diff_z = diff.z.position

        gamma = gamma * np.pi / 180
        delta = delta * np.pi / 180
        beta = 89.337 * np.pi / 180

        z_yaw = 574.668 + 581.20 + diff_z
        z1 = 574.668 + 395.2 + diff_z
        z2 = z1 + 380
        d = 395.2

        x_yaw = np.sin(gamma) * z_yaw / np.sin(beta + gamma)
        R_yaw = np.sin(beta) * z_yaw / np.sin(beta + gamma)
        R1 = R_yaw - (z_yaw - z1)
        R2 = R_yaw - (z_yaw - z2)
        y1 = np.tan(delta) * R1
        y2 = np.tan(delta) * R2
        R_det = R1 / np.cos(delta) - d
        dz = r - R_det

        print('Make sure all motors are zeroed properly, '
            'otherwise calculation will be wrong.')
        if x_yaw > 825 or x_yaw < -200:
            print('diff_x = ', -x_yaw,
                ' out of range, move diff_z upstream and try again')
        elif dz < -250 or dz > 0:
            print('diff_cz = ', dz,
                ' out of range, move diff_z up or down stream and try again')
        elif y1 > 750:
            print('diff_y1 = ', y1, ' out of range, move diff_z upstream '
                'and try again')
        elif y2 > 1000:
            print('diff_y2 = ', y2, ' out of range, move diff_z upstream '
                'and try again')
        else:
            print('diff_x = ', -x_yaw, ' diff_cz = ', dz,
                ' diff_y1 = ', y1, ' diff_y2 = ', y2)
            if calc == 0:

                print('wait for 3 sec, hit Ctrl+c to quit the operation')
                yield from bps.sleep(3)
                yield from bps.mov(diff.y1,y1,
                                diff.y2,y2,
                                diff.x,-x_yaw,
                                diff.yaw,gamma*180.0/np.pi,
                                diff.cz,dz)
                '''
                diff.y1.move(y1, wait=False)
                sleep(0.5)
                diff.y2.move(y2, wait=False)
                sleep(0.5)
                diff.x.move(-x_yaw, wait=False)
                sleep(0.5)
                diff.yaw.move(gamma * 180. / np.pi, wait=False)
                sleep(0.5)
                diff.cz.move(dz, wait=False)
                '''
                while (diff.x.moving is True or diff.y1.moving is True or diff.y2.moving is True or diff.yaw.moving is True):
                    yield from bps.sleep(2)
            else:
                print('Calculation mode; no motor will be moved')


def wh_diff():
    diff_z = diff.z.position
    diff_yaw = diff.yaw.position * np.pi / 180.0
    diff_cz = diff.cz.position
    diff_x = diff.x.position
    diff_y1 = diff.y1.position
    diff_y2 = diff.y2.position

    gamma = diff_yaw
    beta = 89.337 * np.pi / 180
    z_yaw = 574.668 + 581.20 + diff_z
    z1 = 574.668 + 395.2 + diff_z
    z2 = z1 + 380
    d = 395.2

    x_yaw = np.sin(gamma) * z_yaw / np.sin(beta + gamma)
    R_yaw = np.sin(beta) * z_yaw / np.sin(beta + gamma)
    R1 = R_yaw - (z_yaw - z1)
    R2 = R_yaw - (z_yaw - z2)

    # print('x_yaw = ', x_yaw, ' diff_x = ', diff_x)
    if abs(x_yaw + diff_x) > 3:
        print('Not a pure gamma rotation')
        return -1,-1,-1

    elif abs(diff_y1 / R1 - diff_y2 / R2) > 0.01:
        print('Not a pure delta rotation')
        return -1,-1,-1
    else:
        delta = np.arctan(diff_y1 / R1)
        R_det = R1 / np.cos(delta) - d + diff_cz

        Gamma = gamma * 180 / np.pi
        Delta = delta * 180 / np.pi
        print(f'{Gamma = :.2f}, {Delta  = :.2f} , r = {R_det :.2f}')
        return Gamma, Delta, R_det


def diff_status():

    gma, delt, r = wh_diff()

    if gma>0 or delt>0 or diff.yaw.position>0.5 or diff.y1.position>23:
        return "diff_pos"

    elif (int(gma) == 0 and int(delt) == 0) or (diff.x.position>10 and not diff.y1.position>23):
        return "safe"


def diff_to_home(move_out_later = False):

    #close c shutter
    caput("XF:03IDC-ES{Zeb:2}:SOFT_IN:B0", 0)

    gma, delt, r = wh_diff()

    if diff_status() == "diff_pos":

        if gma>10:

            if r<450:
                try:
                    yield from mov_diff(gma-10, delt, 500)
                except:
                    pass

            yield from bps.mov(diff.z, -1,diff.y1, 0, diff.y2, 0,diff.cz,-1)
            #yield from bps.mov(diff.y1, 0, diff.y2, 0)
            #yield from bps.mov(diff.cz,-1)

            try:
                yield from mov_diff(0, 0, 500)

            except:
                yield from bps.mov(diff.yaw, 0,diff.x, 0)
                #yield from bps.mov(diff.x, 0)


        else:
            yield from mov_diff(0,0,500)

    else:
        yield from go_det("merlin")
    
    #TODO need to be better to avoid moving back and forth
    if move_out_later:
        yield from go_det('out')





def go_det(det):

    if caget("XF:03IDC-ES{Stg:FPDet-Ax:Y}Mtr.RBV")<380:

        raise ValueError("Dexela detector maybe IN, Please move it away and try again!")
        return
    
    else:

        with open("/nsls2/data/hxn/shared/config/bluesky/profile_collection/startup/diff_det_pos.json") as fp:

            diff_pos = json.load(fp)
            #print(diff_pos)

        merlin_pos = diff_pos["merlin_pos"]
        cam11_pos = diff_pos["cam11_pos"]
        telescope_pos = diff_pos["telescope_pos"]
        out_pos = diff_pos["out"]

        if det == 'merlin':
            #while zposa.zposax.position<20:
            #yield from bps.mov(diff.x, -1.5, diff.y1,-12.9, diff.y2,-12.9, diff.z, -50, diff.cz, -24.7)
            yield from bps.mov(diff.x, merlin_pos['diff_x'], 
                            diff.y1,merlin_pos['diff_y1'], 
                            diff.y2,merlin_pos['diff_y2'], 
                            diff.z, merlin_pos['diff_z'], 
                            diff.cz,merlin_pos['diff_cz']
                            )
            #yield from bps.mov(diff.y1,-3.2)
            #yield from bps.mov(diff.y2,-3.2)
        elif det == 'cam11':
            #yield from bps.mov(diff.x,206.83, diff.y1, 19.177, diff.y2, 19.177,diff.z, -50, diff.cz, -24.7)
            yield from bps.mov(diff.x, cam11_pos['diff_x'], 
                            diff.y1,cam11_pos['diff_y1'], 
                            diff.y2,cam11_pos['diff_y2'], 
                            diff.z, cam11_pos['diff_z'], 
                            diff.cz,cam11_pos['diff_cz']
                            )
            #yield from bps.mov(diff.y1,22.65)
            #yield from bps.mov(diff.y2,22.65)
        elif det =='telescope':
            #yield from bps.mov(diff.x,-342, diff.z, -50, diff.cz, -24.7)
            yield from bps.mov(diff.x,telescope_pos['diff_x'], 
                            diff.z,telescope_pos['diff_z'], 
                            diff.cz,telescope_pos['diff_cz'])
            #yield from bps.mov(diff.z,-50)


        elif det=='out':
            yield from bps.mov(diff.x, out_pos['diff_x'], 
                            diff.y1,out_pos['diff_y1'], 
                            diff.y2,out_pos['diff_y2'], 
                            diff.z, out_pos['diff_z'], 
                            diff.cz,out_pos['diff_cz']
                            )


        else:
            print('Inout det is not defined. '
                'Available ones are merlin, cam11, telescope and tpx')


def update_det_pos(det = "merlin"):
    
    json_path = "/nsls2/data/hxn/shared/config/bluesky/profile_collection/startup/diff_det_pos.json"
    
    with open(json_path, "r") as read_file:
        diff_pos = json.load(read_file)

    if det == "merlin":

        diff_pos['merlin_pos']['diff_x'] = np.round(diff.x.position,2)
        diff_pos['merlin_pos']['diff_y1'] = np.round(diff.y1.position,2)
        diff_pos['merlin_pos']['diff_y2'] = np.round(diff.y2.position,2)
        diff_pos['merlin_pos']['diff_z'] = np.round(diff.z.position,2)
        diff_pos['merlin_pos']['diff_cz'] = np.round(diff.cz.position,2)

    elif det == "cam11":

        diff_pos['cam11_pos']['diff_x'] = np.round(diff.x.position,2)
        diff_pos['cam11_pos']['diff_y1'] = np.round(diff.y1.position,2)
        diff_pos['cam11_pos']['diff_y2'] = np.round(diff.y2.position,2)
        diff_pos['cam11_pos']['diff_z'] = np.round(diff.z.position,2)
        diff_pos['cam11_pos']['diff_cz'] = np.round(diff.cz.position,2)

    elif det == "telescope":

        diff_pos['telescope_pos']['diff_x'] = np.round(diff.x.position,2)
        diff_pos['merlin_pos']['diff_z'] = np.round(diff.z.position,2)
        diff_pos['merlin_pos']['diff_cz'] = np.round(diff.cz.position,2)

    elif det == "out":
        diff_pos['out']['diff_x'] = np.round(diff.x.position,2)
        diff_pos['out']['diff_y1'] = np.round(diff.y1.position,2)
        diff_pos['out']['diff_y2'] = np.round(diff.y2.position,2)
        diff_pos['out']['diff_z'] = np.round(diff.z.position,2)
        diff_pos['out']['diff_cz'] = np.round(diff.cz.position,2)

    else:
        raise KeyError ("Undefined detector name")
    
    read_file.close()
    ext = datetime.now().strftime('%Y-%m-%d')
    json_path_backup = f"/data/users/backup_params/diff_pos/{ext}_diff_det_pos.json"
    
    with open(json_path, "w") as out_file:
        json.dump(diff_pos, out_file, indent = 6)

    with open(json_path_backup, "w") as out_file:
        json.dump(diff_pos, out_file, indent = 6)
        
    out_file.close()



def find_45_degree(start_angle,end_angle,num, elem="Pt_L"):


    ''' Absolute angles'''

    yield from bps.mov(dsth,start_angle)
    #num = np.float(num)
    #start_angle = np.float(start_angle)
    #end_angle = np.float(end_angle)
    step = (end_angle-start_angle)/num
    w_x = np.zeros(num+1)
    w_z = np.zeros(num+1)
    th = np.zeros(num+1)
    for i in range(num+1):
        yield from fly1d(dets1,dssx,-10,10,200,0.05)
        l,r,c=square_fit(-1,elem)
        plt.close()
        w_x[i] = r-l
        yield from fly1d(dets1,dssz,-10,10,200,0.05)
        l,r,c=square_fit(-1,elem)
        plt.close()
        w_z[i] = r-l
        th[i]=dsth.position
        yield from bps.sleep(1)
        yield from bps.movr(dsth,step)
    plt.figure()
    plt.plot(th,w_x,'r+',th,w_z,'g-')
    return th,w_x,w_z 





    

def VMS_in():

    print("Please wait...")
    
    #vms

    #caput("XF:03IDA-OP{VMS:1-Ax:Y}Mtr.VAL", -0.07) #mirrorY
    #caput("XF:03IDA-OP{VMS:1-Ax:P}Mtr.VAL", 3.06) #mirror picth
    caput("XF:03IDA-OP{VMS:1-Ax:YU}Mtr.VAL", -1.46) #upstream Y
    caput("XF:03IDA-OP{VMS:1-Ax:YD}Mtr.VAL",0.39) #downstream Y
    caput("XF:03IDA-OP{VMS:1-Ax:TX}Mtr.VAL", 0) #trans. X
    #caput("XF:03IDA-OP{VMS:1-Ax:PF}Mtr.VAL", 7) #pitch fine

    #bbpm
    caput("XF:03IDB-OP{Slt:SSA1-Ax:7}Mtr.VAL", -2.4)
    caput("XF:03IDB-OP{Slt:SSA1-Ax:8}Mtr.VAL",0.25)
    
    #cbpm
    caput("XF:03IDC-ES{BPM:7-Ax:Y}Mtr.VAL", 1.3)

    for i in tqdm.tqdm(range(30)):
        yield from bps.sleep(1)

    
    
    #move ssa2
    yield from bps.mov(ssa2.hgap, 2, ssa2.vgap,2, ssa2.hcen,0.098,ssa2.vcen, 0)
    
    #move ssa1
    yield from bps.mov(ssa1.hgap, 2.0, ssa1.vgap,2.0, ssa1.hcen,0.079,ssa1.vcen, 1.7)

    print(" Aligned to VMS")

def VMS_out():


    print("Please wait...")

    #vms

    #caput("XF:03IDA-OP{VMS:1-Ax:Y}Mtr.VAL", -1.0) #mirrorY
    #caput("XF:03IDA-OP{VMS:1-Ax:P}Mtr.VAL", -0.0082) #mirror picth
    caput("XF:03IDA-OP{VMS:1-Ax:YU}Mtr.VAL", -2.2) #upstream Y
    caput("XF:03IDA-OP{VMS:1-Ax:YD}Mtr.VAL",-2.6) #downstream Y
    #caput("XF:03IDA-OP{VMS:1-Ax:TX}Mtr.VAL", 0.084) #trans. X
    caput("XF:03IDA-OP{VMS:1-Ax:PF}Mtr.VAL", 0) #pitch fine

    #bbpm
    caput("XF:03IDB-OP{Slt:SSA1-Ax:7}Mtr.VAL", -0.16)
    caput("XF:03IDB-OP{Slt:SSA1-Ax:8}Mtr.VAL",0.065)
    
    #cbpm
    caput("XF:03IDC-ES{BPM:7-Ax:Y}Mtr.VAL", 0.4)

    
    for i in tqdm.tqdm(range(0, 30), desc ="Moving..."):
        yield from bps.sleep(1)
    

    #move ssa2
    yield from bps.mov(ssa2.hgap, 2, ssa2.vgap,2, ssa2.hcen,0,ssa2.vcen, -0.6765)
    
    #move ssa1
    yield from bps.mov(ssa1.hgap, 2.64, ssa1.vgap,2.5, ssa1.hcen,-0.092,ssa1.vcen, -0.78)

    print(" VMS out")


def feedback_auto_off(wait_time_sec = 0.5):
    
    beam_current = "SR:C03-BI{DCCT:1}I:Real-I"
    fe_xbpm_current = "SR:C03-BI{XBPM:1}Ampl:CurrTotal-I"
    fe_shutter_status = "XF:03ID-PPS{Sh:FE}Sts:Cls-Sts"
    ugap = "SR:C3-ID:G1{IVU20:1-Ax:Gap}-Mtr.RBV"

    b_feeback_x = "XF:03ID-BI{EM:BPM1}fast_pidX.FBON"
    b_feeback_y = "XF:03ID-BI{EM:BPM1}fast_pidY.FBON"

    while True:
        if caget(beam_current)<10 or caget(fe_xbpm_current)<10 or caget(fe_shutter_status)==1:

            if caget(b_feeback_x) == 1 or caget(b_feeback_y) == 1:
                caput(b_feeback_x,0)
                caput(b_feeback_y,0)
                logger.info(f"feedback was disabled by {os.getlogin()}")

            else:
                pass
            
        time.sleep(wait_time_sec)


def check_for_beam_dump(threshold = 5000):

    while (sclr2_ch2.get() < threshold):
        yield from bps.sleep(60)
        print (f"IC3 is lower than {threshold}, waiting...")


def find_edge_2D(scan_id, elem, left_flag=True):

    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    motors = db[scan_id].start['motors']
    x = np.array(df2[motors[0]])
    y = np.array(df2[motors[1]])
    #I0 = np.asfarray(df2.sclr1_ch4)
    I0 = np.asfarray(df2['sclr1_ch4'])
    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']
    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))
    l = np.linspace(y[0],y[-1],ny)
    s = xrf.sum(1)
    #if axis == 'x':
        #l = np.linspace(x[0],x[-1],nx)
        #s = xrf.sum(0)
    #else:
        #l = np.linspace(y[0],y[-1],ny)
        #s = xrf.sum(1)


	#plt.figure()
	#plt.plot(l,s)
	#plt.show()
	#sd = np.diff(s)
    sd = np.gradient(s)
    if left_flag:
        edge_loc1 = l[np.argmax(sd)]
    else:
        edge_loc1 = l[np.argmin(sd)]
    #plt.plot(l,sd)
	#plt.title('edge at '+np.str(edge_loc1))

    sd2 = np.diff(s)
    ll = l[:-1]
	#plt.plot(ll,sd2)
    if left_flag:
        edge_loc2 = ll[np.argmax(sd2)]
    else:
        edge_loc2 = ll[np.argmin(sd2)]
	#plt.xlabel('edge at '+np.str(edge_loc2))

	#edge_pos=find_edge(l,s,10)
	#pos = l[s == edge_pos]
	#pos = l[s == np.gradient(s).max()]
	#popt,pcov=curve_fit(erfunc1,l,s, p0=[edge_pos,0.05,0.5])
    return edge_loc1,edge_loc2







