from scipy.optimize import curve_fit
import scipy
import pickle
import numpy as np

def erfunc3(z,a,b,c,d,e):
    return d+e*z+c*(scipy.special.erf((z-a)/(b*np.sqrt(2.0)))+1.0)
def erfunc4(z,a,b,c,d,e):
    return d+e*z+c*(1.0-scipy.special.erf((z-a)/(b*np.sqrt(2.0))))

def erfunc1(z,a,b,c):
    return c*(scipy.special.erf((z-a)/(b*np.sqrt(2.0)))+1.0)
def erfunc2(z,a,b,c):
    return c*(1.0-scipy.special.erf((z-a)/(b*np.sqrt(2.0))))
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
    y_min=np.min(ydata)
    y_max=np.max(ydata)
    ydata=(ydata-y_min)/y_max
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
    plt.title('sid= %d edge = %.3f, FWHM = %.2f nm' % (sid,popt[0], popt[1]*2.3548*1000.0))
    return (popt[0],popt[1]*2.3548*1000.0)



def square_fit(sid,elem,mon='sclr1_ch4',linear_flag=True):

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
    ydata=(ydata-y_min)/y_max
    plt.figure()
    plt.plot(xdata,ydata,'bo')
    edge_pos_1, edge_pos_2 = find_double_edge(xdata,ydata,10)
    popt,pcov=curve_fit(squarefunc,xdata,ydata,p0=[0.5,edge_pos_1,0.1,edge_pos_2,0.1])
    fit_data=squarefunc(xdata,popt[0],popt[1],popt[2],popt[3],popt[4]);

    #print('a={} b={} c={}'.format(popt[0],popt[1],popt[2]))
    plt.plot(xdata,fit_data)
    plt.title('sid= %d cen = %.3f e1 = %.3f e2 = %.3f ' % (sid,(popt[1]+popt[3])*0.5, popt[1],popt[3]))
    plt.xlabel(mots[0])
    return (popt[1],popt[3],(popt[1]+popt[3])*0.5)
    #return(xdata, ydata, fit_data)



def data_erf_fit(xdata,ydata,linear_flag=True):

    xdata=np.array(xdata,dtype=float)
    ydata=np.array(ydata,dtype=float)
    y_min=np.min(ydata)
    y_max=np.max(ydata)
    ydata=(ydata-y_min)/y_max
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

    #print('a={} b={} c={}'.format(popt[0],popt[1],popt[2]))
    plt.plot(xdata,fit_data)
    plt.title('sid= %d edge = %.3f, FWHM = %.2f nm' % (sid,popt[0], popt[1]*2.3548*1000.0))
    return (popt[0],popt[1]*2.3548*1000.0)



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
        edge_pos,fwhm=erf_fit(-1,elem,mon)
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
    index=scipy.argmin(zdata)
    index=index+j
    return xdata[index]

def find_double_edge(xdata, ydata, size):
    edge_1 = find_edge(xdata, ydata, size)
    index = scipy.argmax(ydata)
    cen = xdata[index]
    if cen > edge_1:
        edge_2 = (cen-edge_1) + cen
        return(edge_1,edge_2)
    else:
        edge_2 = cen - (edge_1 - cen)
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
        yield from fly1d(dets1, mot, start, end, num, acq_time)
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

def zp_z_alignment(z_start, z_end, z_num, mot, start, end, num, acq_time, elem=' ',mon='sclr1_ch4'):
    print('moves the zone plate incrementally and find the focus with a linescan at each position')
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_sz = zp.zpz1.position

    yield from movr_zpz1(z_start)
    for i in range(z_num + 1):

        yield from fly1d(dets1, mot, start, end, num, acq_time)
        edge_pos,fwhm=erf_fit(-1,elem,mon,linear_flag=False)
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
    pix = 74.8
    R = 2.315e5
    th1 = 0.7617
    phi1 = 3.0366
    th2 = 0.1796
    phi2 = 2.5335
    phi3 = -0.1246
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


def return_line_center(sid,elem='Cr',threshold=0.2):
    h = db[sid]

    df2 = h.table()
    xrf = np.array(df2['Det2_' + elem]+df2['Det1_' + elem] + df2['Det3_' + elem])
    #threshold = np.max(xrf)/10.0
    x_motor = h.start['motors']
    x = np.array(df2[x_motor[0]])

    #xrf = xrf * -1
    #xrf = xrf - np.min(xrf)

    #print(x)
    #print(xrf)
    xrf[xrf<(np.max(xrf)*threshold)] = 0.
    #index = np.where(xrf == 0.)
    #xrf[:index[0][0]] = 0.
    #xrf[xrf>=(np.max(xrf)*threshold)] = 1.
    mc = find_mass_center_1d(xrf,x)
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
    peak_index = np.where(xrf_d == np.max(xrf_d))
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

def zp_rot_alignment(a_start, a_end, a_num, start, end, num, acq_time, elem='Pt_L', move_flag=0):
    a_step = (a_end - a_start)/a_num
    x = np.zeros(a_num+1)
    y = np.zeros(a_num+1)
    orig_th = zps.zpsth.position
    for i in range(a_num+1):
        x[i] = a_start + i*a_step
        yield from bps.mov(zps.zpsth, x[i])
        if np.abs(x[i]) > 45:
            yield from fly1d(dets1,zpssz,start,end,num,acq_time)
            tmp = return_line_center(-1, elem=elem,threshold=0.5)
            #edge,fwhm = erf_fit(-1,elem = elem)
            y[i] = tmp*np.sin(x[i]*np.pi/180.0)
        else:
            yield from fly1d(dets1,zpssx,start,end,num,acq_time)
            tmp = return_line_center(-1,elem=elem,threshold=0.5)
            #edge,fwhm = erf_fit(-1,elem = elem)
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
    a_step = (a_end - a_start)/a_num
    x = np.zeros(a_num+1)
    y = np.zeros(a_num+1)
    v = np.zeros(a_num+1)
    orig_th = smlld.dsth.position
    for i in range(a_num+1):
        x[i] = a_start + i*a_step
        yield from bps.mov(smlld.dsth, x[i])
        if np.abs(x[i]) > 45:
            yield from fly1d(dets1,dssz,start,end,num,acq_time)
            tmp = return_line_center(-1, elem=elem)
            y[i] = tmp*np.sin(x[i]*np.pi/180.0)
        else:
            yield from fly1d(dets1,dssx,start,end,num,acq_time)
            tmp = return_line_center(-1,elem=elem)
            y[i] = -tmp*np.cos(x[i]*np.pi/180.0)
        #yield from fly1d(dets1,dssy,start,end,num,acq_time)
        #tmp = return_line_center(-1, elem=elem)
        #v[i] = tmp
        #print('h_cen= ',y[i],'v_cen = ',v[i])
        plot(-1,elem,'sclr1_ch4')
        #insertFig(note='dsth = {}'.format(check_baseline(-1,'dsth')))
        plt.close()

    y = -1*np.array(y)
    x = np.array(x)
    r0, dr, offset = rot_fit_2(x,y)
    #insertFig(note='dsth: {} {}'.format(a_start,a_end))
    yield from bps.mov(smlld.dsth, 0)
    dx = -dr*np.sin(offset*np.pi/180)
    dz = -dr*np.cos(offset*np.pi/180)

    print('dx=',dx,'   ', 'dz=',dz)

    if move_flag:
        yield from bps.movr(smlld.dsx, dx)
        yield from bps.movr(smlld.dsz, dz)

    #plt.figure()
    #plt.plot(x,v)

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
        return(bl[name].data[0])


def check_info(sid):
    h = db[sid]
    sid = h.start['scan_id']
    scan_time = datetime.fromtimestamp(h.start['time'])
    scan_uid = h.start['uid']
    scan_type = h.start['plan_name']
    scan_motors = h.start['motors']
    num_motors = len(scan_motors)
    det_list = h.start['detectors']
    exp_time = h.start['exposure_time']
    print('sid = {}'.format(sid), 'uid = ', scan_uid, scan_time)
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


