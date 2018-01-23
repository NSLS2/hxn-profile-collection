from scipy.optimize import curve_fit

def erfunc1(z,a,b,c):
    return c*(scipy.special.erf((z-a)/(b*np.sqrt(2.0)))+1.0)
def erfunc2(z,a,b,c):
    return c*(1.0-scipy.special.erf((z-a)/(b*np.sqrt(2.0))))
def erf_fit(sid,mot,elem,mon='sclr1_ch4'):

    h=db[sid]
    sid=h['start']['scan_id']
    df=db.get_table(h)
    xdata=df[mot]
    xdata=np.array(xdata,dtype=float)
    #x_mean=np.mean(xdata)
    #xdata=xdata-x_mean
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
        popt,pcov=curve_fit(erfunc1,xdata,ydata, p0=[edge_pos,0.05,0.5])
        fit_data=erfunc1(xdata,popt[0],popt[1],popt[2]);
    else:
        popt,pcov=curve_fit(erfunc2,xdata,ydata,p0=[edge_pos,0.05,0.5])
        fit_data=erfunc2(xdata,popt[0],popt[1],popt[2]);
    #print('a={} b={} c={}'.format(popt[0],popt[1],popt[2]))
    plt.plot(xdata,fit_data)
    plt.title('sid= %d edge = %.3f, FWHM = %.2f nm' % (sid,popt[0], popt[1]*2.3548*1000.0))
    return (popt[0],popt[1]*2.3548*1000.0)

def mll_z_alignment(z_start, z_end, z_num, mot, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_sz = smlld.sbz.position
    movr(smlld.sbz, z_start)
    for i in range(z_num + 1):
        if mot == 'dssy':
            RE(fly1d(dssy, start, end, num, acq_time))
        elif mot == 'dssx':
            RE(fly1d(dssx, start, end, num, acq_time))
        else:
            raise KeyError('mot has to be dssx or dssy')
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
        edge_pos,fwhm=erf_fit(-1,mot,elem,mon)
        fit_size[i]=fwhm
        z_pos[i]=smlld.sbz.position
        movr(smlld.sbz, z_step)
    mov(smlld.sbz, init_sz)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')

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
def hmll_z_alignment(z_start, z_end, z_num, mot, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_hz = hmll.hz.position
    movr(hmll.hz, z_start)
    for i in range(z_num + 1):
        if mot == 'dssx':
            RE(fly1d(dssx, start, end, num, acq_time))
        else:
            raise KeyError('mot has to be dssx')
        edge_pos,fwhm=erf_fit(-1,mot,elem,mon)
        fit_size[i]=fwhm
        z_pos[i]=hmll.hz.position
        movr(hmll.hz, z_step)
    mov(hmll.hz, init_hz)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')
def vmll_z_alignment(z_start, z_end, z_num, mot, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_vz = vmll.vz.position
    movr(vmll.vz, z_start)
    for i in range(z_num + 1):
        if mot == 'dssy':
            RE(fly1d(dssy, start, end, num, acq_time))
        else:
            raise KeyError('mot has to be dssx')
        edge_pos,fwhm=erf_fit(-1,mot,elem,mon)
        fit_size[i]=fwhm
        z_pos[i]=vmll.vz.position
        movr(vmll.vz, z_step)
    mov(vmll.vz, init_vz)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')

def zp_z_alignment(z_start, z_end, z_num, mot, start, end, num, acq_time, elem='Pt_L',mon='sclr1_ch4'):
    z_pos=np.zeros(z_num+1)
    fit_size=np.zeros(z_num+1)
    z_step = (z_end - z_start)/z_num
    init_zpz1 = zp.zpz1.position
    movr_zpz1(z_start)
    for i in range(z_num + 1):
        if mot == 'zpssx':
            RE(fly1d(zpssx, start, end, num, acq_time))
        elif mot == 'zpssy':
            RE(fly1d(zpssy, start, end, num, acq_time))
        else:
            raise KeyError('mot has to be zpssx or zpssy')
        edge_pos,fwhm=erf_fit(-1,mot,elem,mon)
        fit_size[i]=fwhm
        z_pos[i]=zp.zpz1.position
        movr_zpz1(z_step)
    movr_zpz1(-z_end)
    plt.figure()
    plt.plot(z_pos,fit_size,'bo')

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


def return_line_center(sid,elem='Cr'):
    h = db[sid]

    df2 = h.table()
    xrf = np.array(df2['Det2_' + elem]+df2['Det1_' + elem] + df2['Det3_' + elem])
    threshold = np.max(xrf)/10.0
    x_motor = h.start['motor']
    x = np.array(df2[x_motor])
    #print(x)
    #print(xrf)
    xrf[xrf<(np.max(xrf)*0.25)] = 0.
    xrf[xrf>=(np.max(xrf)*0.25)] = 1.
    mc = find_mass_center_1d(xrf,x)
    return mc



def zp_rot_alignment(a_start, a_end, a_num, start, end, num, acq_time, elem='Cr', move_flag=0):
    a_step = (a_end - a_start)/a_num
    x = np.zeros(a_num+1)
    y = np.zeros(a_num+1)
    orig_th = zps.zpsth.position
    for i in range(a_num+1):
        x[i] = a_start + i*a_step
        yield from bps.mov(zps.zpsth, x[i])
        if np.abs(x[i]) > 45:
            yield from fly1d(dets1,zpssz,start,end,num,acq_time)
            tmp = return_line_center(-1, elem)
            y[i] = tmp*np.sin(x[i]*np.pi/180.0)
        else:
            yield from fly1d(dets1,zpssx,start,end,num,acq_time)
            tmp = return_line_center(-1,elem)
            y[i] = tmp*np.cos(x[i]*np.pi/180.0)
        print('y=',y[i])
    y = -1*np.array(y)
    x = np.array(x)
    r0, dr, offset = rot_fit_2(x,y)
    yield from bps.mov(zps.zpsth, 0)
    dx = -dr*np.sin(offset*np.pi/180)/1000.0
    dz = -dr*np.cos(offset*np.pi/180)/1000.0

    print('dx=',dx,'   ', 'dz=',dz)

    if move_flag:
        yield from bps.movr(zps.smarx, dx)
        yield from bps.movr(zps.smarz, dz)


    return x,y


