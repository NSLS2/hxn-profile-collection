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

