

#osaz = 1500, zpz = 5 mm, angle 15, smaz = 0.44, crl#12, -5 to +4 (increasing E)

pre = np.linspace(8.310,8.322,5)
XANES1 = np.linspace(8.322,8.349,28)
XANES2 = np.linspace(8.350,8.360,21)
XANES3 = np.linspace(8.361,8.380,20)
post = np.linspace(8.381,8.390,4)

energies = np.concatenate([pre,XANES1,XANES2,XANES3,post])

#Ni_energies = np.asarray([8.310,8.350,8.390]) #for test only

ugap_ref = 5802
e_ref = 8.355
ugap_slope = 50/0.1
ugap_list = ugap_ref + (Fe_energies - e_ref)*ugap_slope


crl_ref = 4
crl_slope = (6)/0.107
crl_list = crl_ref + (Fe_energies - e_ref)*crl_slope

zpz1_ref = -5.73
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (Fe_energies - e_ref)*zpz1_slope

e_list = np.column_stack((Fe_energies,ugap_list,zpz1_list,crl_list))

def zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
    num_pts, num_mot = np.shape(e_list)
    for i in range(num_pts):
        energy = e_list[i][0]
        gap_sz = e_list[i][1]
        zpz1_pos = e_list[i][2]
        crl_angle = e_list[i][3]
        yield from bps.mov(e,energy)
        yield from bps.sleep(2)
        yield from bps.mov(ugap, gap_sz)
        yield from bps.sleep(2)
        yield from mov_zpz1(zpz1_pos)
        yield from bps.sleep(2)
        yield from bps.mov(crl.p,crl_angle)
        yield from bps.sleep(2)

        #yield from fly2d(dets1, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        yield from dmesh(dets1, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        yield from bps.sleep(2)

        #yield from bps.movr(zps.smary,(0.5*0.001/(len(e_list)))) #to adjust drift,temp
        #yield from bps.movr(zps.smarx,(0.35*0.001/len(e_list))) #to adjust drift_tem

        while (sclr2_ch4.get() < 5000):
            yield from bps.sleep(60)
            print('IC3 is lower than 5000, waiting...')


def return_center_of_mass(scan_id = -1, elem = 'Cr'):
    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    x = np.asarray(df2.zpssx)
    y = np.asarray(df2.zpssy)
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    b = ndimage.measurements.center_of_mass(xrf)

    iy = np.int(np.round(b[0]))
    ix = np.int(np.round(b[1]))
    i_max = ix + iy * nx

    x_cen = x[i_max]
    y_cen = y[i_max]
    return (y_cen, x_cen)


def zp_xanes2d(param_file, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
    e_list = np.loadtxt(param_file)
    zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)




    '''yield from bps.mov(zpssy,0)
       
		yield from fly1d(dets1,zpssy, -14, 14, 100, 0.1)
        yield from bps.sleep(3)
        yc = return_line_center(-1,'Fe',0.2)
        if abs(yc)<1:
            yield from bps.movr(zps.smary,yc/1000)
        plt.close()
        
        yield from bps.mov(zpssx,0)
        yield from fly1d(dets1, zpssx, -14, 14, 100, 0.1)
        yield from bps.sleep(3)
        xc = return_line_center(-1,'Fe',0.2)
        if abs(xc)<2.5:
            if not xc == NaN:
                yield from bps.movr(zps.smarx,xc/1000)
        plt.close()'''

    '''h = db[-1]
        sid = h.start['scan_id']
        plot2dfly(-1,'Fe')
        insertFig(note ='e = {}'.format(energy),title = ' ')
        plt.close()
        yield from bps.sleep(2)
        
    save_page() '''


