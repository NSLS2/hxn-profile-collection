
#osaz = 2500, zpz = 8 mm, angle , , crl#9+6, No CRL change (increasing E), Si(111)

user_folder = '/data/users/2021Q2/Tyson_2021Q2/'

pre = np.linspace(6.520,6.534,8)
XANES1 = np.linspace(6.5345,6.560,72)#50 for 0.5eV, 33 for 0.75 
#XANES2 = np.linspace(6.561,6.570,10)#changes to 10 fro calibration, 20 otherwise
post = np.linspace(6.561,6.6,40)

Mn_energies = np.concatenate([pre,XANES1,post]) #81 points
#Mn_energies = np.concatenate([XANES2])
#Mn_energies=np.asarray([6.525])

ugap_ref = 7095
e_ref = 6.54
ugap_slope = (7130-7085)/0.055
ugap_list = ugap_ref + (Mn_energies - e_ref)*ugap_slope


crl_ref = -12
crl_slope = (0)/0.066
crl_list = crl_ref + (Mn_energies - e_ref)*crl_slope

zpz1_ref = -6.5628 #zpz = 10
zpz1_slope = -5.91
zpz1_list = zpz1_ref + (Mn_energies - e_ref)*zpz1_slope

e_list = np.column_stack((Mn_energies,ugap_list,zpz1_list,crl_list))
np.savetxt(user_folder+'Mn_e_list.txt', e_list,fmt='%5.5f')

#realE_list = []
#scanid_list = []
def zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
    num_pts, num_mot = np.shape(e_list)
    
    realE_list = []
    scanid_list = []

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
        #yield from bps.mov(crl.p,crl_angle)
        #yield from bps.sleep(2)

        yield from fly2d(dets_fs, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t,dead_time=0.001,return_speed=100)

        #yield from dmesh(dets1, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        #yield from bps.sleep(2)
        #yield from bps.movr(zpssx, -0.002)
        h = db[-1]
        if i == 0:
            first_sid = h.start['scan_id']
        last_sid = h.start['scan_id']
        scanid_list.append(last_sid)
        e_pos = e.position
        realE_list.append(e_pos)
        plot2dfly(-1,'Au_M','sclr1_ch4')
        insertFig(note = 'energy={}'.format(e_pos), title ='LuMnO3 bottom')
        plt.close()
        plot2dfly(-1,'Mn','sclr1_ch4')
        insertFig(note = 'energy={}'.format(e_pos), title ='LuMnO3 bottom')
        plt.close()
        #merlin1.unstage()
        xspress3.unstage()

        while (sclr2_ch2.get() < 100000):
            yield from bps.sleep(60)
            print('IC1 is lower than 200000, waiting...')

    sid_e_list = np.column_stack([scanid_list,realE_list])
    np.savetxt(os.path.join(user_folder, 'elist_{}_{}'.format(first_sid,last_sid)+'.txt'),realE_list, fmt = '%5f')

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
    yield from zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)




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


