

def generate_e_list(elem = 'Fe'):

    if elem = 'Fe':

        pre = np.linspace(7.087,7.102,6)
	XANES1 = np.linspace(7.103,7.129,27)
	XANES2 = np.linspace(7.130,7.140,11)
	XANES3 = np.linspace(7.141,7.160,10)
	post = np.linspace(7.161,7.173,5)
	energies = np.concatenate([pre,XANES1,XANES2,XANES3,post])

	high_e = 7.2
	low_e = 7.1
	high_e_ugap = 7678 # harm = 3, pitch 1.01709
	low_e_ugap = 7588


	high_e_crl= 5  crl#12
	low_e_crl = -2

	zpz1_ref = -9.285 #zpz = 6 mm
	zpz1_slope = -5.9 #dr = 30nm, D = 240um


    if elem = 'Zn':
        Pre = np.linspace(9.645,9.660,6)
        XANES = np.linspace(9.661,9.700,40)
        Post = np.linspace(9.703,9.718,6)
        energies = np.concatenate([pre,XANES,post])

        high_e = 9.7
        low_e = 9.645
        high_e_ugap = 6480
        low_e_ugap = 6455

        high_e_crl= 5
        low_e_crl = -2

        zpz1_ref = -13.61
        zpz1_slope = -5.9





    if elem = 'Ni':
    if elem = 'Mn':
    if elem = 'As':

    ugap_slope = (high_e_ugap - low_e_ugap)/(high_e-low_e)
    ugap_list = high_e_ugap + (energies - high_e)*ugap_slope

    crl_slope = (high_e_crl - low_e_crl)/(high_e-low_e)
    crl_list = high_e_crl + (energies - high_e)*crl_slope

    zpz1_list = zpz1_ref + (energies - high_e)*zpz1_slope

    e_list = np.column_stack((energies,ugap_list,zpz1_list,crl_list))

    return e_list

def zp_list_xanes2d(elem ='Fe',mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
   
    e_list = generate_e_list(elem = elem)
    
    realE_list = []
    scanid_list = []

    ic_0 = sclr2_ch4.get()
    zpssx_i = zpssx.position 
    zpssy_i = zpssy.position

    for i in range (len(e_list)):


        yield from bps.mov(e,e_list[i][0])
        yield from bps.sleep(1)
        yield from bps.mov(ugap, e_list[i][1])
        yield from bps.sleep(1)
        yield from mov_zpz1(e_list[i][2])
        yield from bps.sleep(1)
        yield from bps.mov(crl.p,e_list[i][3])
        yield from bps.sleep(1)

        while (sclr2_ch4.get() < (0.1*ic_0)):
            yield from bps.sleep(60)
            print('IC3 is lower than 10000, waiting...')

        if (sclr2_ch4.get() < (0.9*ic_0)):
        #if i in np.arange(0,80,10):

            yield from peak_bpm_y(-10,10,10)
            yield from peak_bpm_x(-30,30,10)
            yield from peak_bpm_y(-10,10,10)

        yield from bps.mov(zpssx,zpssx_i)
        yield from bps.mov(zpssy,zpssy_i)

        xcen, ycen  = return_center_of_mass_blurr(-1,'Fe',10,1)
        if abs(xcen)<1 and abs(ycen)<1:
            yield from bps.mov(zpssx, xcen)
            yield from bps.mov(zpssy,ycen)
        
        yield from fly2d(dets1, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        yield from bps.sleep(1)

        h = db[-1]
        last_sid = h.start['scan_id']
        scanid_list.append(last_sid)
        e_pos = e.position
        realE_list.append(e_pos)
        
        sid_e_list = np.column_stack([scanid_list,realE_list])
        user_folder = '/data/users/2019Q3/Ajith_2109Q3/'
        np.savetxt(os.path.join(user_folder, 'Xanes_elist_startsid_{}'.format(scanid_list[0])+'.txt'),sid_e_list,fmt = '%5f')
    '''
    yield from bps.mov(e,e_list[0][0])
    yield from bps.mov(ugap, e_list[0][1])
    yield from mov_zpz1(e_list[0][2])
    yield from bps.mov(crl.p,e_list[0][3])
    '''

#zp_list_xanes2d(e_list,zpssx,-3,3,50,zpssy,-3,3,50,0.03)  

