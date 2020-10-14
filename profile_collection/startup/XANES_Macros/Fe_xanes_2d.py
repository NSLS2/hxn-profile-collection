

#Last update:Mar2020
#osaz = 1500, zpz = 5 mm, crl#12, -5 to +4 (increasing E)

pre = np.linspace(7.087,7.102,6)
XANES1 = np.linspace(7.103,7.129,27)
XANES2 = np.linspace(7.130,7.140,11)
XANES3 = np.linspace(7.141,7.160,10)
post = np.linspace(7.161,7.176,6)

energies = np.concatenate([pre,XANES1,XANES2,XANES3,post])

high_e = 7.2
low_e = 7.1
high_e_ugap = 7678
low_e_ugap = 7588
ugap_slope = (high_e_ugap - low_e_ugap)/(high_e-low_e)
ugap_list = high_e_ugap + (energies - high_e)*ugap_slope


high_e_crl= 4
low_e_crl = -5
crl_slope = (high_e_crl - low_e_crl)/(high_e-low_e)
crl_list = high_e_crl + (energies - high_e)*crl_slope

zpz1_ref = -16.1656
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (energies - high_e)*zpz1_slope

e_list = np.column_stack((energies,ugap_list,zpz1_list,crl_list))

#e_list = e_list[::-1]

def zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
    
    realE_list = []
    scanid_list = []

    ic_0 = 533800.0
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
            yield from peak_bpm_y(-5,5,5)
            yield from peak_bpm_x(-15,15,5)
            yield from peak_bpm_y(-5,5,5)
        '''
        yield from fly1d(dets1,zpssx,-2, 2,50,0.1)
        xcen = return_line_center(-1,'P',0.75)
        if abs(xcen)<0.5:
            yield from bps.mov(zpssx, xcen)
        '''
        yield from fly2d(dets1, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        yield from bps.sleep(1)

        h = db[-1]
        last_sid = h.start['scan_id']
        scanid_list.append(last_sid)
        e_pos = e.position
        realE_list.append(e_pos)
        
        insert_xrf_map_to_pdf('Fe', 'sclr1_ch4') 
        insert_xrf_map_to_pdf('P', 'sclr1_ch4')
        
        
        sid_e_list = np.column_stack([scanid_list,realE_list])
        user_folder = '/data/users/2020Q1/Ajith_2020Q1/'
        np.savetxt(os.path.join(user_folder, 'Xanes_elist_startsid_{}'.format(scanid_list[0])+'.txt'),sid_e_list,fmt = '%5f')
    save_page()

#zp_list_xanes2d(e_list,zpssx,-3,3,50,zpssy,-3,3,50,0.03)  

