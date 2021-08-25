

#Hafnium L_I edge
#ZP #1, 244 um dia, 30 nm outmost, crl# 22, 8, 3

pre = np.linspace(11.24,11.34,51)
#XANES1 = np.linspace(9.552,9.580,29)
#post = np.linspace(9.582,9.640,30)

#energies = np.concatenate([pre,XANES1,post])
energies = pre

high_e = 11.340
low_e = 11.240
high_e_ugap = 7335
low_e_ugap = 7280
ugap_slope = (high_e_ugap - low_e_ugap)/(high_e-low_e)
ugap_list = high_e_ugap + (energies - high_e)*ugap_slope


high_e_crl= 10
low_e_crl = 8
crl_slope = (high_e_crl - low_e_crl)/(high_e-low_e)
crl_list = high_e_crl + (energies - high_e)*crl_slope

zpz1_ref = -34.186
energy_ref = 11.34
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (energies - energy_ref)*zpz1_slope

e_list = np.column_stack((energies,ugap_list,zpz1_list,crl_list))

#e_list = e_list[::-1]

def zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
    
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
       
        if np.abs(e_list[i][3]-crl.p.position)>0.2:
            yield from bps.mov(crl.p,e_list[i][3])
        yield from bps.sleep(1)

        while (sclr2_ch4.get() < (0.1*ic_0)):
            yield from bps.sleep(60)
            print('IC3 is lower than 10%, waiting...')

        if (sclr2_ch4.get() < (0.8*ic_0)):
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
        
        insert_xrf_map_to_pdf('Hf_L', 'sclr1_ch4') 
        insert_xrf_map_to_pdf('Ti', 'sclr1_ch4')
        
        
        sid_e_list = np.column_stack([scanid_list,realE_list])
        user_folder = '/data/users/2020Q3/Lavoie_2020Q3/'
        np.savetxt(os.path.join(user_folder, 'Xanes_elist_startsid_{}'.format(scanid_list[0])+'.txt'),sid_e_list,fmt = '%5f')
    save_page()

#zp_list_xanes2d(e_list,zpssx,-3,3,50,zpssy,-3,3,50,0.03)  

