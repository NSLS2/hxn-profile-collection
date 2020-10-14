
'''
pre = np.linspace(9.634,9.655,8)
XANES = np.linspace(9.656,9.700,45)
post = np.linspace(9.705,9.735,6)

'''
#pre = np.linspace(9.645,9.655,11)
XANES = np.arange(9.645,9.7,.0005)
#XANES2 = np.linspace(9.685,9.700,16)
#post = np.linspace(9.705,9.755,11)


energies = XANES
#energies = np.concatenate([pre,XANES])
#energies = np.concatenate([pre,XANES1,XANES2])

#energies = np.asarray([7.093,7.13,7.173]) #for test only
high_e = 9.7
low_e = 9.645
high_e_ugap = 6475
low_e_ugap = 6448
ugap_slope = (high_e_ugap - low_e_ugap)/(high_e-low_e)
ugap_list = high_e_ugap + (energies - high_e)*ugap_slope


high_e_crl= 8  #22+5
low_e_crl = 2
crl_slope = (high_e_crl - low_e_crl)/(high_e-low_e)
crl_list = high_e_crl + (energies - high_e)*crl_slope

zpz1_ref = -18.23
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (energies - high_e)*zpz1_slope

e_list = np.column_stack((energies,ugap_list,zpz1_list,crl_list))
#e_list = e_list[::-1]


def zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
    
    realE_list = []
    scanid_list = []

    ic_0 = sclr2_ch4.get()

    for i in range (len(e_list)):
        '''
        if (sclr2_ch4.get() < (0.1*ic_0)):
            i = i-1

        while (sclr2_ch4.get() < (0.1*ic_0)):
            yield from bps.sleep(60)
            print('IC3 is lower than 10000, waiting...')
        '''

        yield from bps.mov(e,e_list[i][0])
        yield from bps.sleep(1)
        yield from bps.mov(ugap, e_list[i][1])
        yield from bps.sleep(1)
        #yield from mov_zpz1(e_list[i][2])
        yield from bps.sleep(1)
        yield from bps.mov(crl.p,e_list[i][3])
        yield from bps.sleep(1)

        '''

        if (sclr2_ch4.get() < (0.9*ic_0)) or i % 10 == 0:
       

            yield from peak_bpm_y(-10,10,10)
            yield from peak_bpm_x(-30,30,10)
            yield from peak_bpm_y(-10,10,10)

        '''
        print(e.position)
        yield from fly2d(dets1, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        yield from bps.sleep(1)

        '''
        xcen, ycen  = return_center_of_mass_blurr(-1,'S',5,1)
        if abs(xcen) < abs(x_s) and abs(ycen) < abs(y_s):

            yield from bps.mov(zpssx,xcen)
            yield from bps.mov(zpssy,ycen)  

        else:

            yield from bps.mov(zpssx,0)
            yield from bps.mov(zpssy,0)

            yield from bps.movr(smarx, (xcen*0.001))
            yield from bps.movr(smary,(ycen*0.001))   
        
        '''
        h = db[-1]
        last_sid = h.start['scan_id']
        scanid_list.append(last_sid)
        e_pos = e.position
        realE_list.append(e_pos)


        sid_e_list = np.column_stack([scanid_list,realE_list])
        user_folder = '/data/users/2019Q3/Ajith_2109Q3/Myneni2019Q3/'
        np.savetxt(os.path.join(user_folder, 'Calib_Xanes_elist_startsid_{}'.format(scanid_list[0])+'.txt'),sid_e_list,fmt = '%5f')
        print('sid_list saved')
        yield from bps.sleep(10)

#zp_list_xanes2d(e_list,zpssx,-3,3,50,zpssy,-3,3,50,0.03)  
#d2scan(dets2,100,e,0,0.055,ugap,0,27,0.5) 

#


