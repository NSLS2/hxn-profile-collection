#Last Update 07/15/2021


import numpy as np
from datetime import datetime

user_folder = '/data/users/2021Q2/Huang_2021Q2/Cu2O' #This is to add scan details in .txt format to your folder

pre = np.linspace(8.96,8.975,4)
XANES1 = np.arange(8.976,9.025,0.001)
post = np.linspace(9.030,9.05,5)
#energies = np.concatenate([pre,XANES1,post])

#pre = np.linspace(8.96,8.975,4)
#XANES1 = np.arange(8.976,9.020,0.0005)

#calib_energies = np.concatenate([pre,XANES1])


energies = np.array([8.95,8.94,8.93, 8.92])

#energies = calib_energies


high_e = 9.05
low_e = 8.96
high_e_ugap = 6150
low_e_ugap = 6102

ugap_slope = (high_e_ugap - low_e_ugap)/(high_e-low_e)
ugap_list = high_e_ugap + (energies - high_e)*ugap_slope


high_e_crl= 10 #22 best theta is 18 but motor won't stay 
low_e_crl = 8
crl_slope = (high_e_crl - low_e_crl)/(high_e-low_e)
crl_list = high_e_crl + (energies - high_e)*crl_slope

zpz1_ref = 54.8
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (energies - high_e)*zpz1_slope

e_list = np.column_stack((energies,ugap_list,zpz1_list,crl_list))
#e_list = e_list[::-1]

def zp_list_xanes2d(e_list,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, 
                    xcen = 0,ycen = 0):

    realE_list = []
    scanid_list = []

    ic_0 = sclr2_ch2.get()

    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1)
    yield from bps.sleep(2)

    ic_3 =  sclr2_ch4.get()
    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0)
    zpssx_i = zpssx.position
    zpssy_i = zpssy.position

    #particle location, I keep the beam away (0,0) from it when so scans done
    #xcen = -4.7
    #ycen = 2.06


    for i in range (len(e_list)):

        yield from bps.sleep(1)

        caput('XF:03IDC-ES{Status}ScanRunning-I', 1)  #tuning the scanning pv on to dispable c bpms

        yield from bps.mov(e,e_list[i][0])
        yield from bps.sleep(1)
        yield from bps.mov(ugap, e_list[i][1])
        yield from bps.sleep(2)
        yield from mov_zpz1(e_list[i][2])
        yield from bps.sleep(1)
        yield from bps.mov(crl.p,e_list[i][3])
        yield from bps.sleep(1)

        caput('XF:03IDC-ES{Status}ScanRunning-I', 0)

        yield from bps.sleep(2)

        while (sclr2_ch2.get() < (0.1*ic_0)):
            yield from bps.sleep(60)
            print('IC3 is lower than 10000, waiting...')
        
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) #opening fast shutter

        yield from bps.sleep(2)
        
        

        if  sclr2_ch4.get() < (0.85*ic_3):
            yield from peak_bpm_y(-2,2,4)
            yield from peak_bpm_x(-15,15,5)
            yield from peak_bpm_y(-2,2,4)
            

        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) #closing fast shutter
        
        yield from bps.sleep(2)

        # move to particle location for alignemnt scan
        yield from bps.mov(zpssx, xcen)
        yield from bps.mov(zpssy, ycen)
        
        # do the alignemnt scan on the xanes elem after it excited , otherwise skip or use another element

        if i%2 == 0 and i<47:
            yield from fly1d(dets,zpssx,-10,10,100,0.02)
            xcen = return_line_center(-1,'Cu',0.6)
            yield from bps.mov(zpssx, xcen)
            yield from fly1d(dets,zpssy,-10,10 ,100,0.02)
            ycen = return_line_center(-1,'Cu',0.6)
            yield from bps.mov(zpssy, ycen)
       
        
        print(f'Current scan: {i+1}/{len(e_list)}')

        # do the fly2d scan

        if dets == dets_fs: #for fast xanes scan, no transmission (merlin) in the list

            yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, dead_time=0.001) #dead_time = 0.001 for 0.015 dwell

        else:

            yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
        yield from bps.sleep(1)
        
        # after scan done go to 0,0 to rest
        
        #yield from bps.mov(zpssx, 0)
        #yield from bps.mov(zpssy, 0)

        #ycen, xcen = return_center_of_mass_blurr(-1,'S') # some cases use 2D mass center for alignemnt
        #print(ycen,xcen)

        # get some scan details and add to the list of scan id and energy

        h = db[-1]
        last_sid = h.start['scan_id']
        scanid_list.append(last_sid)
        e_pos = e.position
        realE_list.append(e_pos)

        #insert_xrf_map_to_pdf(-1,'Fe')# plot data and add to pdf
        insert_xrf_map_to_pdf(-1,'Cu')# plot data and add to pdf
        sid_e_list = np.column_stack([scanid_list,realE_list])
        
        '''

        e_list_df['E Readback'][i] = e_pos #add real energy to the dataframe
        e_list_df['Scan ID'][i] = int(last_sid) #add scan id to the dataframe
        e_list_df['TimeStamp'][i] = pd.Timestamp.now()
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) #opening fast shutter
        yield from bps.sleep(2)

        e_list_df['IC0'][i] = sclr2_ch2.get() #Ic values are useful for calibration
        e_list_df['IC3'][i] = sclr2_ch4.get() #Ic values are useful for calibration   
        #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) #close fast shutter
        
        filename = f"HXN_nanoXANES_{e_list_df['Scan ID'][0]}_to_{e_list_df['Scan ID'][0]}.csv"
        e_list_df.to_csv(os.path.join(user_folder, filename))
        '''

        np.savetxt(os.path.join(user_folder, 'Xanes_elist_startsid_{}'.format(scanid_list[0])+'.txt'),sid_e_list,fmt = '%5f') # saves in the loop so, quitting a scan won't affect
    
    np.savetxt(os.path.join(user_folder, 'elist_{}'.format(scanid_list[0])+'_to'+'_{}'.format(scanid_list[-1])+'.txt'),sid_e_list[:,1],fmt = '%5f')
    

    save_page() #save the pdf

#zp_list_xanes2d(e_list,dets1, zpssx,-3,3,50,zpssy,-3,3,50,0.03)

