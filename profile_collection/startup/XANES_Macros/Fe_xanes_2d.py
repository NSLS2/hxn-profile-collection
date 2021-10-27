#Last Update 07/15/2021


import numpy as np
from datetime import datetime

#osaz = 5000, zpz = 10.465 mm, crl#12, -5 to +4 (increasing E)

pre = np.linspace(7.08,7.11,11)
XANES1 = np.linspace(7.111,7.131,41)
XANES2 = np.linspace(7.132,7.148,33)
XANES3 = np.linspace(7.15,7.18,11)
post = np.linspace(7.182,7.232,6)

#energies = np.concatenate([pre,XANES1,XANES2,XANES3,post])

#energies = np.array([7.08,7.113,7.121,7.132,7.136,7.18])

energies = np.concatenate([pre,XANES1])

user_folder = '/data/users/2021Q2/Ajith_2021Q2/' #This is to add scan details in .txt format to your folder


high_e = 7.2
low_e = 7.1
high_e_ugap = 7675
low_e_ugap = 7585
ugap_slope = (high_e_ugap - low_e_ugap)/(high_e-low_e)
ugap_list = high_e_ugap + (energies - high_e)*ugap_slope


high_e_crl= 4
low_e_crl = -5
crl_slope = (high_e_crl - low_e_crl)/(high_e-low_e)
crl_list = high_e_crl + (energies - high_e)*crl_slope

zpz1_ref = 65.858-0.37
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (energies - high_e)*zpz1_slope

e_list = np.column_stack((energies,ugap_list,zpz1_list,crl_list))
e_list_df = pd.DataFrame()
e_list_df['E Targets'] = energies
e_list_df['E Readback'] = np.nan #add real energy to the dataframe
e_list_df['Scan ID'] = np.nan #add scan id to the dataframe
e_list_df['TimeStamp'] = pd.Timestamp.now()
e_list_df['IC3'] = sclr2_ch4.get() #Ic values are useful for calibration
e_list_df['IC0'] = sclr2_ch2.get() #Ic values are useful for calibration   

#e_list = e_list[::-1]

def zp_list_xanes2d(e_list,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, 
                    xcen = -4.7,ycen = 2.06):

    realE_list = []
    scanid_list = []

    ic_0 = sclr2_ch2.get()+20000

    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1)
    yield from bps.sleep(2)

    ic_3 =  sclr2_ch4.get()+10000
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
        #yield from mov_zpz1(e_list[i][2])
        yield from bps.sleep(1)
        #yield from bps.mov(crl.p,e_list[i][3])
        yield from bps.sleep(1)

        caput('XF:03IDC-ES{Status}ScanRunning-I', 0)

        yield from bps.sleep(2)

        while (sclr2_ch2.get() < (0.1*ic_0)):
            yield from bps.sleep(60)
            print('IC3 is lower than 10000, waiting...')
        '''
        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) #opening fast shutter

        yield from bps.sleep(2)
        
        

        if  sclr2_ch4.get() < (0.9*ic_3):
            yield from peak_bpm_y(-2,2,4)
            yield from peak_bpm_x(-15,15,5)
            yield from peak_bpm_y(-2,2,4)
            
       

        caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) #closing fast shutter
        
        yield from bps.sleep(2)

        # move to particle location for alignemnt scan
        yield from bps.mov(zpssx, xcen)
        yield from bps.mov(zpssy, ycen)
        
        # do the alignemnt scan on the xanes elem after it excited , otherwise skip or use another element

        if e_list[i][0]<7.120:
            pass
            
            yield from fly1d(dets,zpssx,-1,1,100,0.1)
            xcen = return_line_center(-1,'Cl',0.6)
            yield from bps.mov(zpssx, xcen)
            yield from fly1d(dets,zpssy,-1,1 ,100,0.1)
            ycen = return_line_center(-1,'Cl',0.6)
            yield from bps.mov(zpssy, ycen)
            

        else:
            yield from fly1d(dets,zpssx,-2,2,100,0.03)
            xcen = return_line_center(-1,'Fe',0.6)
            yield from bps.mov(zpssx, xcen)
            yield from fly1d(dets,zpssy,-2,2 ,100,0.03)
            ycen = return_line_center(-1,'Fe',0.6)
            yield from bps.mov(zpssy, ycen)
        '''

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
        insert_xrf_map_to_pdf(-1,'Cl')# plot data and add to pdf
        sid_e_list = np.column_stack([scanid_list,realE_list])
        

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

        np.savetxt(os.path.join(user_folder, 'Xanes_elist_startsid_{}'.format(scanid_list[0])+'.txt'),sid_e_list,fmt = '%5f') # saves in the loop so, quitting a scan won't affect
    
    np.savetxt(os.path.join(user_folder, 'elist_{}'.format(scanid_list[0])+'_to'+'_{}'.format(scanid_list[-1])+'.txt'),sid_e_list[:,1],fmt = '%5f')
    

    save_page() #save the pdf

#zp_list_xanes2d(e_list,dets1, zpssx,-3,3,50,zpssy,-3,3,50,0.03)

