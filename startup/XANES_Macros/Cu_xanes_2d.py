
import numpy as np
from datetime import datetime



pre = np.linspace(8.96,8.975,4)
XANES1 = np.arange(8.976,9.025,0.001)
post = np.linspace(9.030,9.05,5)

#calib = np.arange(9.650,9.68,0.0005)

energies = np.concatenate([pre,XANES1,post])

#energies = np.array([7.08,7.113,7.121,7.132,7.136,7.18])

#energies = calib

user_folder = '/data/users/2021Q2/Ajith_2021Q2/PM_LA_025' #This is to add scan details in .txt format to your folder


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

zpz1_ref = -20.5283
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (energies - high_e)*zpz1_slope

e_list = np.column_stack((energies,ugap_list,zpz1_list,crl_list))
#e_list = e_list[::-1]


def zp_list_xanes2d(e_list,dets,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):

    realE_list = []
    scanid_list = []

    ic_3 = sclr2_ch4.get()+10000
    ic_0 = sclr2_ch2.get()

    caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1)
    yield from bps.sleep(3)

    ic_3 =  sclr2_ch4.get()
    #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0)
    zpssx_i = zpssx.position
    zpssy_i = zpssy.position

    #xcen = 0
    #ycen = 0


    for i in range (len(e_list)):

        yield from bps.sleep(1)

        caput('XF:03IDC-ES{Status}ScanRunning-I', 1)  #tuning the scanning pv on to dispable c bpms

        yield from bps.mov(e,e_list[i][0])
        yield from bps.sleep(1)
        yield from bps.mov(ugap, e_list[i][1])
        yield from bps.sleep(5)
        yield from mov_zpz1(e_list[i][2])
        yield from bps.sleep(3)
        yield from bps.mov(crl.p,e_list[i][3])
        yield from bps.sleep(3)

        caput('XF:03IDC-ES{Status}ScanRunning-I', 0)

        #yield from bps.sleep(3)

        while (sclr2_ch2.get() < (0.1*ic_0)):
            yield from bps.sleep(60)
            print('IC3 is lower than 10000, waiting...')

        #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1) #opening fast shutter

        #yield from bps.sleep(5)

        if  sclr2_ch4.get() < (0.90*ic_3):
            yield from peak_bpm_y(-2,2,4)
            yield from peak_bpm_x(-10,10,5)
            yield from peak_bpm_y(-2,2,4)

        #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0) #closing fast shutter

        yield from bps.sleep(2)

        if i>20:
            yield from fly1d(dets,zpssx,-1,1,50,0.25)
            xcen = return_line_center(-1,'Cu',0.7)
            yield from bps.mov(zpssx, xcen)
            yield from fly1d(dets,zpssy,-1,1 ,50,0.25)
            ycen = return_line_center(-1,'Cu',0.7)
            yield from bps.mov(zpssy, ycen)
        else:
            yield from fly1d(dets,zpssx,-1,1,50,0.25)
            xcen = return_line_center(-1,'Ca',0.7)
            yield from bps.mov(zpssx, xcen)
            yield from fly1d(dets,zpssy,-1,1 ,50,0.25)
            ycen = return_line_center(-1,'Ca',0.7)
            yield from bps.mov(zpssy, ycen)



        print(f'Current scan: {i+1}/{len(e_list)}')

        if dets == dets_fs:

            yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t, dead_time=0.001) #dead_time = 0.001 for 0.015 dwell

        else:

            yield from fly2d(dets, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
            yield from bps.sleep(1)

        #ycen, xcen = return_center_of_mass_blurr(-1,'S')

        #print(ycen,xcen)

        h = db[-1]
        last_sid = h.start['scan_id']
        scanid_list.append(last_sid)
        e_pos = e.position
        realE_list.append(e_pos)

        insert_xrf_map_to_pdf(-1,'Cu')
        insert_xrf_map_to_pdf(-1,'Fe')
        sid_e_list = np.column_stack([scanid_list,realE_list])

        np.savetxt(os.path.join(user_folder, 'Xanes_elist_startsid_{}'.format(scanid_list[0])+'.txt'),sid_e_list,fmt = '%5f')
    np.savetxt(os.path.join(user_folder, 'elist_{}'.format(scanid_list[0])+'_to'+'_{}'.format(scanid_list[-1])+'.txt'),sid_e_list[:,1],fmt = '%5f')
    save_page()

#zp_list_xanes2d(e_list,dets1,zpssx,-3,3,50,zpssy,-3,3,50,0.03)

