
PreAs = np.linspace(11845,11860,6)
As_XANES = np.linspace(11861,11885,49)
PostAs = np.linspace(11886,11901,6)

As_energies = np.append(PreAs,As_XANES)
As_energies = (np.append(As_energies,PostAs))/1000

#np.savetxt('As_e_list_61pts.txt',As_energies, fmt='%f')

ugap_ref = 7680
e_ref = 7.2
ugap_slope = (7680 - 7590)/0.1
ugap_list = ugap_ref + (As_energies - e_ref)*ugap_slope

crl_ref = 5
crl_slope = (5 + 5)/0.1
crl_list = crl_ref + (As_energies - e_ref)*crl_slope

zpz1_ref = -4.15
zpz1_slope = -5.9
zpz1_list = zpz1_ref + (As_energies - e_ref)*zpz1_slope


e_list = np.vstack((As_energies,ugap_list,zpz1_list,crl_list))
e_list = np.transpose(e_list)


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
       
       yield from fly2d(dets1, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)
       yield from bps.sleep(2)
       #h = db[-1]
       #sid = h.start['scan_id']
       plot2dfly(-1,'As')
       insertFig(note ='e = {}'.format(energy),title = ' ')
       plt.close()
       yield from bps.sleep(2)
		
		while (sclr2_ch4.get() < 5000):
            yield from bps.sleep(60)
            print('IC3 is lower than 5000, waiting...')
        
    save_page()

def zp_xanes2d(param_file, mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t):
    e_list = np.loadtxt(param_file)
    zp_list_xanes2d(e_list,mot1,x_s,x_e,x_num,mot2,y_s,y_e,y_num,accq_t)



