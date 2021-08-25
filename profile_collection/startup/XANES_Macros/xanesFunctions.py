

def peak_the_flux():
    
    yield from bps.sleep(2)
    yield from peak_bpm_y(-2,2,4)
    yield from peak_bpm_x(-15,15,5)
    yield from peak_bpm_y(-2,2,4)
    

    
def move_energy(e_,ugap_,zpz_,crl_th_, ignoreCRL= False, ignoreZPZ = False):
    yield from bps.sleep(1)
            
    #tuning the scanning pv on to dispable c bpms
    caput('XF:03IDC-ES{Status}ScanRunning-I', 1)  

    yield from bps.mov(e,e_)
    yield from bps.sleep(1)
    yield from bps.mov(ugap, ugap_)
    yield from bps.sleep(2)
    if not ignoreZPZ: yield from mov_zpz1(zpz_)
    yield from bps.sleep(1)
    if not ignoreCRL: yield from bps.mov(crl.p,crl_th_)
    yield from bps.sleep(1)
            
    caput('XF:03IDC-ES{Status}ScanRunning-I', 0) #scan status off
