import numpy as np

def ptycho_z_scan(x_start,x_end,x_num,y_start,y_end,y_num,exposure):
    ssz_steps = np.arange(-10,11,2)
    num_steps = np.size(ssz_steps)
    for i in range(num_steps):
        print(i,ssz_steps[i])
        mov(smlld.dssz,ssz_steps[i])
        RE(fly2d(smlld.dssx,x_start,x_end,x_num,smlld.dssy,y_start,y_end,y_num,exposure))
        sleep(0.1)
    mov(smlld.dssz,0)
    print('repeat twice at ssz = 0')
    RE(fly2d(smlld.dssx,x_start,x_end,x_num,smlld.dssy,y_start,y_end,y_num,exposure))
    RE(fly2d(smlld.dssx,x_start,x_end,x_num,smlld.dssy,y_start,y_end,y_num,exposure))

def ptycho_partial_tomo(x_start,x_end,x_num,y_start,y_end,y_num,exposure):
    '''
    angle_list = np.arange(-30,31,5)
    angle_list = angle_list[::-1]
    tomo_scan_list(angle_list,x_start,x_end,x_num,y_start,y_end,y_num,exposure)
    '''
    mov(smlld.dsx,-413.492)
    mov(smlld.dsy,1.874)
    mov(smlld.dsz,-955.4346)
    mov(smlld.dsth,90)
    sleep(0.1)
    RE(fly2d(smlld.dssz,x_start,x_end,x_num,smlld.dssy,y_start,y_end,y_num,exposure))
    RE(fly2d(smlld.dssz,x_start,x_end,x_num,smlld.dssy,y_start,y_end,y_num,exposure))
    movr(smlld.sbz,50)
    RE(fly2d(smlld.dssz,-2,2,80,smll.dssy,-1,6,140,0.1))
