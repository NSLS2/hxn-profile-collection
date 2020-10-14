
from epics import caput,caget


def save_scan_info(sid):#,export_folder):
    # save baseline, detector angle and roi setting for a scan
    bl = db[sid].table('baseline')
    sid, df = _load_scan(sid, fill_events=False)
    hd = db[sid].start
    exp_s = hd['exposure_time']
    #path = os.path.join(export_folder, 'scan_{}_info.txt'.format(sid))

    #bl.to_csv(path, float_format='%1.5e', sep='\t')


    #header = db[sid].start
    #header.to_csv(path, float_format='%1.5e', sep='\t')

    roi_x0 = caget('XF:03IDC-ES{Merlin:1}ROI1:MinX')
    roi_y0 = caget('XF:03IDC-ES{Merlin:1}ROI1:MinY')
    roi_nx = caget('XF:03IDC-ES{Merlin:1}ROI1:SizeX')
    roi_ny = caget('XF:03IDC-ES{Merlin:1}ROI1:SizeY')

    diff_z = np.array(bl['diff_z'])[0]
    diff_yaw = np.array(bl['diff_yaw'])[0] * np.pi / 180.0
    diff_cz = np.array(bl['diff_cz'])[0]
    diff_x = np.array(bl['diff_x'])[0]
    diff_y1 = np.array(bl['diff_y1'])[0]
    diff_y2 = np.array(bl['diff_y2'])[0]

    gamma = diff_yaw
    beta = 89.337 * np.pi / 180
    z_yaw = 574.668 + 581.20 + diff_z
    z1 = 574.668 + 395.2 + diff_z
    z2 = z1 + 380
    d = 395.2

    x_yaw = sin(gamma) * z_yaw / sin(beta + gamma)
    R_yaw = sin(beta) * z_yaw / sin(beta + gamma)
    R1 = R_yaw - (z_yaw - z1)
    R2 = R_yaw - (z_yaw - z2)

    if abs(x_yaw + diff_x) > 3:
        gamma = 0
        delta = 0
        #R_det = 500

        beta = 89.337 * np.pi / 180
        R_yaw = np.sin(beta) * z_yaw / np.sin(beta + gamma)
        R1 = R_yaw - (z_yaw - z1)
        R_det = R1 / np.cos(delta) - d + diff_cz

    elif abs(diff_y1 / R1 - diff_y2 / R2) > 0.01:
        gamma = 0
        delta = 0
        #R_det = 500

        beta = 89.337 * np.pi / 180
        R_yaw = np.sin(beta) * z_yaw / np.sin(beta + gamma)
        R1 = R_yaw - (z_yaw - z1)
        R_det = R1 / np.cos(delta) - d + diff_cz

    else:
        delta = arctan(diff_y1 / R1)
        R_det = R1 / cos(delta) - d + diff_cz

    print('gamma, delta, dist:', gamma*180/np.pi, delta*180/np.pi, R_det)
    print('ROI: ', roi_x0, roi_y0, roi_nx, roi_ny )
    '''
    with open(path, 'a') as file:
        file.write('\n')
        file.write('det info: \n')
        file.write('gamma = %f \n' % (gamma * 180 / np.pi))
        file.write('delta = %f \n' % (delta * 180 / np.pi))
        file.write('r = %f \n' % R_det)
        file.write('\n')
        file.write('ROI x start = %f \n' % roi_x0)
        file.write('ROI y start = %f \n' % roi_y0)
        file.write('ROI x size = %f \n' % roi_nx)
        file.write('ROI y size = %f \n' % roi_ny)
        file.write('\n')
        file.write('exposure time = %f \n' % exp_s)
        file.write('\n')
    '''

def fly1d_user(motor,start,end,num_pos,exp):
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"R1")
    RE(fly1d(motor,start,end,num_pos,exp))
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"I1")


def fly2d_user(motor1,start1,end1,num1,motor2,start2,end2,num2,exp):
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"R1")
    RE(fly2d(motor1,start1,end1,num1,motor2,start2,end2,num2,exp))
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"I1")


def theta_fly2d_mll(angle_start,angle_end,num_angle,start1,end1,num1,start2,end2,num2,exp):
    angle_current = smlld.dsth.position
    dssx_current = smlld.dssx.position
    dssy_current = smlld.dssy.position
    dssz_current = smlld.dssz.position

    angle_step = (angle_end - angle_start) / num_angle
    start_angle = angle_current + angle_start
    yield from bps.mov(smlld.dsth,start_angle)
    #RE(fly1d(zps.zpssz,-1,1,50,0.1))
    #mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)



    for i in range(num_angle+1):
        while (sclr2_ch4.get() < 400000):
            yield from bps.sleep(60)
            print('IC3 is lower than 400000, waiting...')

        yield from fly1d(dets1,smlld.dssy,-2,2,100,0.05)
        yc = return_line_center(-1,'Ni',0.2)
        yield from bps.mov(smlld.dssy,yc)


        #yield from bps.mov(smlld.dssy,0.7)
        yield from fly1d(dets1,smlld.dssx,-5,5,100,0.05)
        #a,b = erf_fit(-1,'Au_L')
        a = return_line_center(-1,'Ni',0.2)
        #plt.close()
        yield from bps.mov(smlld.dssx,a)

        #yield from bps.mov(smlld.dssy,0)
        yield from fly2d(dets1,smlld.dssx,start1,end1,num1,smlld.dssy,start2,end2,num2,exp,dead_time=0.002)
        #yield from bps.movr(dssx,0.06)
        #mov_to_image_cen_zpss(-1,'Ga',1)
        #export(-1)
        yield from bps.movr(smlld.dsth,angle_step)
        #RE(fly1d(zps.zpssz,-.25,.25,20,0.1))
        #mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
        #xspress3.unstage()
        print('wait 2 sec ...')
        yield from bps.sleep(2)

    yield from bps.mov(smlld.dsth,angle_current)

def theta_fly2d_mll_blank(angle_start,angle_end,num_angle,start1,end1,num1,start2,end2,num2,exp):
    angle_current = smlld.dsth.position
    dssx_current = smlld.dssx.position
    dssy_current = smlld.dssy.position
    dssz_current = smlld.dssz.position

    angle_step = (angle_end - angle_start) / num_angle
    start_angle = angle_current + angle_start
    yield from bps.mov(smlld.dsth,start_angle)
    #RE(fly1d(zps.zpssz,-1,1,50,0.1))
    #mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)

    for i in range(num_angle+1):
        while (sclr2_ch4.get() < 70000):
            yield from bps.sleep(60)
            print('IC3 is lower than 70000, waiting...')

        yield from fly1d(dets1,smlld.dssx,9,11,200,0.1)
        a,b = erf_fit(-1,'Au_L')
        plt.close()
        yield from bps.mov(smlld.dssz,a-10)

        yield from fly1d(dets1,smlld.dssy,9,11,200,0.1)
        a,b = erf_fit(-1,'Au_L')
        plt.close()
        yield from bps.mov(smlld.dssy,a-10)


        yield from fly2d(dets2,smlld.dssz,start1,end1,num1,smlld.dssy,start2,end2,num2,exp,return_speed=40)
        #mov_to_image_cen_zpss(-1,'Ga',1)
        #export(-1)
        yield from bps.movr(smlld.dsth,angle_step)
        #RE(fly1d(zps.zpssz,-.25,.25,20,0.1))
        #mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
        xspress3.unstage()
        print('wait 2 sec ...')
        yield from bps.sleep(2)



def searchtheta(angle_start,angle_end,num_angle,motor1,start1,end1,num1):
    angle_current = zps.zpsth.position
    smarx_current = zps.smarx.position
    smary_current = zps.smary.position
    smarz_current = zps.smarz.position
    ssx_current = zps.zpssx.position
    ssy_current = zps.zpssy.position
    ssz_current = zps.zpssz.position

    angle_step = (angle_end - angle_start) / num_angle
    start_angle = angle_current + angle_start
    mov(zps.zpsth,start_angle)
    RE(fly1d(zps.zpssz,-2,2,50,0.1))
    mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)

    for i in range(num_angle+1):
        movr(zps.zpsth,angle_step)
        RE(fly1d(zps.zpssz,-2,2,40,0.1))
        mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
        print('wait 2 sec ...')
        sleep(2)


def theta_fly2d_ce(angle_start,angle_end,num_angle,motor1,start1,end1,num1,motor2,start2,end2,num2,exp):
    angle_current = zps.zpsth.position
    smarx_current = zps.smarx.position
    smary_current = zps.smary.position
    smarz_current = zps.smarz.position
    ssx_current = zps.zpssx.position
    ssy_current = zps.zpssy.position
    ssz_current = zps.zpssz.position

    angle_step = np.int((angle_end - angle_start) / num_angle)
    start_angle = angle_current + angle_start
    yield from bps.mov(zps.zpsth,start_angle)
   # RE(fly1d(zps.zpssz,-1,1,50,0.1))
   # mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)

    for i in range(num_angle+1):
        yield from fly2d(dets1,motor1,start1,end1,num1,motor2,start2,end2,num2,exp)
    #    mov_to_image_cen_zpss(-1,'Ga',1)
       # export(-1)
        yield from bps.movr(zps.zpsth,angle_step)
        #yield from bps.movr(zps.zpssz,-2*angle_step)
    #    RE(fly1d(zps.zpssz,-.25,.25,20,0.1))
    #    mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
        print('wait 2 sec ...')
        sleep(2)

    yield from bps.mov(zps.zpsth,angle_current)

def search_smarx():
    movr(zps.smarx, -.1)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)
    print('10sec pause')
    sleep(10)
    RE(fly2d(zpssz,-15,15,60, zpssy, -15,15,60, .1))
    movr(zps.smarx, -.1)

def vlm_tomo(angle_start,angle_end,angle_step):
    angle_0 = dsth.position
    angle_list = np.arange(angle_start,angle_end+0.1,angle_step)
    num_angle = np.size(angle_list)
    mov(dsth,angle_list[0])
    for i in range(num_angle):
        mov(dsth,angle_list[i])
        sleep(1)
        caput('XF:03IDC-ES{CAM:10}TIFF1:WriteFile',1)
    mov(dsth,angle_0)


def vlm_orth(start,end,step,angle):
    x_list = np.arange(start,end+1,step)
    y_list = np.arange(start,end+1,step)
    num = np.size(x_list)
    mov(dsth,angle)
    for i in range(num):
        mov(dssy,y_list[i])
        for j in range(num):
            if angle == 0:
                mov(dssx,x_list[j])
            else:
                mov(dssz,x_list[j])
            sleep(2)
            caput('XF:03IDC-ES{CAM:10}TIFF1:WriteFile',1)
        sleep(2)

    mov(dssy,0)
    mov(dssx,0)
    mov(dssz,0)

def vlm_orth_coarse(start,end,step,angle):
    x_list = np.arange(start,end+1,step)
    y_list = np.arange(start,end+1,step)
    num = np.size(x_list)
    dsx_0 = dsx.position
    dsy_0 = dsy.position
    dsz_0 = dsz.position
    mov(dsth,angle)
    movr(dsy,start/1000.)

    for i in range(num):

        if angle == 0:
            mov(dsx,dsx_0+start)
        else:
            mov(dsz,dsz_0+start)

        for j in range(num):
            if angle == 0:
                movr(dsx,step)
            else:
                movr(dsz,step)
            sleep(1)
            caput('XF:03IDC-ES{CAM:10}TIFF1:WriteFile',1)
        movr(dsy,step/1000.)

    mov(dsx,dsx_0)
    mov(dsy,dsy_0)
    mov(dsz,dsz_0)

def peak_up_xbpm(start_relative,end_relative,step_size,direction='y'):
    print('peaking up xbmp set point in ', direction, 'direction')
    num_points = np.int((end_relative - start_relative) // step_size + 1)
    ic = np.zeros(num_points)

    if direction == 'x':
        setpoint_pv = 'XF:03ID{XBPM:17}Fdbk:A-SP'
    elif direction == 'y':
        setpoint_pv = 'XF:03ID{XBPM:17}Fdbk:B-SP'

    init_position = caget(setpoint_pv)

    position = init_position + start_relative + np.arange(num_points) * step_size

    for i in range(num_points):
        caput(setpoint_pv,position[i])
        sleep(2)
        ic[i] = caget('XF:03IDC-ES{Sclr:2}_cts1.D')
        print(direction, 'set point: ', position[i], ', IC 3 counts: ', ic[i])

    peak_position = position[ic == np.max(ic)]
    print('peak count at ',peak_position)
    caput(setpoint_pv,peak_position)


def nightscan():
    RE(fly2d(dssx,-0.5, 0.5, 100, dssy, -0.5, 0.5, 100, 0.1))
    sleep(5)
    RE(fly2d(dssx,-0.5, 0.5, 100, dssy, -0.5, 0.5, 100, 0.05))
    sleep(5)
    movr(dssy,1)
    RE(fly2d(dssx,-0.5, 0.5, 100, dssy, -0.5, 0.5, 100, 0.1))
    sleep(5)
    RE(fly2d(dssx,-0.5, 0.5, 100, dssy, -0.5, 0.5, 100, 0.1))
    movr(dssx, 1.)
    sleep(5)
    RE(fly2d(dssx,-0.5, 0.5, 100, dssy, -0.5, 0.5, 100, 0.1))
    sleep(5)
    RE(fly2d(dssx,-0.5, 0.5, 100, dssy, -0.5, 0.5, 100, 0.05))
    RE(fly2d(dssx,-0.5, 0.5, 100, dssy, -0.5, 0.5, 100, 0.5))
    movr(dssy, 2)
    movr(dssx,-2)
    RE(fly2d(dssx,-1,1,200,dssy,-1,1,200,0.2))
    shutter('close')


def recover_and_scan():
    yield from recover_mll_scan_pos(54029, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 120, dssy, -3, 3, 120, 0.05, return_speed=40)

    yield from recover_mll_scan_pos(54032, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 120, dssy, -3, 3, 120, 0.05, return_speed=40)

    yield from recover_mll_scan_pos(54035, True, False, False)
    yield from fly2d(dets1, dssx, -4, 4, 160, dssy, -4, 4, 160, 0.05, return_speed=40)

    yield from recover_mll_scan_pos(54038, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 120, dssy, -3, 3, 120, 0.05, return_speed=40)

    yield from recover_mll_scan_pos(54042, True, False, False)
    yield from fly2d(dets1, dssx, -5, 5, 200, dssy, -3, 3, 120, 0.05, return_speed=40)

    yield from recover_mll_scan_pos(54045, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 120, dssy, -3, 3, 120, 0.05, return_speed=40)



    yield from recover_mll_scan_pos(54029, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 240, dssy, -3, 3, 240, 0.05, return_speed=40)
    yield from fly2d(dets1, dssx, -3, 3, 120, dssy, -3, 3, 120, 0.1, return_speed=40)

    yield from recover_mll_scan_pos(54032, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 240, dssy, -3, 3, 240, 0.05, return_speed=40)
    yield from fly2d(dets1, dssx, -3, 3, 120, dssy, -3, 3, 120, 0.1, return_speed=40)

    #yield from recover_mll_scan_pos(54035, True, False, False)
    #yield from fly2d(dets1, dssx, -4, 4, 160, dssy, -4, 4, 160, 0.05, return_speed=40)

    yield from recover_mll_scan_pos(54038, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 240, dssy, -3, 3, 240, 0.05, return_speed=40)

    #yield from recover_mll_scan_pos(54042, True, False, False)
    #yield from fly2d(dets1, dssx, -5, 5, 200, dssy, -3, 3, 120, 0.05, return_speed=40)

    yield from recover_mll_scan_pos(54045, True, False, False)
    yield from fly2d(dets1, dssx, -3, 3, 240, dssy, -3, 3, 240, 0.05, return_speed=40)

import epics
def zp_mesh_scan(xs,xe,xn,ys,ye,yn,exposure):
    zpssx_0 = zpssx.position
    zpssy_0 = zpssy.position
    x_step_size = (xe-xs) / xn
    y_step_size = (ye-ys) / yn
    print(xs,xe,xn,x_step_size,ys,ye,yn,y_step_size)

    caput('XF:03IDC-ES{Merlin:2}cam1:Acquire',0)
    caput('XF:03IDC-ES{Merlin:1}cam1:Acquire',0)

    caput('XF:03IDC-ES{Merlin:1}cam1:AcquireTime',exposure)
    yield from bps.sleep(5)
    caput('XF:03IDC-ES{Merlin:2}cam1:AcquireTime',exposure)
    yield from bps.sleep(5)

    caput('XF:03IDC-ES{Merlin:1}cam1:AcquirePeriod',exposure+0.1)
    yield from bps.sleep(5)
    caput('XF:03IDC-ES{Merlin:2}cam1:AcquirePeriod',exposure+0.1)
    yield from bps.sleep(5)
    yield from bps.movr(zpssx,xs)
    yield from bps.movr(zpssy,ys)

    for i in range(xn+1):
        for j in range(yn+1):
            print(i,j,zpssx.position,zpssy.position)

            caput('XF:03IDC-ES{Merlin:2}TIFF1:Capture',1)
            yield from bps.sleep(2)
            caput('XF:03IDC-ES{Merlin:1}TIFF1:Capture',1)
            yield from bps.sleep(2)

            print('start exposure ...')
            caput('XF:03IDC-ES{Merlin:2}cam1:Acquire',1)
            yield from bps.sleep(2)
            caput('XF:03IDC-ES{Merlin:1}cam1:Acquire',1)
            yield from bps.sleep(2)

            print('wait ...')
            yield from bps.sleep(exposure+5)


            yield from bps.movr(zpssx,x_step_size)

        yield from bps.mov(zpssx,zpssx_0+xs)
        yield from bps.movr(zpssy,y_step_size)

    yield from bps.mov(zpssx,zpssx_0)
    yield from bps.mov(zpssy,zpssy_0)

def makeup_scan(sid_list,x_start, x_end, x_num, y_start, y_end, y_num, exposure, elem):
    num_scan = np.size(sid_list)
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    y_num = np.int(y_num)
    exposure = np.float(exposure)

    for i in range(num_scan):
        sid = np.int(sid_list[i])
        yield from recover_mll_scan_pos(sid,1,0,0)
        angle = dsth.position

        if np.abs(angle) <= 45:
            yield from fly1d(dets1,dssx, -10, 10, 200, 0.03)
            xc = return_line_center(-1,elem,0.1)
            yield from bps.mov(dssx,xc)
        else:
            yield from fly1d(dets1,dssz, -10, 10, 200, 0.03)
            xc = return_line_center(-1,elem,0.1)
            yield from bps.mov(dssz,xc)

        yield from fly1d(dets1,dssy, -8, 8, 160, 0.03)
        yc = return_line_center(-1,elem,0.1)
        yield from bps.mov(dssy,yc)

        merlin1.unstage()
        xspress3.unstage()

        while (sclr2_ch4.get() < 400000):
            yield from bps.sleep(60)
            print('IC3 is lower than 400000, waiting...')

        if np.abs(angle) <= 45:
            x_start_real = x_start / np.cos(angle * np.pi / 180.)
            x_end_real = x_end / np.cos(angle * np.pi / 180.)
            yield from fly2d(dets1, smlld.dssx,x_start_real,x_end_real,x_num,smlld.dssy,
                     y_start, y_end, y_num, exposure, return_speed=40)
        else:
            x_start_real = x_start / np.abs(np.sin(angle * np.pi / 180.))
            x_end_real = x_end / np.abs(np.sin(angle * np.pi / 180.))
            yield from fly2d(dets1, smlld.dssz,x_start_real,x_end_real,x_num, smlld.dssy,
                     y_start, y_end, y_num, exposure, return_speed = 40)

        merlin1.unstage()
        xspress3.unstage()
        print('waiting for 2 sec...')
        yield from bps.sleep(2)

def night_mosaic_xh(nx,ny):
    ic_0 = sclr2_ch4.get()

    yield from fly2d(dets1,zpssx,-15,15,nx,zpssy,-15,15,ny,0.05,return_speed=40)
    while (sclr2_ch4.get() < (0.9*ic_0)):
        yield from peak_bpm_y(-5,5,10)
    yield from bps.sleep(1)
    yield from bps.movr(smarx,0.02)
    yield from fly2d(dets1,zpssx,-15,15,nx,zpssy,-15,15,ny,0.05,return_speed=40)
    while (sclr2_ch4.get() < (0.9*ic_0)):
        yield from peak_bpm_y(-5,5,10)
    yield from bps.sleep(1)
    yield from bps.movr(smary,0.02)
    yield from fly2d(dets1,zpssx,-15,15,nx,zpssy,-15,15,ny,0.05,return_speed=40)
    while (sclr2_ch4.get() < (0.9*ic_0)):
        yield from peak_bpm_y(-5,5,10)
    yield from bps.sleep(1)
    yield from bps.movr(smarx,-0.02)
    yield from fly2d(dets1,zpssx,-15,15,nx,zpssy,-15,15,ny,0.05,return_speed=40)
    yield from bps.movr(smary,-0.02)

def night_mosaic_mp(nx,ny):
     yield from bps.mov(smarx,-1.648)
     yield from bps.mov(smary,1.5175)
     yield from night_mosaic_xh(200,200)

     yield from bps.mov(smarx,-1.803)
     yield from bps.mov(smary,1.5875)
     yield from night_mosaic_xh(200,200)

     #yield from bps.mov(smarx,-1.838)
     #yield from bps.mov(smary,0.646)
     #yield from night_mosaic_xh(200,200)



import epics
def theta_dexela(angle_start,angle_end,angle_step,exposure_time):
    theta_zero = dsth.position
    angle_step_num = np.int((angle_end - angle_start) / angle_step)
    print('number of steps:', angle_step_num)

    caput('XF:03IDC-ES{Dexela:1}cam1:AcquireTime',exposure_time)
    caput('XF:03IDC-ES{Dexela:1}cam1:AcquirePeriod',exposure_time+0.2)

    yield from bps.mov(dsth,angle_start)
    for i in range(angle_step_num+1):
        print('dsth angle:', dsth.position)
        caput('XF:03IDC-ES{Dexela:1}TIFF1:Capture',1)
        yield from bps.movr(dsth,angle_step)
        yield from bps.sleep(exposure_time+0.2)

    yield from bps.mov(dsth,theta_zero)


def repeat_2d(zs,ze,z_num):
    z_0 = zps.zpsz.position
    z_step = (ze - zs) / z_num
    yield from bps.mov(zps.zpsz,(z_0+zs))
    for i in range(z_num+1):
        yield from fly2d(dets1,zpssx,2,8,120,zpssy,-7,-4,60,0.05,return_speed=40)
        yield from bps.sleep(2)
        plot2dfly(-1,'Ni')
        yield from bps.movr(zps.zpsz,z_step)

    yield from bps.mov(zps.zpsz,z_0)

