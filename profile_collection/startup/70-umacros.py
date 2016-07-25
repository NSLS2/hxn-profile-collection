import matplotlib.pyplot as plt
import numpy as np
from time import sleep


from scipy.optimize import curve_fit
from databroker import db, get_table
from ophyd import mov, movr

import sys
from datetime import datetime
import shutil

def focusmerlin(cnttime):
    merlin1.cam.acquire.put(0)
    merlin1.cam.acquire_time.put(cnttime)
    merlin1.cam.acquire_period.put(cnttime)
    merlin1.cam.trigger_mode.put(0)
    merlin1.cam.image_mode.put(2)
    sleep(.2)
    merlin1.cam.acquire.put(1)


def printfig():
    plt.savefig('/home/xf03id/Desktop/temp.png', bbox_inches='tight',
                pad_inches=4)
    os.system("lp -d HXN-printer-1 /home/xf03id/Desktop/temp.png")


def shutter(cmd):
    if cmd == 'open':
        shutter_open.put(1)
        sleep(5)
    elif cmd == 'close':
        shutter_close.put(1)
        sleep(5)


def mll_z_linecan(z_start, z_end, z_num, *,mot='ssy', start, end, num, acq_time, elem='Pt'):
    z_step = (z_end - z_start)/z_num
    c_ssz = smll.ssz.position
    movr_sz(z_start)
    for i in range(z_num + 1):
        if mot == 'ssy':
            RE(fly1d(ssy, start, end, num, acq_time))
        elif mot == 'ssx':
            RE(fly1d(ssy, start, end, num, acq_time))
        else:
            raise KeyError('mot has to be ssx or ssy')
        plot(-1, elem, 'sclr1_ch4')
        plt.title('ssz = %.3f' % smll.ssz.position)
    mov_sz(c_ssz)


def go_det(det):
    if det == 'merlin':
        diff.x.move(0, wait=False)
        sleep(0.5)
        diff.y1.move(0, wait=False)
        sleep(0.5)
        diff.y2.move(0, wait=False)
    elif det == 'cam11':
        diff.x.move(217.7, wait=False)
        sleep(0.5)
        diff.y1.move(20.25, wait=False)
        sleep(0.5)
        diff.y2.move(20.25, wait=False)
    elif det == 'tpx':
        mov(diff.x, -112)
        mov(diff.y1, -50)
        mov(diff.y2, -50)
    else:
        print('Inout det is not defined. '
              'Available ones are merlin, cam11 and tpx')


def go_energy(energy):
    energy = np.float(energy)
    if energy == 11950:
        mov(dcm_th, 9.52237)
        mov(dcm_p, 0.75671)
        mov(dcm_r, -0.12587)
        mov(zpx, -5297.64)
        mov(zpy, 4388.46)
        mov(zpz1, -35.2668)
        mov(ugap, 5.9182, wait=False)
        sleep(5)
    elif energy == 11874:
        mov(dcm_th, 9.58419)
        mov(dcm_p, 0.75485)
        mov(dcm_r, -0.12833)
        mov(zpx, -5292.64)
        mov(zpy, 4387.46)
        mov(zpz1, -34.1997)
        mov(ugap, 5.8803, wait=False)
        sleep(5)
    elif energy == 11870:
        mov(dcm_th, 9.58778)
        mov(dcm_p, 0.75547)
        mov(dcm_r, -0.12956)
        mov(zpx, -5290.64)
        mov(zpy, 4383.46)
        mov(zpz1, -34.1376)
        mov(ugap, 5.8803, wait=False)
        sleep(5)
    elif energy == 11869:
        mov(dcm_th, 9.58875)
        mov(dcm_p, 0.75547)
        mov(dcm_r, -0.12956)
        mov(zpx, -5290.64)
        mov(zpy, 4383.46)
        mov(zpz1, -34.1207)
        mov(ugap, 5.885, wait=False)
        sleep(5)
    elif energy == 11860:
        mov(dcm_th, 9.59594)
        mov(dcm_p, 0.75578)
        mov(dcm_r, -0.12915)
        mov(zpx, -5294.64)
        mov(zpy, 4384.46)
        mov(zpz1, -33.996)
        mov(ugap, 5.875, wait=False)
        sleep(5)
    else:
        print('energy not defined')

def sample_to_lab(xp, zp, alpha):
    x = np.cos(alpha)*xp + np.sin(alpha)*zp
    z = -np.sin(alpha)*xp + np.cos(alpha)*zp
    return(x, z)

def lab_to_sample(x, z, alpha):
    xp = np.cos(alpha)*x - np.sin(alpha)*z
    zp = np.sin(alpha)*x + np.cos(alpha)*z
    return(xp, zp)

def mll_mosaic_scan(x_start, x_end, x_num, x_block, y_start, y_end, y_num, y_block, acq_time, elem=None):

    max_travel = 500
    angle = 15.0*np.pi/180.0

    #initialize parameters
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_num = np.int(y_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    x_block = np.int(x_block)
    y_block = np.int(y_block)

    #read initial position
    pre_ssx = smll.ssx.position
    pre_ssy = smll.ssy.position
    pre_ssz = smll.ssz.position

    print('Initial ssx = ', pre_ssx)
    print('Initial ssy = ', pre_ssy)
    print('Initial ssz = ', pre_ssz)

    #calculate block size
    x_block_size = (x_end - x_start)
    y_block_size = (y_end - y_start)

    #move to first block
    dx = -(x_block*x_block_size/2.0 - 0.5*x_block_size)
    dy = -(y_block*y_block_size/2.0 - 0.5*y_block_size)

    if np.abs(dx) < max_travel and np.abs(dy) < max_travel:
        movr_sx(dx)
        movr_sy(dy)
    else:
        raise KeyError('Too large travel range')

    #start mosaic scan
    for i in range(y_block):
        for j in range(x_block):
            smll_sync_piezos()
            RE(fly2d(smll.ssx, x_start, x_end, x_num, smll.ssy, y_start, y_end, y_num, acq_time, return_speed=40))
            dx = x_block_size
            movr_sx(dx)
            if elem is not None:
                plot2dfly(-1, elem, 'sclr1_ch4')
        dx = -x_block_size*(x_block-1)
        dy = y_block_size
        movr_sx(dx)
        movr_sy(dy)

    #return to initial position
    print('Return to prior positions')
    mov_sx(pre_ssx)
    mov_sy(pre_ssy)

    print('%d x %d mosaic scan finished' % (x_block, y_block))


def mosaic_scan(x_start, x_end, x_num, y_start, y_end, y_num):
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_num = np.int(y_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)

    # kill close-loop
    zps.zp_kill_piezos.put(1)
    sleep(5)
    if x_num == 1:
        x_step = 0.0
    else:
        x_step = (x_end - x_start) / (x_num - 1)
    if y_num == 1:
        y_step = 0.0
    else:
        y_step = (y_end - y_start) / (y_num - 1)

    # read initial positions
    pre_x = zps.smarx.position
    pre_y = zps.smary.position
    pre_ssx = zps.zpssx.position
    pre_ssy = zps.zpssy.position

    print('x_step = ', x_step)
    print('y_step = ', y_step)
    print('Original smarx = ', pre_x)
    print('Original smary = ', pre_y)
    print('Original zpssx = ', pre_ssx)
    print('Original zpssy = ', pre_ssy)

    # move to start position
    movr(smarx, x_start)
    movr(smary, y_start)
    x = pre_ssx + (x_start * 1000)
    y = pre_ssy + (y_start * 1000)
    sleep(5)

    for i in range(y_num):
        for j in range(x_num):
            #            x = zps.zpssx.position
            #            print('zpssx =', x)
            #            sleep(5)
            mov(zpssx, x)
            sleep(2)
#            y = zps.zpssy.position
#            sleep(5)
#            print('zpssy =', y)
            mov(zpssy, y)
            sleep(2)
#            z = zps.zpssz.position
#            print('zpssz =', z)
#            sleep(5)
            mov(zpssz, 0)
            sleep(5)
            fly2d(zpssx, -10, 10, 100, zpssy, -10,
                  10, 100, 0.2, return_speed=40)
            print('scan finished, waiting for 120s...')
            sleep(120)
            zps.zp_kill_piezos.put(1)
            sleep(15)
            movr(smarx, x_step)
            x = x + (x_step * 1000)
        mov(smarx, pre_x)
        movr(smarx, x_start)
        x = pre_ssx + (x_start * 1000)
        movr(smary, y_step)
        y = y + (y_step * 1000)

    print('mosaic scan finished, move back to prior positions')
    mov(smarx, pre_x)
    mov(smary, pre_y)


def sin_offset(x, p0, p1, p2):
    return (p0 + p1 * np.sin((x + p2) * np.pi / 180)) / np.cos(x * np.pi / 180)


def sin_offset_fit(x, y, para):
    para = np.array(para)
    popt, pcov = curve_fit(sin_offset, x, y, para)
    # print(popt)
    y_fit = sin_offset(x, popt[0], popt[1], popt[2])
    return popt, pcov, y_fit


def rot_fit(x, y):
    x = np.array(x)
    y = -1 * np.array(y)

    para = [1, 1, 0]
    #para = [4.23077509,   0.58659241,   0.21648658, -19.8329533]
    popt, pcov, y_fit = sin_offset_fit(x, y, para)

    print(popt)
    plt.figure()
    plt.plot(x, y, label='data')
    plt.plot(x, y, 'go')
    plt.plot(x, y_fit, label='fit')
    plt.legend(loc='best')
    plt.title(np.str(popt[0]) + '+' + np.str(popt[1]) +
              '*sin(x+' + np.str(popt[2]) + ')')
    plt.xlabel('x:' + np.str(-1 * popt[1] * np.sin(popt[2] * np.pi / 180.)) +
               '  z:' + np.str(-1 * popt[1] * np.cos(popt[2] * np.pi / 180.)))
    plt.show()
    return(popt)

def coarse_align_rot(x, y, pix_size):
    r0, dr, offset = rot_fit(x,y)
    zps_kill_piezos()
    mov(zpsth, 0)
    dx = -dr*np.sin(offset)*pix_size/1000.0
    dz = -dr*np.cos(offset)*pix_size/1000.0

    movr(zps.smarx, dx)
    movr(zps.smarz, dz)



def linear_fit(x, y):
    x = np.array(x)
    y = np.array(y)

    p = np.polyfit(x, y, 1)
    y_fit = p[1] + p[0] * x
    plt.figure()
    plt.plot(x, y, 'go', label='points')
    plt.plot(x, y_fit, label='fit')
    plt.title(np.str(p[1]) + '+' + np.str(p[0]) + 'x')
    plt.show()


def inplane_angle(x, p0, p1):
    return p0 * np.sin(x * np.pi / 180) / np.sin((180 - x - p1) * np.pi / 180)


def inplane_angle_fit(x, y, para):
    para = np.array(para)
    popt, pcov = curve_fit(inplane_angle, x, y, para)
    # print(popt)
    y_fit = inplane_angle(x, popt[0], popt[1])
    return popt, pcov, y_fit


def inplane_fit(x, y):
    x = np.array(x)
    y = np.array(y)

    para = [572.7, 90]
    popt, pcov, y_fit = inplane_angle_fit(x, y, para)

    print(popt)
    plt.figure()
    plt.plot(x, y, label='data')
    plt.plot(x, y, 'go')
    plt.plot(x, y_fit, label='fit')
    plt.legend(loc='best')
    plt.title('r=' + np.str(popt[0]) + ' theta=' + np.str(popt[1]))
    plt.show()


def sin_offset_2(x, p0, p1, p2):
    return p0 + p1 * np.sin((x + p2) * np.pi / 180)


def sin_offset_fit_2(x, y, para):
    para = np.array(para)
    popt, pcov = curve_fit(sin_offset_2, x, y, para)
    # print(popt)
    y_fit = sin_offset_2(x, popt[0], popt[1], popt[2])
    return popt, pcov, y_fit


def rot_fit_2(x, y):
    x = np.array(x)
    y = np.array(y)

    para = [1, 1, -1]
    popt, pcov, y_fit = sin_offset_fit_2(x, y, para)

    print(popt)
    plt.figure()
    plt.plot(x, y, label='data')
    plt.plot(x, y, 'go')
    plt.plot(x, y_fit, label='fit')
    plt.legend(loc='best')
    plt.title(np.str(popt[0]) + '+' + np.str(popt[1]) +
              '*sin(x+' + np.str(popt[2]) + ')')
    plt.show()


def find_mass_center(array):
    n = np.size(array)
    tmp = 0
    for i in range(n):
        tmp += i * array[i]
    mc = np.round(tmp / np.sum(array))
    return mc


def tomo_scan(angle_start, angle_end, angle_num, x_start, x_end, x_num,
              y_start, y_end, y_num, exposure):

    angle_start = np.float(angle_start)
    angle_end = np.float(angle_end)
    angle_num = np.int(angle_num)
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    y_num = np.int(y_num)
    exposure = np.float(exposure)

    angle_step = (angle_end - angle_start) / angle_num

    for i in range(angle_num + 1):

        #while beamline_status.beam_current.get() <= 245:
        #    sleep(60)

        angle = angle_start + i * angle_step
        mov(zps.zpsth, angle)
        #x_start_real = x_start / np.cos(angle*np.pi/180.)
        #x_end_real = x_end / np.cos(angle*np.pi/180.)

        RE(dscan(zps.zpsx, -0.01, 0.01, 50, 0.5))
        scan_id, df = _load_scan(-1, fill_events=False)
        x = df['zpsx']
        data = df['Det1_Ce'] + df['Det2_Ce'] + df['Det3_Ce'] +\
            df['Det1_Gd'] + df['Det2_Gd'] + df['Det3_Gd'] +\
            df['Det1_Fe'] + df['Det2_Fe'] + df['Det3_Fe'] +\
            df['Det1_Co'] + df['Det2_Co'] + df['Det3_Co']
        x = np.asarray(x)
        data = np.asarray(data)
        #diff = np.diff(data)
        #i_max = np.where(diff == np.max(diff))
        #i_min = np.where(diff == np.min(diff))
        #i_center = np.round((i_max[0][0]+i_min[0][0])/2)+1

        mc = find_mass_center(data)
        mov(zps.zpsx, x[mc]-0.001)

        if np.abs(angle) <= 45:
            x_start_real = x_start / np.cos(angle * np.pi / 180.)
            x_end_real = x_end / np.cos(angle * np.pi / 180.)
            RE(fly2d(zpssx, x_start_real, x_end_real, x_num, zpssy,
                     y_start, y_end, y_num, exposure, return_speed=40))

        else:
            x_start_real = x_start / np.abs(np.sin(angle * np.pi / 180.))
            x_end_real = x_end / np.abs(np.sin(angle * np.pi / 180.))
            RE(fly2d(zpssz, x_start_real, x_end_real, x_num, zpssy,
                     y_start, y_end, y_num, exposure, return_speed=40))

        print('waiting for 60 sec...')
        sleep(10)

    mov(zps.zpsth, 0)


def tomo_slice_scan(angle_start, angle_end, angle_num, x_start, x_end, x_num,
                    y_start, y_end, y_num, exposure):
    angle_start = np.float(angle_start)
    angle_end = np.float(angle_end)
    angle_num = np.int(angle_num)
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    y_num = np.int(y_num)
    exposure = np.float(exposure)
    angle_step = (angle_end - angle_start) / angle_num

    y_step = (y_end - y_start) / y_num

    for j in range(y_num + 1):
        for i in range(angle_num + 1):
            angle = angle_start + i * angle_step
            mov(zpsth, angle)
            sleep(1)
            # mesh(zpssy,y_start,y_en,y_num,zpssx_lab,x_start,x_end,x_num,exposure)
            dscan(zpssx_lab, x_start, x_end, x_num, exposure)

            #print('waiting for 10 sec...')
            # sleep(10)
        movr(zpssy, y_step)

    mov(zpsth, 0)


def movr_zpz1(dz):
    movr(zp.zpz1, dz)
    movr(zp.zpx, dz * 3.75)
    movr(zp.zpy, -dz * 2.89)


def reset_tpx(num):
    for i in range(1000):
        timepix2.cam.num_images.put(num, wait=False)
        sleep(0.5)



def th_fly1d(th_start, th_end, num, m_start, m_end, m_num, sec):
    th_step = (th_end - th_start) / num
    movr(zpsth, th_start)
    for i in range(num + 1):
        fly1d(zpssy, m_start, m_end, m_num, sec)
        movr(zpsth, th_step)
    movr(zpsth, -(th_end + th_step))


def th_fly2d(th_start, th_end, num, x_start, x_end, x_num, y_start, y_end,
             y_num, sec):
    shutter('open')
    th_step = (th_end - th_start) / num
    movr(zps.zpsth, th_start)
    for i in range(num + 1):
        RE(fly2d(zpssz, x_start, x_end, x_num, zpssy, y_start, y_end, y_num, sec, return_speed=40))
        movr(zps.zpsth, th_step)
    movr(zps.zpsth, -(th_end + th_step))
    shutter('close')


def mov_diff(gamma, delta, r=500):
    diff_z = diff.z.position

    gamma = gamma * np.pi / 180
    delta = delta * np.pi / 180
    beta = 89.337 * np.pi / 180

    z_yaw = 574.668 + 581.20 + diff_z
    z1 = 574.668 + 395.2 + diff_z
    z2 = z1 + 380
    d = 395.2

    x_yaw = sin(gamma) * z_yaw / sin(beta + gamma)
    R_yaw = sin(beta) * z_yaw / sin(beta + gamma)
    R1 = R_yaw - (z_yaw - z1)
    R2 = R_yaw - (z_yaw - z2)
    y1 = tan(delta) * R1
    y2 = tan(delta) * R2
    R_det = R1 / cos(delta) - d
    dz = r - R_det

    print('Make sure all motors are zeroed properly, '
          'otherwise calculation will be wrong.')
    if x_yaw > 747 or x_yaw < -200:
        print('diff_x = ', -x_yaw,
              ' out of range, move diff_z upstream and try again')
    elif dz < -250 or dz > 0:
        print('diff_cz = ', dz,
              ' out of range, move diff_z up or down stream and try again')
    elif y1 > 450:
        print('diff_y1 = ', y1, ' out of range, move diff_z upstream '
              'and try again')
    elif y2 > 600:
        print('diff_y2 = ', y2, ' out of range, move diff_z upstream '
              'and try again')
    else:
        print('diff_x = ', -x_yaw, ' diff_cz = ', dz,
              ' diff_y1 = ', y1, ' diff_y2 = ', y2)
        print('wait for 1 sec, hit Ctrl+c to quit the operation')
        sleep(1)
        diff.y1.move(y1, wait=False)
        sleep(0.5)
        diff.y2.move(y2, wait=False)
        sleep(0.5)
        diff.x.move(-x_yaw, wait=False)
        sleep(0.5)
        diff.yaw.move(gamma * 180. / np.pi, wait=False)
        sleep(0.5)
        diff.cz.move(dz, wait=False)


def wh_diff():
    diff_z = diff.z.position
    diff_yaw = diff.yaw.position * np.pi / 180.0
    diff_cz = diff.cz.position
    diff_x = diff.x.position
    diff_y1 = diff.y1.position
    diff_y2 = diff.y2.position

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

    # print('x_yaw = ', x_yaw, ' diff_x = ', diff_x)
    if abs(x_yaw + diff_x) > 3:
        print('Not a pure gamma rotation')
    elif abs(diff_y1 / R1 - diff_y2 / R2) > 0.01:
        print('Not a pure delta rotation')
    else:
        delta = arctan(diff_y1 / R1)
        R_det = R1 / cos(delta) - d + diff_cz
        print('gamma = ', gamma * 180 / np.pi, ' delta = ',
              delta * 180 / np.pi, ' r = ', R_det)


def xanes_scan(bragg_list,x_start,x_end,x_num,y_start,y_end,y_num,exposure,sign='max'):
    bragg_list = np.array(bragg_list)
    num_bragg = np.size(bragg_list)
    current_det = gs.PLOT_Y
    gs.PLOT_Y='sclr1_ch4'
    for i in range(num_bragg):
        mov(dcm.th,bragg_list[i])

        RE(dscan(dcm.rf,-1, 1, 40, 1))
        if sign == 'max':
            mov(dcm.rf,gs.PS.max[0])
        elif sign == 'cen':
            mov(dcm.rf,gs.PS.cen)
        #df=get_table(db[-1],fill=False)
        #ic = np.asarray(df['sclr1_ch4'])
        #x = np.asarray(df['dcm_rf'])
        ##i_max = find_mass_center(ic)
        #mov(dcm.rf,x[ic == np.max(ic)])

        RE(dscan(m2.pf,-1, 1, 40, 1))
        if sign == 'max':
            mov(m2.pf,gs.PS.max[0])
        elif sign == 'cen':
            mov(m2.pf,gs.PS.cen)
        #df = get_table(db[-1],fill=False)
        #ic = np.asarray(df['sclr1_ch4'])
        #x = np.asarray(df[-1],'m2_pf')
        ##i_max = find_mass_center(ic)
        #mov(m2.pf,x[ic == np.max(ic)])

        RE(fly2d(ssx, x_start, x_end, x_num, ssy, y_start, y_end, y_num, exposure, return_speed=50))
    gs.PLOT_Y = current_det

def peak_ic():
    current_det = gs.PLOT_Y
    gs.PLOT_Y = 'sclr1_ch4'
    RE(dscan(dcm.rf,-2,2,80,1))
    mov(dcm.rf,gs.PS.max[0])
    RE(dscan(m2.pf,-2,2,80,1))
    mov(m2.pf,gs.PS.max[0])
    gs.PLOT_Y = current_det

def mono_m1(pf_start, pf_end, pf_num, b_start, b_end, b_num):
    b_step = (b_end - b_start)/b_num
    b_current = m1.b.position
    current_det = gs.PLOT_Y
    gs.PLOT_Y = 'sclr1_ch2'

    peak_int = np.zeros(b_num + 1);
    b_pos = np.zeros(b_num + 1)

    movr(m1.b, b_start)
    for i in range(b_num + 1):
        RE(dscan(dcm.pf, pf_start, pf_end, pf_num, 1))
        peak_int[i] = gs.PS.max[0]
        b_pos[i] = m1.b.position
        if i < b_num:
            movr(m1.b, b_step)
    mov(m1.b, b_current)
    gs.PLOT_Y = current_det
    plt.figure()
    plt.plot(b_pos, peak_int)
    plt.show()


def smll_kill_piezos():
    smll.kill.put(1)
    sleep(5)

def smll_zero_piezos():
    smll.zero.put(1)
    sleep(3)

def smll_sync_piezos():
    #sync positions
    mov(ssx, smll.ssx.position + 0.0001)
    mov(ssy, smll.ssy.position + 0.0001)
    mov(ssz, smll.ssz.position + 0.0001)

def movr_sx(dist):
    alpha = 15*np.pi/180.0
    c_ssx = smll.ssx.position
    c_ssy = smll.ssy.position
    c_ssz = smll.ssz.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    smll_kill_piezos()

    t_ssx = c_ssx + dist

    dxp = t_ssx - smll.ssx.position
    dzp = c_ssz - smll.ssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)

    movr(sx, dx)
    movr(sz, dz)

    dy = c_ssy -smll.ssy.position

    movr(sy, dy)

    sleep(5)

    smll_sync_piezos()

    mov(ssx, t_ssx)
    mov(ssy, c_ssy)
    mov(ssz, c_ssz)

    print('Post-move x = %.3f' % smll.ssx.position)
    print('Post-move y = %.3f' % smll.ssy.position)
    print('Post-move z = %.3f' % smll.ssz.position)

def mov_sx(t_pos):
    alpha = 15*np.pi/180.0
    c_ssx = smll.ssx.position
    c_ssy = smll.ssy.position
    c_ssz = smll.ssz.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    smll_kill_piezos()

    t_ssx = t_pos

    dxp = t_ssx - smll.ssx.position
    dzp = c_ssz - smll.ssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)

    movr(sx, dx)
    movr(sz, dz)

    dy = c_ssy - smll.ssy.position

    movr(sy, dy)

    sleep(5)

    smll_sync_piezos()
    mov(ssx, t_ssx)
    mov(ssy, c_ssy)
    mov(ssz, c_ssz)

    print('Post-move x = %.3f' % (smll.ssx.position))
    print('Post-move y = %.3f' % (smll.ssy.position))
    print('Post-move z = %.3f' % (smll.ssz.position))


def movr_sy(dist):
    alpha = 15*np.pi/180.0
    c_ssx = smll.ssx.position
    c_ssy = smll.ssy.position
    c_ssz = smll.ssz.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    smll_kill_piezos()

    t_ssy = c_ssy + dist
    dy = t_ssy - smll.ssy.position
    movr(sy, dy)

    dxp = c_ssx - smll.ssx.position
    dzp = c_ssz - smll.ssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)

    movr(sx, dx)
    movr(sz, dz)

    sleep(5)

    smll_sync_piezos()
    mov(ssx, c_ssx)
    mov(ssy, t_ssy)
    mov(ssz, c_ssz)

    print('Post-move x = %.3f' % smll.ssx.position)
    print('Post-move y = %.3f' % smll.ssy.position)
    print('Post-move z = %.3f' % smll.ssz.position)

def mov_sy(t_pos):
    alpha = 15*np.pi/180.0
    c_ssx = smll.ssx.position
    c_ssy = smll.ssy.position
    c_ssz = smll.ssz.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    smll_kill_piezos()

    t_ssy = t_pos
    dy = t_ssy - smll.ssy.position
    movr(sy, dy)

    dxp = c_ssx - smll.ssx.position
    dzp = c_ssz - smll.ssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)

    movr(sx, dx)
    movr(sz, dz)

    sleep(5)

    smll_sync_piezos()
    mov(ssx, c_ssx)
    mov(ssy, t_ssy)
    mov(ssz, c_ssz)

    print('Post-move x = %.3f' % (smll.ssx.position))
    print('Post-move y = %.3f' % (smll.ssy.position))
    print('Post-move z = %.3f' % (smll.ssz.position))

def movr_sz(dist):
    alpha = 15.0*np.pi/180

    c_ssz = smll.ssz.position
    c_ssy = smll.ssy.position
    c_ssx = smll.ssx.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    smll_kill_piezos()

    t_ssz = c_ssz + dist*np.cos(alpha)
    dz = t_ssz - smll.ssz.position
    dy = c_ssy - smll.ssy.position

    movr(sz, dz)
    movr(sy, dy)

    sleep(5)

    smll_sync_piezos()
    mov(ssy, c_ssy)
    mov(ssz, t_ssz)

    print('post-move x = %.3f' % smll.ssx.position)
    print('Post-move y = %.3f' % smll.ssy.position)
    print('Post-move z = %.3f' % smll.ssz.position)

def mov_sz(t_pos):
    alpha = 15.0*np.pi/180

    c_ssz = smll.ssz.position
    c_ssy = smll.ssy.position
    c_ssx = smll.ssx.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    smll_kill_piezos()

    t_ssz = t_pos
    dz = (t_ssz - smll.ssz.position)/np.cos(alpha)
    dy = c_ssy - smll.ssy.position

    movr(sz, dz)
    movr(sy, dy)

    sleep(5)

    smll_sync_piezos()
    mov(ssy, c_ssy)
    mov(ssz, t_ssz)

    print('Post-move x = %.3f' % smll.ssx.position)
    print('Post-move y = %.3f' % smll.ssy.position)
    print('Post-move z = %.3f' % smll.ssz.position)

def list_fly2d(x_pos, y_pos, scan_p):
    if np.size(x_pos) != np.size(y_pos):
        raise KeyError('size of x_pos list is not equal to that of y_pos list')
    else:
        num_pos = np.size(x_pos)
    if np.size(scan_p) != 7:
        raise KeyError('Last argument needs 7 numbers')
    else:
        for i in range(num_pos):
            mov_sx(x_pos[i])
            mov_sy(y_pos[i])
            RE(fly2d(ssx, scan_p[0], scan_p[1], scan_p[2], ssy, scan_p[3], scan_p[4], scan_p[5], scan_p[6], return_speed=20))
            plot2dfly(-1, 'Ca', 'sclr1_ch4')
            peak_ic()
            sleep(2)


def save_wh_pos(print_flag=False):
    class Tee(object):
        def __init__(self, *files):
            self.files = files
        def write(self, obj):
            for f in self.files:
                f.write(obj)
                f.flush() # If you want the output to be visible immediately
        def flush(self) :
            for f in self.files:
                f.flush()

    now = datetime.datetime.now()
    fn = '/data/wh_pos_printout/log-'+np.str(now.year)+'-'+np.str(now.month)+'-'+np.str(now.day)+'-'+np.str(now.hour)+'-'+np.str(now.minute)+'.log'
    f = open(fn,'w')
    original = sys.stdout
    sys.stdout = Tee(sys.stdout, f)
    wh_pos()
    sys.stdout = original
    f.close()
    if print_flag:
        shutil.copyfile(fn,'/data/wh_pos_printout/tmp.log')
        os.system("lp -o cpi=20 -o lpi=8 -o media='letter' -d HXN-printer-1 /data/wh_pos_printout/tmp.log")


def zps_kill_piezos():
    zps.kill.put(1)
    sleep(5)

def zps_zero_piezos():
    zps.zero.put(1)
    sleep(3)

def zps_sync_piezos():
    #sync positions
    mov(zps.zpssx, zps.zpssx.position + 0.0001)
    mov(zps.zpssy, zps.zpssy.position + 0.0001)
    mov(zps.zpssz, zps.zpssz.position + 0.0001)


def movr_zpsx(dist):
    alpha = 0.0*np.pi/180.0
    c_ssx = zps.zpssx.position
    c_ssy = zps.zpssy.position
    c_ssz = zps.zpssz.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    zps_kill_piezos()

    t_ssx = c_ssx + dist

    dxp = t_ssx - zps.zpssx.position
    dzp = c_ssz - zps.zpssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)
    dx = dx/1000.
    dz = dz/1000.

    movr(zps.zpsx, dx)
    movr(zps.zpsz, dz)

    dy = c_ssy -zps.zpssy.position
    dy = dy/1000.

    movr(smary, dy)

    sleep(5)

    zps_sync_piezos()

    mov(zps.zpssx, t_ssx)
    mov(zps.zpssy, c_ssy)
    mov(zps.zpssz, c_ssz)

    print('Post-move x = %.3f' % zps.zpssx.position)
    print('Post-move y = %.3f' % zps.zpssy.position)
    print('Post-move z = %.3f' % zps.zpssz.position)

def mov_zpsx(t_pos):
    alpha = 0.0*np.pi/180.0
    c_ssx = zps.zpssx.position
    c_ssy = zps.zpssy.position
    c_ssz = zps.zpssz.position

    print('Current zpssx = %.3f' % c_ssx)
    print('Current zpssy = %.3f' % c_ssy)
    print('Current zpssz = %.3f' % c_ssz)

    zps_kill_piezos()

    t_ssx = t_pos

    dxp = t_ssx - zps.zpssx.position
    dzp = c_ssz - zps.zpssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)
    dx = dx/1000.0
    dz = dz/1000.0

    movr(zps.zpsx, dx)
    movr(zps.zpsz, dz)

    dy = c_ssy - zps.zpssy.position
    dy = dy/1000.0

    movr(smary, dy)

    sleep(5)

    zps_sync_piezos()
    mov(zps.zpssx, t_ssx)
    mov(zps.zpssy, c_ssy)
    mov(zps.zpssz, c_ssz)

    print('Post-move x = %.3f' % (zps.zpssx.position))
    print('Post-move y = %.3f' % (zps.zpssy.position))
    print('Post-move z = %.3f' % (zps.zpssz.position))


def movr_zpsy(dist):
    alpha = 0.0*np.pi/180.0
    c_ssx = zps.zpssx.position
    c_ssy = zps.zpssy.position
    c_ssz = zps.zpssz.position

    print('Current zpssx = %.3f' % c_ssx)
    print('Current zpssy = %.3f' % c_ssy)
    print('Current zpssz = %.3f' % c_ssz)

    zps_kill_piezos()

    t_ssy = c_ssy + dist
    dy = t_ssy - zps.zpssy.position
    dy = dy/1000.0

    movr(smary, dy)

    dxp = c_ssx - zps.zpssx.position
    dzp = c_ssz - zps.zpssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)
    dx = dx/1000.0
    dz = dz/1000.0

    movr(zps.zpsx, dx)
    movr(zps.zpsz, dz)

    sleep(5)

    zps_sync_piezos()
    mov(zps.zpssx, c_ssx)
    mov(zps.zpssy, t_ssy)
    mov(zps.zpssz, c_ssz)

    print('Post-move x = %.3f' % zps.zpssx.position)
    print('Post-move y = %.3f' % zps.zpssy.position)
    print('Post-move z = %.3f' % zps.zpssz.position)

def mov_zpsy(t_pos):
    alpha = 0.0*np.pi/180.0
    c_ssx = zps.zpssx.position
    c_ssy = zps.zpssy.position
    c_ssz = zps.zpssz.position

    print('Current zpssx = %.3f' % c_ssx)
    print('Current zpssy = %.3f' % c_ssy)
    print('Current zpssz = %.3f' % c_ssz)

    zps_kill_piezos()

    t_ssy = t_pos
    dy = t_ssy - zps.zpssy.position
    dy = dy/1000.0

    movr(smary, dy)

    dxp = c_ssx - zps.zpssx.position
    dzp = c_ssz - zps.zpssz.position

    dx, dz = sample_to_lab(dxp, dzp, alpha)
    dx = dx/1000.0
    dz = dz/1000.0

    movr(zps.zpsx, dx)
    movr(zps.zpsz, dz)

    sleep(5)

    zps_sync_piezos()
    mov(zps.zpssx, c_ssx)
    mov(zps.zpssy, t_ssy)
    mov(zps.zpssz, c_ssz)

    print('Post-move x = %.3f' % (zps.zpssx.position))
    print('Post-move y = %.3f' % (zps.zpssy.position))
    print('Post-move z = %.3f' % (zps.zpssz.position))

def movr_zpsz(dist):
    alpha = 0.0*np.pi/180

    c_ssz = zps.zpssz.position
    c_ssy = zps.zpssy.position
    c_ssx = zps.zpssx.position

    print('Current zpssx = %.3f' % c_ssx)
    print('Current zpssy = %.3f' % c_ssy)
    print('Current zpssz = %.3f' % c_ssz)

    zps_kill_piezos()

    t_ssz = c_ssz + dist*np.cos(alpha)
    dz = t_ssz - zps.zpssz.position
    dy = c_ssy - zps.zpssy.position
    dx = c_ssx - zps.zpssx.position

    dz = dz/1000.0
    dy = dy/1000.0
    dx = dx/1000.0

    movr(zps.zpsz, dz)
    movr(smary, dy)
    movr(zps.zpsx, dx)

    sleep(5)

    zps_sync_piezos()
    mov(zpssy, c_ssy)
    mov(zpssz, t_ssz)
    mov(zpssx, c_ssx)

    print('post-move x = %.3f' % zps.zpssx.position)
    print('Post-move y = %.3f' % zps.zpssy.position)
    print('Post-move z = %.3f' % zps.zpssz.position)

def mov_zpsz(t_pos):
    alpha = 0.0*np.pi/180

    c_ssz = zps.zpssz.position
    c_ssy = zps.zpssy.position
    c_ssx = zps.zpssx.position

    print('Current zpssx = %.3f' % c_ssx)
    print('Current zpssy = %.3f' % c_ssy)
    print('Current zpssz = %.3f' % c_ssz)

    zps_kill_piezos()

    t_ssz = t_pos
    dz = (t_ssz - zps.zpssz.position)/np.cos(alpha)
    dy = c_ssy - zps.zpssy.position
    dx = c_ssx - zps.zpssx.position

    dx = dx/1000.0
    dy = dy/1000.0
    dz = dz/1000.0

    movr(zps.zpsz, dz)
    movr(smary, dy)
    movr(zps.zpsx, dx)

    sleep(5)

    zps_sync_piezos()

    mov(zpssy, c_ssy)
    mov(zpssz, t_ssz)
    mov(zpssx, c_ssx)

    print('Post-move x = %.3f' % zps.zpssx.position)
    print('Post-move y = %.3f' % zps.zpssy.position)
    print('Post-move z = %.3f' % zps.zpssz.position)
