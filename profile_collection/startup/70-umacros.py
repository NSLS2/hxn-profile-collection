import matplotlib.pyplot as plt
import numpy as np
from time import sleep
import os

from scipy.optimize import curve_fit
from ophyd import mov, movr
from scipy import ndimage

import sys
from datetime import datetime
import shutil

import scipy
from scipy import signal

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
        shutter_open.put(1)
        sleep(5)
    elif cmd == 'close':
        shutter_close.put(1)
        sleep(5)
        shutter_close.put(1)
        sleep(5)


def mll_z_linescan(z_start, z_end, z_num, mot, start, end, num, acq_time, elem='Pt_L'):
    z_step = (z_end - z_start)/z_num
    init_sz = smlld.sbz.position
    movr(smlld.sbz, z_start)
    for i in range(z_num + 1):
        if mot == 'dssy':
            RE(fly1d(dssy, start, end, num, acq_time))
        elif mot == 'dssx':
            RE(fly1d(dssx, start, end, num, acq_time))
        else:
            raise KeyError('mot has to be dssx or dssy')
        plot(-1, elem, 'sclr1_ch4')
        plt.title('sbz = %.3f' % smlld.sbz.position)
        movr(smlld.sbz, z_step)
    mov(smlld.sbz, init_sz)

def mll_z_2dscan(z_start, z_end, z_num, mot1, start1, end1, num1, mot2, start2, end2, num2, acq_time, elem='Au'):
    z_step = (z_end - z_start)/z_num
    init_sz = smlld.sbz.position
    movr(smlld.sbz, z_start)
    for i in range(z_num + 1):
        RE(fly2d(mot1, start1, end1, num1, mot2, start2, end2, num2, acq_time))
        plot2dfly(-1, elem, 'sclr1_ch4')
        plt.title('sbz = %.3f' % smlld.sbz.position)
        movr(smlld.sbz, z_step)
    mov(smlld.sbz, init_sz)

def zp_z_2dscan(z_start, z_end, z_num, mot1, start1, end1, num1, mot2, start2, end2, num2, acq_time, elem='Co'):
    z_step = (z_end - z_start)/z_num
    init_sz = zps.zpsz.position
    movr(zps.zpsz, z_start/1000)
    for i in range(z_num + 1):
        RE(fly2d(mot1, start1, end1, num1, mot2, start2, end2, num2, acq_time))
        plot2dfly(-1, elem, 'sclr1_ch4')
        plt.title('zpsz = %.3f' % zps.zpsz.position)
        movr(zps.zpsz, z_step/1000)
    mov(zps.zpsz, init_sz)

def go_det(det):
    if det == 'merlin':
        diff.x.move(-2.9, wait=False)
        sleep(0.5)
        diff.y1.move(-3.2, wait=False)
        sleep(0.5)
        diff.y2.move(-3.2, wait=False)
    elif det == 'cam11':
        diff.x.move(211.62, wait=False)
#        diff.x.move(211.32, wait=False)
        sleep(0.5)
        diff.y1.move(22.5, wait=False)
        sleep(0.5)
        diff.y2.move(22.5, wait=False)
    elif det == 'tpx':
        mov(diff.x, -112)
        mov(diff.y1, -50)
        mov(diff.y2, -50)
    elif det =='telescope':
        mov(diff.x, -342.)
        mov(diff.z,-50)
    else:
        print('Inout det is not defined. '
              'Available ones are merlin, cam11, telescope and tpx')


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
    #zps.zp_kill_piezos.put(1)
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
    x_ini=pre_x + x_start
    y_ini = pre_y + y_start
    movr(smarx, x_start)
    movr(smary, y_start)
    x = pre_ssx + (x_start * 1000)
    y = pre_ssy + (y_start * 1000)
    sleep(5)


    for i in range(y_num):

        for j in range(x_num):
            print(i,j,zps.smarx.position,zps.smary.position)
            RE(fly2d(zpssx, -15, 15, 30, zpssy, -15, 15, 30, 0.05, return_speed=40))
            scan_id,df=_load_scan(-1,fill_events=False)
            current_smarx = zps.smarx.position
            current_smary = zps.smary.position
            '''
            plot2dfly(-1,'Zn',norm='sclr1_ch4')
            plt.title('#'+np.str(scan_id)+', smarx '+np.str(current_smarx)+', smary '+np.str(current_smary))
            printfig()
            plot2dfly(-1, 'sclr1_ch5')
            plt.title('#'+np.str(scan_id))
            printfig()
            print('scan finished, waiting for 2s...')
            '''
           # zps.zp_kill_piezos.put(1)
            sleep(2)
            movr(smarx, x_step)

        mov(smarx, x_ini)
        movr(smary, y_step)

    print('mosaic scan finished, move back to prior positions')
    mov(smarx, pre_x)
    mov(smary, pre_y)
    shutter('close')


def sin_offset(x, p0, p1, p2):
    return (p0 + p1 * np.sin((x + p2) * np.pi / 180.)) / np.cos(x * np.pi / 180.)


def sin_offset_fit(x, y, para):
    para = np.array(para)
    popt, pcov = curve_fit(sin_offset, x, y, para)
    print(popt)
    y_fit = sin_offset(x, popt[0], popt[1], popt[2])
    print(x,y_fit)
    return popt, pcov, y_fit


def rot_fit(x, y):
    x = np.array(x)
    y = -1 * np.array(y)

    para = [0.1, -0.1, 0]
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
    r0, dr, offset = rot_fit_2(x,y)
    zps_kill_piezos()
    mov(zps.zpsth, 0)
    dx = -dr*np.sin(offset*np.pi/180)*pix_size/1000.0
    dz = -dr*np.cos(offset*np.pi/180)*pix_size/1000.0
    print(dx,dz)
    #movr(zps.smarx, dx)
    #movr(zps.smarz, dz)
    movr(smlld.dsx,dx*1000.)
    movr(smlld.dsz,dz*1000.)


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
    return popt[0], popt[1], popt[2]

def find_mass_center(array):
    n = np.size(array)
    tmp = 0
    for i in range(n):
        tmp += i * array[i]
    mc = np.round(tmp / np.sum(array))
    return mc

def find_mass_center_1d(array,x):
    n = np.size(array)
    tmp = 0
    for i in range(n):
        tmp += x[i] *array[i]
    mc = tmp / np.sum(array)
    return mc

def mov_to_image_center_tmp(scan_id=-1, elem='W_L', bitflag=1, moveflag=1,piezomoveflag=1):
    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    hdr = db[scan_id]['start']
    x_motor = hdr['motor1']
    y_motor = hdr['motor2']
    x = np.asarray(df2[x_motor])
    y = np.asarray(df2[y_motor])
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    if bitflag:
        xrf[xrf <= 0.2*np.max(xrf)] = 0.
        xrf[xrf > 0.2*np.max(xrf)] = 1.

    #plt.figure()
    #plt.imshow(xrf)

    b = ndimage.measurements.center_of_mass(xrf)

    iy = np.int(np.round(b[0]))
    ix = np.int(np.round(b[1]))
    i_max = ix + iy * nx

    x_cen = x[i_max]
    y_cen = y[i_max]

    xrf_proj = np.sum(xrf,axis=1)
    xrf_proj_d = xrf_proj - np.roll(xrf_proj,1,0)
    i_tip = np.where(xrf_proj_d > 0.1)
    #print(i_tip[0][0])
    y_cen = y[(i_tip[0][0])*nx]

    #plt.figure()
    #plt.plot(xrf_proj)
    #plt.plot(xrf_proj_d)

    if moveflag:
        if x_motor == 'dssx':
            if piezomoveflag:
                print('move dssx to', x_cen)
                mov(smlld.dssx,x_cen)
            else:
                print('move dsx by', x_cen)
                movr(smlld.dsx, -1.*x_cen)
            sleep(.1)

        elif x_motor == 'dssz':
            if piezomoveflag:
                print('move dssz to', x_cen)
                mov(smlld.dssz,x_cen)
            else:
                print('move dsz by', x_cen)
                movr(smlld.dsz, x_cen)
            sleep(.1)

    if moveflag:
        print('y center', y_cen)
        if piezomoveflag:
            print('move dssy to:', y_cen)
            mov(smlld.dssy,y_cen)
        else:
            #movr(smlld.dsy, y_cen*0.001)
            print('move dsy by:', y_cen*0.001)
            movr(smlld.dsy, y_cen*0.001)
    else:
        print(x_cen,y_cen)

def tomo_scan_list(angle_list, x_start, x_end, x_num,
              y_start, y_end, y_num, exposure):
    x_0 = smlld.dssx.position
    z_0 = smlld.dssz.position
    y_0 = smlld.dssy.position
    th_0 = smlld.dsth.position
    dx_0 = smlld.dsx.position
    dz_0 = smlld.dsz.position
    dy_0 = smlld.dsy.position

    angle_list = np.array(angle_list)
    angle_num = np.size(angle_list)
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    y_num = np.int(y_num)
    exposure = np.float(exposure)
    offset_y = 0.32
    offset_x = 0
    #RE(fly2d(dssz, -8, 8, 80, dssy, -2, 2, 20, 0.05, return_speed=40))
    #RE(fly2d(dssz, -1, 1, 10, dssy, -1, 1, 10, 0.05, return_speed=40)) #dummy scan

    for i in range(angle_num):

        #while beamline_status.beam_current.get() <= 245:
        #    sleep(60)

        mov(smlld.dsth, angle_list[i])
        #x_start_real = x_start / np.cos(angle*np.pi/180.)
        #x_end_real = x_end / np.cos(angle*np.pi/180.)

        while (sclr2_ch4.get() < 20000):
            sleep(60)
            print('IC3 is lower than 20000, waiting...')

        if np.abs(angle_list[i]) <= 45:
            x_start_real = x_start / np.cos(angle_list[i] * np.pi / 180.)
            x_end_real = x_end / np.cos(angle_list[i] * np.pi / 180.)
            #x_range = (x_end - x_start) / np.cos(angle_list[i] * np.pi / 180.)
            #y_range = y_end - y_start
            RE(fly2d(smlld.dssx, -10, 10, 100, smlld.dssy,-2.5, 2.5, 25, 0.05, return_speed=40))
            #mov_to_image_cen_dsx(scan_id=-1, elem='Au_L', bitflag=0, moveflag=1,piezomoveflag=1)
            mov_to_image_center_tmp(scan_id=-1, elem='W_L', bitflag=1, moveflag=1,piezomoveflag=0)
            movr(smlld.dssy, offset_y)
            #RE(fly1d(smlld.dssx,-15,15,100,0.2))
            #mov_to_image_cen_corr_dsx(-1)
            #mov_to_line_center_mll(scan_id=-1,elem='Au_L',threshold=0.2,moveflag=1,movepiezoflag=1)


            sleep(1)
            RE(fly2d(smlld.dssx,x_start_real,x_end_real,x_num,smlld.dssy,
                     y_start, y_end, y_num, exposure, return_speed=40))
            #RE(fermat(smlld.dssx,smlld.dssy,x_range,y_range,0.05,1,exposure))
        else:
            x_start_real = x_start / np.abs(np.sin(angle_list[i] * np.pi / 180.))
            x_end_real = x_end / np.abs(np.sin(angle_list[i] * np.pi / 180.))
            #x_range = (x_end - x_start) / np.abs(np.sin(angle_list[i] * np.pi / 180.))
            #y_range = y_end - y_start
            RE(fly2d(smlld.dssz, -10, 10, 100, smlld.dssy,-2.5, 2.5, 25, 0.05, return_speed=40))
            #mov_to_image_cen_dsx(scan_id=-1, elem='Au_L', bitflag=0, moveflag=1,piezomoveflag=1)
            mov_to_image_center_tmp(scan_id=-1, elem='W_L', bitflag=1, moveflag=1,piezomoveflag=0)
            movr(smlld.dssy, offset_y)
            #RE(fly1d(smlld.dssz,-15,15,100,0.2))
            #mov_to_line_center_mll(scan_id=-1,elem='Au_L',threshold=0.2,moveflag=1,movepiezoflag=1)
            sleep(1)
            RE(fly2d(smlld.dssz,x_start_real,x_end_real,x_num, smlld.dssy,
                     y_start, y_end, y_num, exposure, return_speed = 40))
            #RE(fermat(smlld.dssz,smlld.dssy,x_range,y_range,0.05,1, exposure))
        #mov_to_image_cen_smar(-1)
        #mov_to_image_cen_dsx(-1)
        #merlin1.unstage()
        movr(smlld.dssy,-1*offset_y)
        print('waiting for 2 sec...')
        sleep(2)
    mov(smlld.dsth, th_0)
    mov(smlld.dssx, x_0)
    mov(smlld.dssy, y_0)
    mov(smlld.dssz, z_0)
    mov(smlld.dsx, dx_0)
    mov(smlld.dsy, dy_0)
    mov(smlld.dsz, dz_0)


def scan_translate(step_size,step_num, x_start, x_end, x_num,
              y_start, y_end, y_num, exposure):
    x_0 = zps.smarx.position
    z_0 = zps.smarz.position
    y_0 = zps.smary.position
    theta_0 = zps.zpsth.position

    step_size = np.float(step_size)
    step_num = np.int(step_num)
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    y_num = np.int(y_num)
    exposure = np.float(exposure)


    for i in range(step_num):

        while (sclr2_ch4.get() < 100000):
            sleep(60)
            print('IC3 is lower than 100000, waiting...')

        RE(fly1d(zps.zpssx,-10,10,100,0.5))
        mov_to_line_center(scan_id=-1,elem='Ba',threshold=0.2,moveflag=1,movepiezoflag=0)

        RE(fly2d(zps.zpssx,x_start,x_end,x_num,zps.zpssy,
                     y_start, y_end, y_num, exposure, return_speed=40))

        movr(zps.smary,step_size/1000.)

        merlin1.unstage()
        print('waiting for 2 sec...')
        sleep(2)

    mov(zps.smarx, x_0)
    mov(zps.smary, y_0)
    mov(zps.smarz, z_0)

def tomo_scan_list_zp(angle_list, x_start, x_end, x_num,
              y_start, y_end, y_num, exposure):
    x_0 = zps.smarx.position
    z_0 = zps.smarz.position
    y_0 = zps.smary.position
    theta_0 = zps.zpsth.position

    angle_list = np.array(angle_list)
    angle_num = np.size(angle_list)
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    y_num = np.int(y_num)
    exposure = np.float(exposure)

    for i in range(angle_num):

        mov(zps.zpsth, angle_list[i])

        while (sclr2_ch4.get() < 100000):
            sleep(60)
            print('IC3 is lower than 100000, waiting...')

        if np.abs(angle_list[i]) <= 45:

            x_start_real = x_start / np.abs(np.cos(angle_list[i] * np.pi / 180.))
            x_end_real = x_end / np.abs(np.cos(angle_list[i] * np.pi / 180.))
            print(x_start_real,x_end_real)
#            if angle_list[i] < -45:
#                x_start_real = (x_start+1.5) / np.cos(angle_list[i] * np.pi / 180.)
#                x_end_real = (x_end+1.5) / np.cos(angle_list[i] * np.pi / 180.)


#            RE(fly2d(zps.zpssx,-1.5,1.5,30,zps.zpssy,-1,1,20,0.1,return_speed=40))
#            mov_to_image_cen_smar(-1)

            RE(fly1d(zps.zpssx,-10,10,100,0.5))
            mov_to_line_center(scan_id=-1,elem='Ba',threshold=0.2,moveflag=1,movepiezoflag=0)

            RE(fly2d(zps.zpssx,x_start_real,x_end_real,x_num,zps.zpssy,
                     y_start, y_end, y_num, exposure, return_speed=40))

        else:
            x_start_real = x_start / np.abs(np.sin(angle_list[i] * np.pi / 180.))
            x_end_real = x_end / np.abs(np.sin(angle_list[i] * np.pi / 180.))
            print(x_start_real,x_end_real)
#            RE(fly2d(zps.zpssz,-1.5,1.5,30,zps.zpssy,-1,1,20,0.1,return_speed=40))
#            mov_to_image_cen_smar(-1)

            RE(fly1d(zps.zpssz,-10,10,100,0.5))
            mov_to_line_center(scan_id=-1,elem='Ba',threshold=0.2,moveflag=1,movepiezoflag=0)

            RE(fly2d(zps.zpssz,x_start_real,x_end_real,x_num, zps.zpssy,
                     y_start, y_end, y_num, exposure, return_speed = 40))

        merlin1.unstage()
        print('waiting for 2 sec...')
        sleep(2)


    mov(zps.zpsth, theta_0)
    mov(zps.smarx, x_0)
    mov(zps.smary, y_0)
    mov(zps.smarz, z_0)




def tomo_scan_list_zp_no_move(angle_list, x_start, x_end, x_num,
              y_start, y_end, y_num, exposure):
    x_0 = zps.smarx.position
    z_0 = zps.smarz.position
    y_0 = zps.smary.position

    angle_list = np.array(angle_list)
    angle_num = np.size(angle_list)
    x_start = np.float(x_start)
    x_end = np.float(x_end)
    x_num = np.int(x_num)
    y_start = np.float(y_start)
    y_end = np.float(y_end)
    y_num = np.int(y_num)
    exposure = np.float(exposure)


    for i in range(angle_num):

        mov(zps.zpsth, angle_list[i])

        while (sclr2_ch4.get() < 100000):
            sleep(60)
            print('IC3 is lower than 100000, waiting...')

        if np.abs(angle_list[i]) <= 45:
            x_start_real = x_start / np.cos(angle_list[i] * np.pi / 180.)
            x_end_real = x_end / np.cos(angle_list[i] * np.pi / 180.)

            RE(fly2d(zps.zpssx,-5,5,40,zps.zpssy,-4,4,40,0.05,return_speed=40))
#            mov_to_image_cen_smar(-1)
            cen = calc_image_cen_smar(-1,'Ca')
            xn_s = cen[0] - (x_end_real - x_start_real)/2
            xn_f = cen[0] + (x_end_real - x_start_real)/2
            yn_s = cen[1] - (y_end - y_start)/2
            yn_f = cen[1] + (y_end - y_start)/2
            print('x_cen=',cen[0])
            print('y_cen=',cen[1])
            print('xn_s=',xn_s)
            print('xn_f=',xn_f)
            print('yn_s=',y_start+cen[1])
            print('yn_f=',y_end+cen[1])

            RE(fly2d(zps.zpssx,xn_s,xn_f,x_num,zps.zpssy,yn_s, yn_f, y_num, exposure, return_speed=40))

        else:
            x_start_real = x_start / np.abs(np.sin(angle_list[i] * np.pi / 180.))
            x_end_real = x_end / np.abs(np.sin(angle_list[i] * np.pi / 180.))
            RE(fly2d(zps.zpssz,-5,5,40,zps.zpssy,-4,4,40,0.05,return_speed=40))
#            mov_to_image_cen_smar(-1)
            cen = calc_image_cen_smar(-1,'Ca')
            xn_s = cen[0] - (x_end_real - x_start_real)/2
            xn_f = cen[0] + (x_end_real - x_start_real)/2
            yn_s = cen[1] - (y_end - y_start)/2
            yn_f = cen[1] + (y_end - y_start)/2
            print('x_cen=',cen[0])
            print('y_cen=',cen[1])
            print('xn_s=',xn_s)
            print('xn_f=',xn_f)
            print('yn_s=',y_start+cen[1])
            print('yn_f=',y_end+cen[1])

            RE(fly2d(zps.zpssz,xn_s,xn_f,x_num, zps.zpssy,yn_s, yn_f, y_num, exposure, return_speed = 40))

        merlin1.unstage()
        print('waiting for 2 sec...')
        sleep(2)
    mov(zps.zpsth, angle_list[0])
    mov(zps.smarx, x_0)
    mov(zps.smary, y_0)
    mov(zps.smarz, z_0)





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
        mov(smlld.dsth, angle)
        #x_start_real = x_start / np.cos(angle*np.pi/180.)
        #x_end_real = x_end / np.cos(angle*np.pi/180.)

        '''
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
        '''
        while (sclr2_ch4.get() < 150000):
            sleep(60)
            print('IC3 is lower than 150000, waiting...')

        if np.abs(angle) <= 45:
            x_start_real = x_start / np.cos(angle * np.pi / 180.)
            x_end_real = x_end / np.cos(angle * np.pi / 180.)
            #RE(fly2d(zpssx, x_start_real, x_end_real, x_num, zpssy,
            #         y_start, y_end, y_num, exposure, return_speed=40))
            RE(fly2d(smlld.dssz,x_start_real,x_end_real,x_num,smlld.dssy,
                     y_start, y_end, y_num, exposure, return_speed=40))

        else:
            x_start_real = x_start / np.abs(np.sin(angle * np.pi / 180.))
            x_end_real = x_end / np.abs(np.sin(angle * np.pi / 180.))
            #RE(fly2d(zpssz, x_start_real, x_end_real, x_num, zpssy,
            #         y_start, y_end, y_num, exposure, return_speed=40))
            Re(fly2d(smlld.dssx,x_start_real,x_end_real,x_num, smlld.dssy,
                     y_start, y_end, y_num, exposure, return_speed = 40))

        #mov_to_image_cen_smar(-1)
        mov_to_image_cen_dsx(-1)
        merlin1.unstage()
        print('waiting for 5 sec...')
        sleep(5)
    #mov(zps.zpsth, 0)


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
    #movr(zp.zpx, dz * 3.75)
    movr(zp.zpy, -dz*0.003091258)
    movr(zp.zpx, (dz*0.003/40.33))


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


def move_fly_center(elem):
    scan_id, df = _load_scan(-1, fill_events=False)
    hdr = db[scan_id]['start']
    if elem in df:
        roi_data = np.asarray(df[elem])
    else:
        channels = [1, 2, 3]
        roi_keys = ['Det%d_%s' % (chan, elem) for chan in channels]
        for key in roi_keys:
            if key not in df:
                raise KeyError('ROI %s not found' % (key, ))
        roi_data = np.sum([getattr(df, roi) for roi in roi_keys], axis=0)

    scanned_axis = hdr['motor1']
    x = np.asarray(df[scanned_axis])
    nx = hdr['num1']
    ny = hdr['num2']
    roi_data = roi_data.reshape(ny,nx)
    x = x.reshape(ny,nx)
    ix,iy = ndimage.measurements.center_of_mass(roi_data)
    ix = np.int(ix)
    iy = np.int(iy)
    #i_max = find_mass_center(roi_data)
    #i_max = np.int(i_max)
    #print(ix,iy,x[ix,iy])
    #print(x)
    print('mass center:', x[ix,iy])
    #i_max = np.where(roi_data == np.max(roi_data))
    #mov(eval(scanned_axis),x[i_max[0]][0])

def th_fly2d_mll(th_start, th_end, num, x_start, x_end, x_num, y_start, y_end, y_num, sec):
    shutter('open')
    th_step = (th_end - th_start) / num
    movr(smlld.dsth, th_start)

    for i in range(num + 1):
        sleep(5)
        #RE(fly1d(zpssx,-5,5,100,0.1))
        #move_fly_center('Ge')
        RE(fly2d(dssx, x_start, x_end, x_num, dssy, y_start, y_end, y_num, sec, return_speed=40))
        movr(smlld.dsth, th_step)
    movr(smlld.dsth, -(th_end + th_step))
    shutter('close')


def th_fly2d(th_start, th_end, num, x_start, x_end, x_num, y_start, y_end,
             y_num, sec):
    shutter('open')
    th_step = (th_end - th_start) / num
    movr(zps.zpsth, th_start)

    for i in range(num + 1):
        sleep(5)
        #RE(fly1d(zpssx,-5,5,100,0.1))
        #move_fly_center('Ge')
        RE(fly2d(zpssx, x_start, x_end, x_num, zpssy, y_start, y_end, y_num, sec, return_speed=40))
        movr(zps.zpsth, th_step)
    movr(zps.zpsth, -(th_end + th_step))
    shutter('close')


def mov_diff(gamma, delta, r=500, calc=0):
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
    if x_yaw > 787 or x_yaw < -200:
        print('diff_x = ', -x_yaw,
              ' out of range, move diff_z upstream and try again')
    elif dz < -250 or dz > 0:
        print('diff_cz = ', dz,
              ' out of range, move diff_z up or down stream and try again')
    elif y1 > 750:
        print('diff_y1 = ', y1, ' out of range, move diff_z upstream '
              'and try again')
    elif y2 > 1000:
        print('diff_y2 = ', y2, ' out of range, move diff_z upstream '
              'and try again')
    else:
        print('diff_x = ', -x_yaw, ' diff_cz = ', dz,
              ' diff_y1 = ', y1, ' diff_y2 = ', y2)
        if calc == 0:
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
            while (diff.x.moving is True or diff.y1.moving is True or diff.y2.moving is True or diff.yaw.moving is True):
                sleep(2)
        else:
            print('Calculation mode; no motor will be moved')


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

def movr_zpsz_new(dist):
    movr(zps.zpsz,dist)
    movr(zps.smarx,dist*5.4161/1000.)
    movr(zps.smary,dist*1.8905/1000.)
def movr_smarz(dist):
    movr(zps.smarz, dist)
    movr(zps.smarx, dist*5.4161/1000.)
    movr(zps.smary, dist*1.8905/1000.)
def mll_movr_samp(angle, dx, dz):
    angle = angle*np.pi/180.0
    delta_x = (dx*np.cos(angle) - dz*np.sin(angle))
    delta_z = (dx*np.sin(angle) + dz*np.cos(angle))
    movr(smlld.dsx,delta_x)
    movr(smlld.dsz,delta_z)
def mll_movr_lab(dx, dz):
    angle = smlld.dsth.position
    angle = angle*np.pi/180.0
    delta_x = dx*np.cos(angle) - dz*np.sin(angle)
    delta_z = dx*np.sin(angle) + dz*np.cos(angle)
    movr(smlld.dsx, delta_x)
    movr(smlld.dsz, delta_z)

def mll_movr_samp_test(angle_offset,dist):
    angle_offset = -1 * angle_offset
    angle = smlld.dsth.position

    if np.abs(angle) <= 45.:
        alpha = (90 - angle - angle_offset) * np.pi / 180.
        beta = (90 + angle) * np.pi / 180.
        delta_x = -1 * np.sin(beta) * dist / np.sin(alpha) * np.cos(angle_offset)
        delta_z = np.sin(beta) * dist / np.sin(alpha) * np.sin(angle_offset)
    else:
        alpha = -1*(90 - angle - angle_offset) * np.pi / 180.
        beta = (180 - angle) * np.pi / 180.

        delta_x = -1 * np.sin(beta) * dist / np.sin(alpha) * np.cos(angle_offset)
        delta_z = np.sin(beta) * dist / np.sin(alpha) * np.sin(angle_offset)

    movr(smlld.dsx, delta_x)
    movr(smlld.dsz, delta_z)
def zp_movr_samp(th, dx, dz):
    th = th*np.pi/180.0
    delta_x = dx*np.cos(th) - dz*np.sin(th)
    delta_z = dx*np.sin(th) + dz*np.cos(th)
    movr(smarx, delta_x)
    movr(smarz, delta_z)


def zp_movr_lab(dx, dz):
    angle = zps.zpsth.position
    angle = angle*np.pi/180.0
    #angle = 14.2*np.pi/180
    delta_x = dx*np.cos(angle) - dz*np.sin(angle)
    delta_z = dx*np.sin(angle) + dz*np.cos(angle)
    movr(smarx, delta_x)
    movr(smarz, delta_z)


def discharge_scan():
    go_to_energy(7.131)
    sleep(1)
    RE(fly2d(zpssx, -10, 10, 100, zpssy, -10, 10, 100, 0.025, return_speed=40))
#    sleep(1)
#    merlin1.unstage()
    go_to_energy(7.151)
    sleep(1)
    RE(fly2d(zpssx, -10, 10, 100, zpssy, -10, 10, 100, 0.025, return_speed=40))
#    sleep(1)
#    merlin1.unstage()
    go_to_energy(7.131)

def go_to_energy(energy_kev=7.13):
    current_bragg = dcm.th.position
    current_energy = 12.39842 / (2.*3.1355893*np.sin(current_bragg*np.pi/180.))

    bragg = np.arcsin(12.39842/(2.*3.1355893*energy_kev)) * 180. / np.pi
    mov(dcm.th,bragg)

    current_ugap = ugap.position
    ugap_value = 7.61 + (energy_kev - 7.114)
    if np.abs(ugap_value-current_ugap) > 0.005:
        mov(ugap,ugap_value)

    dz = (energy_kev - current_energy) * 14.1129
    movr_zpsz_new(dz)


def multi_pos_scan(scan_list,
                   x_range_list, x_num_list,
                   y_range_list, y_num_list,
                   exp_list):
    """
    Parameters
    ----------
    scan_list : list
         list of pre-idetifed locations as scan_id

    x_range_list, y_range_list : list
         list of scanwidth, scan +/- value/2 for each point

    x_num_list, y_num_list : list
         Number of points at each scan location in

    exp_list : list
         Exposure time per location
    """
    for i, (scan, x_range, x_num, y_range, y_num, exposure) in enumerate(
            zip(scan_list, x_range_list, x_num_list, y_range_list, y_num_list, exp_list)):
        print('scan ', i, ' move to #', scan, 'position')
        # TODO make this a plan
        recover_mll_scan_pos(int(scan))
        sleep(0.5)
        RE(fly2d(dssx,
                 -x_range/2, x_range/2, x_num,
                 dssy,
                 -y_range/2, y_range/2, y_num,
                 exposure))
        sleep(0.5)

def multi_pos_scan_plan(scan_list,
                        x_range_list, x_num_list,
                        y_range_list, y_num_list,
                        exp_list):
    """
    Parameters
    ----------
    scan_list : list
         list of pre-idetifed locations as scan_id

    x_range_list, y_range_list : list
         list of scanwidth, scan +/- value/2 for each point

    x_num_list, y_num_list : list
         Number of points at each scan location in

    exp_list : list
         Exposure time per location
    """
    for i, (scan, x_range, x_num, y_range, y_num, exposure) in enumerate(
            zip(scan_list, x_range_list, x_num_list, y_range_list, y_num_list, exp_list)):
        print('scan ', i, ' move to #', scan, 'position')
        yield from recover_mll_scan_pos_plan(int(scan))
        yield from fly2d(dssx,
                         -x_range/2, x_range/2, x_num,
                         dssy,
                         -y_range/2, y_range/2, y_num,
                         exposure)


def multi_pos_xanes():
    x_list = np.array([-2.836,-2.926,-2.925,-2.969])
    y_list = np.array([4.0728,4.098,4.1175,4.103])
    num_p = np.size(x_list)

    x_range_list = np.array([6,5,6,3])
    y_range_list = np.array([6,3,4,3])

    x_num_list = np.round(2*x_range_list/0.5)
    y_num_list = np.round(2*y_range_list/0.5)

    #e_list = np.array([7.11,7.1106,7.1111,7.1116,7.1121,7.113,7.114,7.115,7.116,7.118,7.12,7.125,7.13])
    #exp_list = np.array([0.1,0.1,0.1,0.1,0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05])

    e_list = np.array([7.115,7.116])
    exp_list = np.array([0.02,0.02])
    x_zero = 0

    for i in range(num_p):
        mov(zps.smarx,x_list[i])
        mov(zps.smary,y_list[i])
        sleep(1)
        print(-1*x_range_list[i], x_range_list[i], x_num_list[i], -1*y_range_list[i], y_range_list[i], y_num_list[i])
        xanes_scan(e_list, x_zero, x_zero, -1*x_range_list[i], x_range_list[i], x_num_list[i], -1*y_range_list[i], y_range_list[i], y_num_list[i], exp_list)

def test_scan():
    e_list = np.array([7.11,7.1106,7.1111,7.1116,7.1121,7.113,7.114,7.115,7.116,7.118,7.12,7.125,7.13])
    exp_list = np.array([0.1,0.1,0.1,0.1,0.1, 0.1, 0.1, 0.1, 0.1, 0.05, 0.05, 0.05, 0.05])
    print(e_list.size,exp_list.size)

    mov(zps.smarx,-2.836)
    mov(zps.smary,4.0728)
    sleep(1)
    xanes_scan(e_list,0,0,-6,6,120,-6,6,120,exp_list)


    mov(zps.smarx,-2.926)
    mov(zps.smary,4.098)
    sleep(1)
    xanes_scan(e_list,0,0,-5,5,100,-3,3,60,exp_list)

    mov(zps.smarx,-2.925)
    mov(zps.smary,4.1175)
    sleep(1)
    xanes_scan(e_list,0,0,-6,6,120,-4,4,80,exp_list)

    mov(zps.smarx,-2.969)
    mov(zps.smary,4.103)
    sleep(1)
    xanes_scan(e_list,0,0,-3,3,60,-3,3,60,exp_list)

    '''
    x_start = -6
    x_end = 6
    x_num = 24.
    y_start = -6
    y_end = 6
    y_num =24.
    exp = np.array([0.02,0.02])
    RE(fly2d(zps.zpssx,x_start,x_end,x_num,zps.zpssy,y_start,y_end,y_num,exp[1]))
    '''

def xanes_scan_sim(x_start,x_end,x_num,y_start,y_end,y_num,exposure,scan_flag = True):
    energy_list = np.array([16.118,16.135])  # keV unit
    bragg_0 = dcm.th.position
    ugap_0 = ugap.position
    dsx_0 = smlld.dsx.position
    dsy_0 = smlld.dsy.position
    sbx_0 = smlld.sbx.position
    sbz_0 = smlld.sbz.position

    factor = 350 # um / keV

    energy_0 = 16.118#12.39842 / (2.*3.1355893 * np.sin(bragg_0 * np.pi / 180.))
    bragg_0 = np.arcsin(12.39842/(2.*3.1355893*energy_0)) * 180. / np.pi
    bragg_list =np.arcsin(12.39842/(2.*3.1355893*energy_list)) * 180. / np.pi
    num_point = np.size(energy_list)

    print('energy list:', energy_list)
    print('bragg list:',bragg_list)
    print('current bragg angle',bragg_0)
    print('current sbz',sbz_0)
    #print(factor*(energy_list[1] - energy_0))
    if scan_flag:
        for i in range(num_point):
            print('move to energy ', energy_list[i],'keV')
            mov(dcm.th,bragg_list[i])
            if i == 0:
                print('move ugap to 7.419 for', energy_list[i], 'keV')
                mov(ugap,7.419)
                sleep(1)
                print('move sbz by', factor*(energy_list[i] - energy_0))
                movr(smlld.sbz, factor*(energy_list[i] - energy_0))
            else:
                print('move ugap to 7.429 for ',energy_list[i], 'keV')
                mov(ugap,7.429)
                sleep(1)
                print('move sbz by', factor*(energy_list[i] - energy_list[i-1]))
                movr(smlld.sbz, factor*(energy_list[i] - energy_list[i-1]))
            RE(fly2d(smlld.dssx,x_start,x_end,x_num,smlld.dssy,y_start,y_end,y_num,exposure))
        mov(dcm.th,bragg_0)
        mov(smlld.sbz,sbz_0)
        mov(ugap,ugap_0)

def xanes_scan(energy_list,x_start,x_end,x_num,y_start,y_end,y_num,exposure,peak_flag=0,sign='max',elem='Fe',printflag=True):


    # fit ugap curve
#    x = [9.6482,9.6532,9.6582,9.6687,9.6706,9.6737,9.6752,9.6817,9.7022]
#    y = [6.462,6.462,6.465,6.468,6.47,6.472,6.472,6.478,6.488]
    x=[7.05, 7.1, 7.12, 7.13, 7.142, 7.15, 7.2, 7.25] # for Fe edge
    y=[7.54448, 7.5894, 7.6073, 7.6163, 7.6271, 7.6343, 7.6791, 7.7240]
    fit_para = np.polyfit(x,y,1)
    fit_func = np.poly1d(fit_para)

    energy_list = np.array(energy_list) # unit keV
    bragg_list = np.arcsin(12.39842/(2.*3.1355893*energy_list)) * 180. / np.pi
    num_bragg = np.size(bragg_list)
    current_det = gs.PLOT_Y
    #gs.PLOT_Y= elem

    exposure = np.float(exposure)

    start_bragg = dcm.th.position
    start_zpz1 = zp.zpz1.position
    start_zpx = zp.zpx.position
    start_zpy = zp.zpy.position
    start_smarx = zps.smarx.position
    start_smary = zps.smary.position
    start_gap = ugap.position
    start_energy = 12.39842 / (2.*3.1355893 * np.sin(start_bragg * np.pi / 180.))

    for i in range(num_bragg):
        current_bragg = dcm.th.position
        current_energy = 12.39842 / (2.*3.1355893 * np.sin(current_bragg * np.pi / 180.))
        mov(dcm.th,bragg_list[i])

        new_ugap = fit_func(energy_list[i])
        print('New ugap:'+np.str(new_ugap))
        #mov(ugap, new_ugap)


        #current_ugap = 7.615 + (energy_list[i] - 7.114) * 0.035 / 0.035
    #    current_ugap = 8.795 + (energy_list[i] - 8.34)
        gap_tmp = ugap.position
        #print(energy_list[i],current_ugap,gap_tmp,exposure[i])

        if (np.abs(new_ugap - gap_tmp) > 0.005):
            mov(ugap,new_ugap)
            sleep(10)

        delta_kev = energy_list[i] - current_energy
        dist_zpz1 = -1*delta_kev*7.7419
        movr_zpz1(dist_zpz1)


        #movr(zps.zpsz,delta_kev*14.11290323)
        #movr(zps.smarx,delta_kev*14.11290323*5.4161/1000.)
        #movr(zps.smary,delta_kev*14.11290323*1.8905/1000.)
        if peak_flag:
            peak_ic()
        '''
        RE(dscan(dcm.rf,-1, 1, 40, 1))
        if sign == 'max':
            mov(dcm.rf,gs.PS.max[0])
        elif sign == 'cen':
            mov(dcm.rf,gs.PS.cen)

        RE(dscan(m2.pf,-1, 1, 40, 1))
        if sign == 'max':
            mov(m2.pf,gs.PS.max[0])
        elif sign == 'cen':
            mov(m2.pf,gs.PS.cen)
        '''
        print(x_start,x_end,x_num,y_start,y_end,y_num,exposure)
        RE(fly2d(zps.zpssx,x_start,x_end,x_num,zps.zpssy,y_start,y_end,y_num,exposure))
        scan_id,df = _load_scan(-1,fill_events=False)
        if printflag:
            plot2dfly(-1,elem,norm='sclr1_ch4')
            plt.title('#'+np.str(scan_id)+' '+elem+', at'+np.str(energy_list[i])+' keV')
            printfig()
#        RE(fly2d(zps.zpssx, x_start, x_end, x_num, zps.zpssy, y_start, y_end, y_num, exposure[i], return_speed=50))
#        RE(fly2d(zps.zpssx, x_start+x_list[i], x_end+x_list[i], x_num, zps.zpssy, y_start+y_list[i], y_end+y_list[i], y_num, exposure, return_speed=50))

#        sleep(1)
        #merlin1.unstage()
    gs.PLOT_Y = current_det
    mov(dcm.th,start_bragg)
    mov(zp.zpz1,start_zpz1)
    mov(zp.zpx,start_zpx)
    mov(zp.zpy,start_zpy)
    mov(zps.smarx,start_smarx)
    mov(zps.smary,start_smary)
    mov(ugap,start_gap)

def xanes_scan_bp(bragg_list,x_start,x_end,x_num,y_start,y_end,y_num,exposure,sign='max'):
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

    movr(sbx, dx)
    movr(sbz, dz)

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

    c_ssz = smll.dssz.position
    c_ssy = smll.dssy.position
    c_ssx = smll.dssx.position

    print('Current ssx = %.3f' % c_ssx)
    print('Current ssy = %.3f' % c_ssy)
    print('Current ssz = %.3f' % c_ssz)

    smll_kill_piezos()

    t_ssz = c_ssz + dist*np.cos(alpha)
    dz = t_ssz - smll.dssz.position
    dy = c_ssy - smll.dssy.position

    movr(sbz, dz)
    movr(dsy, dy/1000.0)

    sleep(5)

    smll_sync_piezos()
    mov(dssy, c_ssy)
    mov(dssz, t_ssz)

    print('post-move x = %.3f' % smll.dssx.position)
    print('Post-move y = %.3f' % smll.dssy.position)
    print('Post-move z = %.3f' % smll.dssz.position)

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

    now = datetime.now()
    fn = '/data/motor_positions/log-'+np.str(now.year)+'-'+np.str(now.month)+'-'+np.str(now.day)+'-'+np.str(now.hour)+'-'+np.str(now.minute)+'.log'
    f = open(fn,'w')
    original = sys.stdout
    sys.stdout = Tee(sys.stdout, f)
    wh_pos()
    sys.stdout = original
    f.close()
    if print_flag:
        shutil.copyfile(fn,'/data/motor_positions/tmp.log')
        os.system("lp -o cpi=20 -o lpi=8 -o media='letter' -d HXN-printer-1 /data/motor_positions/tmp.log")


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


def plot_fermat(scan_id,elem='Ga',norm=1):
    df = db.get_table(db[scan_id],fill=False)
    x = np.asarray(df.zpssx)
    y = np.asarray(df.zpssy)
    io = np.asfarray(df.sclr1_ch4)
    #if elem == 'Ga':
    xrf = np.asfarray(eval('df.Det1_'+elem)) + np.asfarray(eval('df.Det2_'+elem)) + np.asfarray(eval('df.Det3_'+elem))
    #elif elem == 'K':
    #    xrf = np.asfarray(df.Det1_K) + np.asfarray(df.Det2_K) + np.asfarray(df.Det3_K)

    if norm:
        xrf /= (io+1.e-8)
        #print(xrf.dtype)
    props = dict(alpha=0.8, edgecolors='none' )
    plt.figure()
    plt.scatter(x,y,c=xrf,s=50,marker='s',**props)
    plt.xlim([np.min(x),np.max(x)])
    plt.ylim([np.min(y),np.max(y)])
    plt.title('scan '+ np.str(scan_id))
    plt.axes().set_aspect('equal')
    plt.gca().invert_yaxis()
    plt.colorbar()

    plt.show()


def mov_to_image_cen_zpsx(scan_id=-1, elem='Ni', bitflag=1):

    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    x = np.asarray(df2.zpssx)
    y = np.asarray(df2.zpssy)
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    if bitflag:
        xrf[xrf <= 0.25*np.max(xrf)] = 0.
        xrf[xrf > 0.25*np.max(xrf)] = 1.

    b = ndimage.measurements.center_of_mass(xrf)

    iy = np.int(np.round(b[0]))
    ix = np.int(np.round(b[1]))
    i_max = ix + iy * nx

    x_cen = x[i_max]
    y_cen = y[i_max]
    print('move zpsx by', x_cen)
    #print('move zpssx, zpssy to ',0, 0)

    #movr(zps.zpsx, x_cen*0.001)
    #mov(zps.zpssx,0)
    #mov(zps.zpssy,0)
    sleep(.1)


def retreat_xrf_roi(scan_id = -1, elem='Au', bitflag=1):
    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))
    max_xrf=np.max(xrf)

    if bitflag:
        xrf[xrf<(0.2*max_xrf)] = 0
        #xrf[xrf>=(0.2*max_xrf)] = 1

    return xrf


def mov_to_image_cen_corr_dsx(scan_id=-1, elem='Pt',bitflag=1, moveflag=1):
    print(scan_id)
    image_ref = retreat_xrf_roi(scan_id-2, elem,bitflag)
    image = retreat_xrf_roi(scan_id,elem,bitflag)
    corr = signal.correlate2d(image_ref, image, boundary='symm', mode='same')
    #nx,ny = np.shape(image)
    max_y,max_x = np.where(corr == np.max(corr))

    df2 = db.get_table(db[scan_id],fill=False)
    hdr = db[scan_id]['start']
    x_motor = hdr['motor1']
    y_motor = hdr['motor2']
    x = np.asarray(df2[x_motor])
    y = np.asarray(df2[y_motor])

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    step_x_um = (hdr['scan_end1'] - hdr['scan_start1']) / nx
    step_y_um = (hdr['scan_end2'] - hdr['scan_start2']) / ny


    dx_um = -1*(max_x - nx/2) * step_x_um
    dy_um = -1*(max_y - ny/2) * step_y_um

    if x_motor == 'dssx':
        print('move dsx by', -dx_um)
        if moveflag:
            if np.abs(dx_um)>step_x_um:
                movr(smlld.dsx, -dx_um)

    if x_motor == 'dssz':
        print('move dsz by', dx_um)
        if moveflag:
            if np.abs(dx_um)>step_x_um:
                movr(smlld.dsz, dx_um)

    '''
    image_ref_crop = image_ref[:,nx/4:nx*3/4]
    image_crop = image[:,nx/4:nx*3/4]
    corr_crop = signal.correlate2d(image_ref_crop, image_crop, boundary='symm', mode='same')
    max_y_crop,max_x_crop = np.where(corr_crop == np.max(corr_crop))
    dy_um_crop = -1*(max_y_crop - ny/2) * step_y_um
    plt.figure()
    plt.subplot(221)
    plt.imshow(image_ref_crop)
    plt.subplot(222)
    plt.imshow(image_crop)
    plt.subplot(223)
    plt.imshow(corr_crop)
    plt.show()
    '''

    ly_ref = np.sum(image_ref,axis=1)
    ly = np.sum(image,axis=1)
    corr_1d = np.correlate(np.squeeze(ly_ref),np.squeeze(ly),'same')
    max_y= np.where(corr_1d == np.max(corr_1d))
    #print(max_y)
    dy_um = -1*(max_y[0] - ny/2) * step_y_um
    '''
    plt.figure()
    plt.subplot(311)
    plt.plot(ly_ref)
    plt.subplot(312)
    plt.plot(ly)
    plt.subplot(313)
    plt.plot(corr_1d)
    plt.show()
    '''
    print('move y by', dy_um*0.001)
    if moveflag:
        #if np.abs(dy_um)>step_y_um:
        movr(smlld.dsy, dy_um*0.001)



def mov_to_image_cen_dsx(scan_id=-1, elem='Au', bitflag=1, moveflag=1,piezomoveflag=1,x_offset=0,y_offset=0):

    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    #xrf_Pt = np.asfarray(eval('df2.Det2_' + 'Ni')) + np.asfarray(eval('df2.Det1_' + 'Ni')) + np.asfarray(eval('df2.Det3_' + 'Ni'))
    hdr = db[scan_id]['start']
    x_motor = hdr['motor1']
    y_motor = hdr['motor2']
    x = np.asarray(df2[x_motor])
    y = np.asarray(df2[y_motor])
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    '''
    xrf_Pt = xrf_Pt/I0
    index_Pt = np.where(xrf_Pt == np.max(xrf_Pt))
    Pt_max_y = y[index_Pt]
    Pt_max_x = x[index_Pt]
    Pt_max_target_y = -0.5
    Pt_max_target_x = 0
    '''
    #x_offset = 0
    #y_offset = -2.5
 #   print('move dsy by', Pt_max_target_y-Pt_max_y)

 #   movr(smlld.dsy,-1.*(Pt_max_target_y-Pt_max_y)*1.e-3)


    #plt.figure()
    #plt.imshow(xrf)

    if bitflag:
        xrf[xrf <= 0.2*np.max(xrf)] = 0.
        xrf[xrf > 0.2*np.max(xrf)] = 1.

    b = ndimage.measurements.center_of_mass(xrf)

    iy = np.int(np.round(b[0]))
    ix = np.int(np.round(b[1]))
    i_max = ix + iy * nx

    x_cen = x[i_max]
    y_cen = y[i_max]

    print(b,ix,iy,i_max,x_cen,y_cen)
    # if moveflag:
        # movr(smlld.dsy,y_cen/1000.)

    if x_motor == 'dssx':
        # print('move dsz,by', Pt_max_x-Pt_max_target_x, 'um')
        #print('move dsx by', 1*(x_cen - target_x))
        if moveflag:
            if piezomoveflag:
                print('move dssx to', (x_cen+x_offset))
                mov(smlld.dssx,(x_cen+x_offset))
            else:
                movr(smlld.dsx, -1.*(x_cen+x_offset))
            #movr(smlld.dsz, Pt_max_x-Pt_max_target_x)
            #mov(zps.zpssx,0)
        sleep(.1)

    elif x_motor == 'dssz':
        #print('move dsx,by', -1*(Pt_max_x-Pt_max_target_x) , 'um')
        print('x center ',x_cen)
        #print('move dsz by', (x_cen+x_offset))
        if moveflag:
            if piezomoveflag:
                print('move dssz to', (x_cen+x_offset))
                mov(smlld.dssz,(x_cen + x_offset))
            else:
                print('move dsz by', (x_cen+x_offset))
                movr(smlld.dsz, (x_cen + x_offset))
#        movr(smlld.dsx, -1*(Pt_max_x-Pt_max_target_x))
        #mov(zps.zpssx,0)
    sleep(.1)

    if moveflag:
        print('y center', y_cen)
        if piezomoveflag:
            print('move dssy to:', (y_cen +y_offset)*0.001)
            mov(smlld.dssy,(y_cen - y_offset))
        else:
            #movr(smlld.dsy, y_cen*0.001)
            print('move dsy by:', (y_cen +y_offset)*0.001)
            movr(smlld.dsy, (y_cen + y_offset)*0.001)
    #mov(zps.zpssy,0)
    sleep(.1)



def calc_image_cen_smar(scan_id=-1, elem='Er', bitflag=1, movflag=1):

    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    hdr = db[scan_id]['start']
    x_motor = hdr['motor1']
    y_motor = hdr['motor2']
    x = np.asarray(df2[x_motor])
    y = np.asarray(df2[y_motor])
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0

    if bitflag:
        coe=0.6
        xrf[xrf <= coe*np.max(xrf)] = 0.
        xrf[xrf > coe*np.max(xrf)] = 1.
        xrf = np.asarray(np.reshape(xrf,(ny,nx)))


    else:
        xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    b = ndimage.measurements.center_of_mass(xrf)

    iy = np.int(np.round(b[0]))
    ix = np.int(np.round(b[1]))
    i_max = ix + iy * nx

    x_cen = x[i_max]
    y_cen = y[i_max]

    cen=[x_cen,y_cen]
    return cen




def mov_to_image_cen_smar(scan_id=-1, elem='Er', bitflag=1, movflag=1):

    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    hdr = db[scan_id]['start']
    x_motor = hdr['motor1']
    y_motor = hdr['motor2']
    x = np.asarray(df2[x_motor])
    y = np.asarray(df2[y_motor])
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0
#    xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    #plt.figure()
    #plt.imshow(xrf)

    if bitflag:
        coe=0.6
        xrf[xrf <= coe*np.max(xrf)] = 0.
        xrf[xrf > coe*np.max(xrf)] = 1.
        xrf = np.asarray(np.reshape(xrf,(ny,nx)))


    else:
        xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    b = ndimage.measurements.center_of_mass(xrf)

    iy = np.int(np.round(b[0]))
    ix = np.int(np.round(b[1]))
    i_max = ix + iy * nx

    x_cen = x[i_max]
    y_cen = y[i_max]



    #print(b,ix,iy,i_max,x_cen,y_cen)

    if x_motor == 'zpssx':
        print('move smarx,by', x_cen, 'um')
        #print('move zpssx, zpssy to ',0, 0)
        if movflag:
            movr(zps.smarx, x_cen*0.001)
        #mov(zps.zpssx,0)
        sleep(.1)

    elif x_motor == 'zpssz':
        print('move smarz,by', x_cen, 'um')
        if movflag:
            movr(zps.smarz, x_cen*0.001)
        #mov(zps.zpssx,0)
        sleep(.1)
#    print('move smary,by', y_cen, 'um')
#    if movflag:
#        movr(zps.smary, y_cen*0.001)
    #mov(zps.zpssy,0)
    sleep(.1)


def mov_to_line_center(scan_id=-1,elem='Ga',threshold=0,moveflag=1,movepiezoflag=0):
    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    hdr=db[scan_id]['start']
    x_motor = hdr['motor']
    x = np.asarray(df2[x_motor])
    xrf[xrf<(np.max(xrf)*threshold)] = 0.
    xrf[xrf>=(np.max(xrf)*threshold)] = 1.
    #print(x)
    #print(xrf)
    mc = find_mass_center_1d(xrf,x)
    print(mc)
    if moveflag:
        if x_motor == 'zpssx':
            #if((mc < .2) and movepiezoflag):
            mov(zps.zpssx,mc)
            #else:
            #    movr(zps.smarx,(mc)/1000.)
        if x_motor == 'zpssy':
            #if((mc < .2) and movepiezoflag):
            mov(zps.zpssy,mc)
            #else:
            #    movr(zps.smary,(mc-zps.zpssy.position)/1000.)
        if x_motor == 'zpssz':
            #if((mc < .2) and movepiezoflag):
            mov(zps.zpssz,mc)
            #else:
            #    movr(zps.smarz,(mc)/1000.)
    else:
        if x_motor == 'zpssx':
            print('move smarx by '+np.str(mc/1000.))
        if x_motor == 'zpssy':
            print('move smary by '+np.str(mc/1000.))


def mov_to_line_center_mll(scan_id=-1,elem='Au',threshold=0,moveflag=1,movepiezoflag=0):
    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    hdr=db[scan_id]['start']
    x_motor = hdr['motor']
    x = np.asarray(df2[x_motor])
    xrf[xrf<(np.max(xrf)*threshold)] = 0.
    xrf[xrf>=(np.max(xrf)*threshold)] = 1.
    #print(x)
    #print(xrf)
    mc = find_mass_center_1d(xrf,x)
    print(mc)
    if moveflag:
        if x_motor == 'dssx':
            if(movepiezoflag):
                mov(smlld.dssx,mc)
            else:
                movr(smlld.dsx,-1*mc)
        if x_motor == 'dssy':
            if(movepiezoflag):
                mov(smlld.dssy,mc)
            else:
                movr(smlld.dsy,mc/1000.)
        if x_motor == 'dssz':
            if(movepiezoflag):
                mov(smlld.dssz,mc)
            else:
                movr(smlld.dsz,mc)
    else:
        if x_motor == 'dssx':
            print('move dssx by '+np.str(mc))
        if x_motor == 'dssz':
            print('move dssz by '+np.str(mc))

def mov_to_image_cen_zpss(scan_id=-1, elem='Ni', bitflag=1):

    df2 = db.get_table(db[scan_id],fill=False)
    xrf = np.asfarray(eval('df2.Det2_' + elem)) + np.asfarray(eval('df2.Det1_' + elem)) + np.asfarray(eval('df2.Det3_' + elem))
    #x = np.asarray(df2.zpssx)
    x = np.asarray(df2.zpssz)
    y = np.asarray(df2.zpssy)
    I0 = np.asfarray(df2.sclr1_ch4)

    scan_info=db[scan_id]
    tmp = scan_info['start']
    nx=tmp['plan_args']['num1']
    ny=tmp['plan_args']['num2']

    xrf = xrf/I0
    xrf = np.asarray(np.reshape(xrf,(ny,nx)))

    if bitflag:
        xrf[xrf <= 0.25*np.max(xrf)] = 0.
        xrf[xrf > 0.25*np.max(xrf)] = 1.


    b = ndimage.measurements.center_of_mass(xrf)

    iy = np.int(np.round(b[0]))
    ix = np.int(np.round(b[1]))
    i_max = ix + iy * nx

    x_cen = x[i_max]
    y_cen = y[i_max]
    #print('move smarx, smary by', x_cen, y_cen)
    print('move zpssz, zpssy to ',x_cen, y_cen)

    #movr(zps.smarx, x_cen*0.001)
    mov(zps.zpssz,x_cen)
    sleep(.1)
    #movr(zps.smary, y_cen*0.001)
    mov(zps.zpssy,y_cen)
    sleep(.1)

def night_scan_nov():
    RE(fly2d(zpssx, -7, 6, 130, zpssy, -6.5, 6.5, 130, 0.1, return_speed=40))
    sleep(10)
    RE(fly2d(zpssx, -3, 3, 120, zpssy, -6, 1, 140, 0.2, return_speed=40))
    sleep(10)

    shutter('close')



def overnight_scan():
    RE(fly2d(dssx, -2.5, 2.5, 200, dssy, -1.5, 1.5, 120, 0.4))

    mov(smlld.dsx,2.02)
    mov(smlld.dsy,-0.00055)

    for i in range(3):
        RE(fly2d(dssx, -5, 5, 200, dssy, -5, 5, 200, 0.4))

    mov(smlld.dsx,2.)
    mov(smlld.dsy,0.0184)

def go_to_grid(grid='top'):
    if grid == 'top':
        mov(smlld.dsx,-2335.72)
        mov(smlld.dsz,1720.98)
        mov(smlld.dsy,0.56015)
        mov(smlld.sz,754.948)
    elif grid == 'bottom':
        mov(smlld.dsx,-2275.72)
        mov(smlld.dsz,2075.98)
        mov(smlld.dsy,2.2576)
        mov(smlld.sz,754.948)
    elif grid == 's1':
        mov(smlld.dsx,-2322.611)
        mov(smlld.dsz,1720.96)
        mov(smlld.dsy,0.5928)
        mov(smlld.sz,754.948)
def mll_tracking(dx, dy):
    th = smlld.dsth.position
    if np.abs(th) <= 45:
        movr(smlld.dsx, -dx)
    elif np.abs(th) > 45:
        movr(smlld.dsz, dx)
    movr(smlld.dsy, dy/1000.0)


def beep(rp=3):
    b = lambda x: os.system("echo -n '\a'; sleep 0.2;" * x)
    b(rp)

def scale_fly2d(x_start,x_end,x_num,y_start,y_end,y_num,exp):
    angle = smlld.dsth.position
    angle_rad = np.abs(angle * 3.14 / 180.)

    if np.abs(angle) <= 45.:
        x_start_new = x_start / np.cos(angle_rad)
        x_end_new = x_end / np.cos(angle_rad)

        print(angle,' deg', 'scan dssx', 'x scan range: ', x_start_new, '--', x_end_new)
        RE(fly2d(smlld.dssx,x_start_new,x_end_new,x_num,smlld.dssy,y_start,y_end,y_num,exp))
        beep()

    else:
        x_start_new = x_start / np.sin(angle_rad)
        x_end_new = x_end / np.sin(angle_rad)

        print(angle,' deg','scan dssz', 'x scan range: ', x_start_new, '--', x_end_new)
        RE(fly2d(smlld.dssz,x_start_new,x_end_new,x_num, smlld.dssy, y_start, y_end, y_num, exp))
        beep()


def over_night_scan_Sima():
    for i in range(20):
        print('scan #: ',i)
        RE(fly2d(zpssx, -2.5, 9.5, 130, zpssy, -8.5, 3, 115, 0.15))
        print('sleeping 2 sec')
        sleep(2)

def extract_mll_scan_pos(scan_id):
    data = db.get_table(db[scan_id], stream_name='baseline')
    dsx_pos = data.dsx[1]
    dsy_pos = data.dsy[1]
    dsz_pos = data.dsz[1]
    dsth_pos = data.dsth[1]
    sbx_pos = data.sbx[1]
    sbz_pos = data.sbz[1]
    dssx_pos = data.dssx[1]
    dssy_pos = data.dssy[1]
    dssz_pos = data.dssz[1]
    #fdet1_x = data.fdet1_x[1]
    #print(ssx,ssy,ssz)

    print('scan '+np.str(scan_id))
    print('dsx:',dsx_pos, ', dsy:',dsy_pos, ', dsz:',dsz_pos,', dsth:',dsth_pos)
    print('sbx:',sbx_pos, ', sbz:',sbz_pos)
    print('dssx:',dssx_pos,', dssy:',dssy_pos,', dssz:',dssz_pos)

    return {
        smlld.dsz: dsz_pos,
        smlld.dsx: dsx_pos,
        smlld.dsy: dsy_pos,
        smlld.dsth: dsth_pos,
        smlld.dssx: dssx_pos,
        smlld.dssy: dssy_pos,
        smlld.dssz: dssz_pos,
        smlld.sbx: sbx_pos,
        smlld.sbz: sbz_pos}


def recover_mll_scan_pos_plan(scan_id, base_moveflag=True):
    coarse_motors = [smlld.dsz, smlld.dsx, smlld.dsy, smlld.dsth]
    piezo_motors = [smlld.dssx, smlld.dssy, smlld.dssz]
    base_motors = [smlld.sbx, smlld.sbz]

    targets = extract_mll_scan_pos(scan_id)
    cur_postions = {}
    for m in coarse_motors + piezo_motors + base_motors:
        cur_pos[m] = (yield from bp.read(m))

    grp_name = 'recover_moves'

    for m in coarse_motors + piezo_motors:
        yield from bp.abs_set(m, targets[m], group=grp_name)
    if base_moveflag:
        for m in base_motors:
            yield from bp.abs_set(m, targets[m], group=grp_name)
    yield from bp.wait(grp_name)

    return cur_pos

def recover_mll_scan_pos(scan_id,moveflag=True,base_moveflag=True,det_moveflag=False):
    data = db.get_table(db[scan_id], stream_name='baseline')
    dsx_pos = data.dsx[1]
    dsy_pos = data.dsy[1]
    dsz_pos = data.dsz[1]
    dsth_pos = data.dsth[1]
    sbx_pos = data.sbx[1]
    sbz_pos = data.sbz[1]
    dssx_pos = data.dssx[1]
    dssy_pos = data.dssy[1]
    dssz_pos = data.dssz[1]
    #fdet1_x = data.fdet1_x[1]
    #print(ssx,ssy,ssz)

    print('scan '+np.str(scan_id))
    print('dsx:',dsx_pos, ', dsy:',dsy_pos, ', dsz:',dsz_pos,', dsth:',dsth_pos)
    print('sbx:',sbx_pos, ', sbz:',sbz_pos)
    print('dssx:',dssx_pos,', dssy:',dssy_pos,', dssz:',dssz_pos)

    if moveflag:
        #if det_moveflag:
        #    mov(fdet1.x,fdet1_x)
        #    print('moving flourescence det, wait 5 sec ...')
        #    sleep(5)
        mov(smlld.dsz,dsz_pos)
        mov(smlld.dsx,dsx_pos)
        mov(smlld.dsy,dsy_pos)
        mov(smlld.dsth,dsth_pos)
        mov(smlld.dssx,dssx_pos)
        mov(smlld.dssy,dssy_pos)
        mov(smlld.dssz,dssz_pos)
        if base_moveflag:
            mov(smlld.sbx,sbx_pos)
            mov(smlld.sbz,sbz_pos)

def recover_zp_scan_pos(scan_id,zp_move_flag=0,smar_move_flag=0):
    data = db.get_table(db[scan_id],stream_name='baseline')
    bragg = data.dcm_th[1]
    zpz1 = data.zpz1[1]
    #zpx = data.zpx[1]
    #zpy = data.zpy[1]
    smarx = data.smarx[1]
    smary = data.smary[1]
    smarz = data.smarz[1]
    ssx = data.zpssx[1]
    ssy = data.zpssy[1]
    ssz = data.zpssz[1]
    #print(ssx,ssy,ssz)

    print('scan '+np.str(scan_id))
    print('dcm_th:'+np.str(bragg))
    #print('zpz1: '+np.str(zpz1)+', zpx:'+np.str(zpx)+', zpy:'+np.str(zpy))
    print('zpz1:', np.str(zpz1))
    print('smarx:'+np.str(smarx)+', smary:'+np.str(smary)+', smarz:'+np.str(smarz))

    if zp_move_flag:
        mov(dcm.th,bragg)
        mov(zp.zpz1,zpz1)
        #mov(zp.zpx,zpx)
        #mov(zp.zpy,zpy)

    if smar_move_flag:
        #mov(dcm.th,bragg)
        mov(zps.smarx,smarx)
        mov(zps.smary,smary)
        mov(zps.smarz,smarz)
        mov(zps.zpssx,ssx)
        mov(zps.zpssy,ssy)
        mov(zps.zpssz,ssz)

def stitch_mosaic(start_scan_id, end_scan_id, nx_mosaic, ny_mosaic,elem,norm=None,clim=None,channels=None,cmap='viridis',fill_events=False):
    num_frame = end_scan_id - start_scan_id + 1

    for i in range(num_frame):

        scan_id = start_scan_id + i
        print('loading scan %d' %scan_id)
        if channels is None:
            channels = [1, 2, 3]

        scan_id, df = _load_scan(scan_id, fill_events=fill_events)
        hdr = db[scan_id]['start']
        data = db.get_table(db[scan_id],stream_name='baseline')

        if i == 0:
            nx_flyscan, ny_flyscan = get_flyscan_dimensions(hdr)
            nx = nx_mosaic * nx_flyscan
            ny = ny_mosaic * ny_flyscan
            array = np.zeros((nx,ny))

            x_motor = hdr['motor1']
            x_data = np.asarray(df[x_motor])
            y_motor = hdr['motor2']
            y_data = np.asarray(df[y_motor])
            x_range = np.nanmax(x_data) - np.nanmin(x_data)
            y_range = np.nanmax(y_data) - np.nanmin(y_data)
            smarx = data.smarx[1]
            smary = data.smary[1]
            #print(smarx,smary)
            extent = ((smarx+x_range*nx_mosaic/1000.), (smarx-x_range/2000.),
                      (smary+y_range*ny_mosaic/1000.), (smary-y_range/2000.))



        if elem in df:
            spectrum = np.asarray(df[elem], dtype=np.float32)
        else:
            roi_keys = ['Det%d_%s' % (chan, elem) for chan in channels]

            for key in roi_keys:
                if key not in df:
                    raise KeyError('ROI %s not found' % (key, ))
            spectrum = np.sum([getattr(df, roi) for roi in roi_keys], axis=0)

        if norm is not None:
            monitor = np.asarray(df[norm],dtype=np.float32)
            spectrum = spectrum/(monitor + 1e-8)

        spectrum2 = fly2d_reshape(hdr, spectrum)

        ix = np.mod(i,nx_mosaic)
        iy = np.int(np.floor(i/nx_mosaic))
        array[nx_flyscan*ix:nx_flyscan*(ix+1),ny_flyscan*iy:ny_flyscan*(iy+1)] = spectrum2.T

    if clim is None:
        clim = (np.nanmin(array), np.nanmax(array))

    title = 'Scan id %s. ' % start_scan_id +'- %s' % end_scan_id +' '+ elem
    fig, ax1 = plt.subplots(ncols=1,figsize=(8,5))
    fig.set_tight_layout(True)
    imshow = ax1.imshow(np.fliplr(array.T), extent=extent, interpolation='None', cmap=cmap,
                        vmin=clim[0], vmax=clim[1])
    ax1.set_title('IMSHOW. ' + title)
    fig.gca().invert_xaxis()
    #fig.gca().invert_yaxis()
    ax1.set_xlabel('smarx')
    ax1.set_ylabel('smary')
    fig.colorbar(imshow)

def export_merlin(sid,num=1):
    for i in range(num):
        sid, df = _load_scan(sid, fill_events=False)
        path = os.path.join('/data/users/2017Q2/Robinson_TaS2/', 'scan_{}.txt'.format(sid))
        print('Scan {}. Saving to {}'.format(sid, path))
        #non_objects = [name for name, col in df.iteritems() if col.dtype.name not in ('object', )]
        #dump all data
        non_objects = [name for name, col in df.iteritems()]
        df.to_csv(path, float_format='%1.5e', sep='\t',
                  columns=sorted(non_objects))

        path = os.path.join('/data/users/2017Q2/Robinson_TaS2/', 'scan_{}_scaler.txt'.format(sid))
        #np.savetxt(path, (df['sclr1_ch3'], df['p_ssx'], df['p_ssy']), fmt='%1.5e')
        #np.savetxt(path, (df['sclr1_ch4'], df['zpssx'], df['zpssy']), fmt='%1.5e')
        filename = get_all_filenames(sid,'merlin1')
        num_subscan = len(filename)
        if num_subscan == 1:
            for fn in filename:
                break
            path = os.path.join('/data/users/2017Q2/Robinson_TaS2/', 'scan_{}.h5'.format(sid))
            mycmd = ''.join(['scp', ' ', fn, ' ', path])
            os.system(mycmd)
        else:
            h = db[sid]
            images = db.get_images(h,name='merlin1')
            path = os.path.join('/data/users/2017Q2/Robinson_TaS2/', 'scan_{}.h5'.format(sid))
            f = h5py.File(path, 'w')
            dset = f.create_dataset('/entry/instrument/detector/data', data=images)
            f.close()
            '''''
            j = 1
            for fn in filename:
                path = os.path.join('/home/hyan/export/', 'scan_{}_{}.h5'.format(sid, j))
                mycmd = ''.join(['scp', ' ', fn, ' ', path])
                os.system(mycmd)
                j = j + 1
            '''''
        sid = sid + 1

def position_scan(dsx_list,dsy_list,x_range_list,x_num_list,y_range_list,y_num_list,exp_list):
    x_list = np.array(dsx_list)
    y_list = np.array(dsy_list)
    x_range_list = np.array(x_range_list)
    y_range_list = np.array(y_range_list)
    x_num_list = np.array(x_num_list)
    y_num_list = np.array(y_num_list)
    exp_list = np.array(exp_list)

    dsx_0 = smlld.dsx.position
    dsy_0 = smlld.dsy.position
    num_scan = np.size(x_list)
    #mov(ssa2.hgap,0.03)
    #mov(ssa2.vgap,0.02)
    for i in range(num_scan):
        print('move to position ',i+1,'/',num_scan)
        mov(smlld.dsx,x_list[i])
        mov(smlld.dsy,y_list[i])
        RE(fly2d(smlld.dssx,-x_range_list[i]/2,x_range_list[i]/2,x_num_list[i]/1,smlld.dssy,-y_range_list[i]/2,y_range_list[i]/2,y_num_list[i]/1,exp_list[i]/1))
        #plot2dfly(-1,'Er_L')
        #printfig()
        print('wait 0.2 sec...')
        sleep(0.2)

    mov(smlld.dsx,dsx_0)
    mov(smlld.dsy,dsy_0)
    #mov(ssa2.hgap,0.15)
    #mov(ssa2.vgap,0.05)
