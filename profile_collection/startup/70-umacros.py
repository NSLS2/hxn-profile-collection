import matplotlib.pyplot as plt
import numpy as np
from time import sleep


from scipy.optimize import curve_fit
from databroker import db, get_table
from ophyd import mov, movr


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
        caput('XF:03IDB-PPS{PSh}Cmd:Opn-Cmd', 1)
        sleep(10)
    elif cmd == 'close':
        caput('XF:03IDB-PPS{PSh}Cmd:Cls-Cmd', 1)
        sleep(10)


def z_scan(z_range, num_step):
    z_range = np.float(z_range)
    num_step = np.int(num_step)
    z_step_size = z_range / num_step
    print(z_range, num_step, z_step_size)
    print(-1 * num_step * z_step_size / 2)
    movr(smarz, -1 * num_step * z_step_size / 2)
#    movr(zpsx,0.2588*num_step*z_step_size/2)

    for i in range(num_step + 1):
        fly2d(zpssx, -5, 5, 100, zpssy, -5, 5, 100, 0.025, return_speed=20)
        movr(smarz, z_step_size)
#        movr(zpsx,-0.2588*z_step_size)
    mov(smarz, 1.716)
    mov(zpsx, -0.9288)


def go_det(det):
    if det == 'merlin':
        mov(diff_x, -177)
        mov(diff_y1, -31)
        mov(diff_y2, -31)
    elif det == 'cam11':
        mov(diff_x, 0)
        mov(diff_y1, -0.)
        mov(diff_y2, -0.)
    elif det == 'tpx':
        mov(diff_x, -112)
        mov(diff_y1, -50)
        mov(diff_y2, -50)
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

        angle = angle_start + i * angle_step
        mov(zpsth, angle)
        #x_start_real = x_start / np.cos(angle*np.pi/180.)
        #x_end_real = x_end / np.cos(angle*np.pi/180.)

        dscan(zpsx, -0.01, 0.01, 50, 0.5)
        scan_id, df = _load_scan(-1, fill_events=False)
        x = df['zpsx']
        data = df['Det2_Co']
        x = np.asarray(x)
        data = np.asarray(data)
        #diff = np.diff(data)
        #i_max = np.where(diff == np.max(diff))
        #i_min = np.where(diff == np.min(diff))
        #i_center = np.round((i_max[0][0]+i_min[0][0])/2)+1

        mc = find_mass_center(data)
        mov(zpsx, x[mc] - 0.001)

        if np.abs(angle) <= 45:
            x_start_real = x_start / np.cos(angle * np.pi / 180.)
            x_end_real = x_end / np.cos(angle * np.pi / 180.)
            fly2d(zpssx, x_start_real, x_end_real, x_num, zpssy,
                  y_start, y_end, y_num, exposure, return_speed=40)

        else:
            x_start_real = x_start / np.abs(np.sin(angle * np.pi / 180.))
            x_end_real = x_end / np.abs(np.sin(angle * np.pi / 180.))
            fly2d(zpssz, x_start_real, x_end_real, x_num, zpssy,
                  y_start, y_end, y_num, exposure, return_speed=40)

        print('waiting for 120 sec...')
        sleep(120)

    mov(zpsth, 0)


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


def _load_scan(scan_id, fill_events=False):
    '''Load scan from databroker by scan id'''
    if scan_id > 0 and scan_id in data_cache:
        df = data_cache[scan_id]
    else:
        hdr = db[scan_id]
        scan_id = hdr['start'].scan_id
        if scan_id in data_cache:
            df = data_cache[scan_id]
        else:
            data_cache[scan_id] = get_table(hdr, fill=fill_events)

    return scan_id, df


def movr_zpz1(dz):
    movr(zpz1, dz)
    movr(zpx, dz * 4.428)
    movr(zpy, -dz * 1.490)


def reset_tpx(num):
    for i in range(1000):
        timepix1.cam.num_images.put(num, wait=False)
        sleep(0.5)


def mov_gamma(gamma):

    diff_z_pos = diff.z.position
    diff_yaw_pos = diff.yaw.position * np.pi / 180.0
    diff_cz_pos = diff.cz.position
    diff_x_pos = diff.x.position

    r_yaw = 581.20

    # when both diff_z and diff_cz are at positive limit:
    r0 = 574.668 + diff_z_pos
    beta = 89.377 * np.pi / 180.0
    gamma = gamma * np.pi / 180.0

    current_r = (sin(beta) / sin(diff_yaw_pos + beta)) * \
        (r0 + r_yaw * (1.0 - cos(diff_yaw_pos))) + diff_cz_pos
    print('current r =', current_r)

    x = r_yaw * sin(gamma) + (sin(gamma) / sin(gamma + beta)) * \
        (r0 + r_yaw * (1.0 - cos(gamma)))
    r = (sin(beta) / sin(gamma + beta)) * (r0 + r_yaw * (1.0 - cos(gamma)))

    print('gamma =', gamma * 180.0 / np.pi, '   x =', x, '   r =', r)
    print('diff_x=', -x, 'diff_cz=', current_r - r)
    sleep(10)
    mov(diff.yaw, gamma * 180. / np.pi, wait=False)
    sleep(1)
    mov(diff.x, -x, wait=False)
    sleep(1)
    mov(diff.cz, current_r - r, wait=False)

#    mov(diff_yaw, gamma*180.0/pi)
#    mov(diff_x, -x)
#    mov(diff_cz, current_r - r)


def mov_delta(delta):
    diff_z_pos = diff.z.position
    diff_yaw_pos = diff.yaw.position * np.pi / 180.0
    diff_cz_pos = diff.cz.position
    diff_x_pos = diff.x.position

    r0 = 574.668 + diff_z_pos
    z1 = r0 + 395.2
    z2 = z1 + 380

    z1 = z1 / cos(diff_yaw_pos)
    z2 = z1 + 380

    delta = delta * np.pi / 180.0
    y1 = tan(delta) * z1
    y2 = tan(delta) * z2

    dr = z1 / cos(delta) - z1

    print('y1 =', y1, '   y2 =', y2, '  dr =', dr)

    sleep(10)

    diff.y1.move(y1, wait=False)
    sleep(0.5)
    diff.y2.move(y2, wait=False)


def th_fly1d(th_start, th_end, num, m_start, m_end, m_num, sec):
    th_step = (th_end - th_start) / num
    movr(zpsth, th_start)
    for i in range(num + 1):
        fly1d(zpssy, m_start, m_end, m_num, sec)
        movr(zpsth, th_step)
    movr(zpsth, -(th_end + th_step))


def th_fly2d(th_start, th_end, num, x_start, x_end, x_num, y_start, y_end,
             y_num, sec):
    th_step = (th_end - th_start) / num
    movr(zpsth, th_start)
    for i in range(num + 1):
        fly2d(zpssx, x_start, x_end, x_num, zpssy, y_start, y_end, y_num, sec)
        movr(zpsth, th_step)
    movr(zpsth, -(th_end + th_step))


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
        mov(diff.x, -x_yaw, wait=False)
        sleep(1)
        mov(diff.yaw, gamma * 180. / np.pi, wait=False)
        sleep(0.5)
        mov(diff.cz, dz, wait=False)


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
