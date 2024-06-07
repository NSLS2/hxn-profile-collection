def plotlastfluo(id=-1,elem = 'Ni'):
    st = db[id].start['scan']
    scan_size = [st['scan_input'][2],st['scan_input'][5]]
    fluo = np.zeros(scan_size[1]*scan_size[0])
    for roi in xspress3.enabled_rois:
        if elem in roi.name:
            f = roi.settings.array_data.get()
            fluo[:len(f)] += f
    f2d = fluo.reshape(scan_size[1],scan_size[0])

    plt.close(110)
    plt.figure(110)
    plt.imshow(f2d)
    return f2d

def cenlastfluo():
    st = db[-1].start['scan']
    if st['fast_axis']['motor_name'] == 'zpssx':
        px = db.reg.retrieve(db[-1].table()['inenc3_val'][1])
        px = (px - np.mean(px))*9.7e-5
    else:
        px = db.reg.retrieve(db[-1].table()['inenc4_val'][1])
        px = -(px - np.mean(px))*1.006e-4
    scan_size = [st['scan_input'][2],st['scan_input'][5]]
    py = np.linspace(st['scan_input'][3],st['scan_input'][4],scan_size[0]*scan_size[1])
    fluo = np.zeros(scan_size[1]*scan_size[0])
    for roi in xspress3.enabled_rois:
        if 'Ni' in roi.name:
            fluo += roi.settings.array_data.get()[:len(fluo)]
    return get_masscenter(px,py,fluo)


def zp_tomo_scan_rapid(angle_start, angle_end, angle_step, x_start, x_end, x_num,
              y_start, y_end, y_num, exposure, elem,save_file,ic_0=None):
    #if os.path.isfile('rotCali'):
    #    caliFile = open('rotCali','rb')
    #    y = pickle.load(caliFile)
    angle_start = float(angle_start)
    angle_end = float(angle_end)
    angle_step = float(angle_step)
    x_start = float(x_start)
    x_end = float(x_end)
    x_num = int(x_num)
    y_start = float(y_start)
    y_end = float(y_end)
    y_num = int(y_num)
    exposure = float(exposure)
    #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0)
    angle_num = int(np.ceil(np.abs((angle_end-angle_start)/angle_step)))

    x_scale_factor = 0.9542
    z_scale_factor = 1.0309

    if ic_0 is None:
        ic_0 = sclr2_ch2.get()
    for i in range(angle_num + 1):
        yield from bps.mov(zpssy,0)

        angle = angle_start + i * angle_step * np.sign(angle_end-angle_start)
        yield from bps.mov(zps.zpsth, angle)

        #yield from bps.mov(zpssx,0)
        #yield from bps.mov(zpssy,0)
        #yield from bps.mov(zpssz,0)

        while (sclr2_ch2.get() < 10000):
            yield from bps.sleep(60)
            print('IC1 is lower than 1000, waiting...')
        #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',1)
        #yield from bps.sleep(3)
        #caput('XF:03IDC-ES{Zeb:2}:SOFT_IN:B0',0)
        #yield from bps.sleep(3)

        angle_offset = -0.2
        RE.md['tomo_angle_offset'] = angle_offset
        RE.md['x_scale_factor'] = x_scale_factor
        RE.md['z_scale_factor'] = z_scale_factor

        if np.abs(angle-angle_offset) <= 45.:
            #yield from bps.movr(zpssy,-1)
            yield from fly2dcontpd([fs,eiger2,xspress3],zpssx,-8,8,20,zpssy,-8,8,20,0.01,dead_time = 0.001)
            yield from bps.sleep(0.5)
            plotlastfluo();
            plt.pause(0.05)
            plt.show()
            cen = cenlastfluo()
            if not np.isnan(cen[0]):
                yield from bps.movr(smarx,cen[0]/1000)
            if not np.isnan(cen[1]):
                yield from bps.movr(smary,(cen[1]+2.0)/1000)
        else:
            yield from fly2dcontpd([fs,eiger2,xspress3],zpssz,-8,8,20,zpssy,-8,8,20,0.01,dead_time = 0.001)
            yield from bps.sleep(0.5)
            plotlastfluo();
            cen = cenlastfluo()
            if not np.isnan(cen[0]):
                yield from bps.movr(smarz,cen[0]/1000)
            if not np.isnan(cen[1]):
                yield from bps.movr(smary,(cen[1]+2.0)/1000)
        #eiger2.hdf5.warmup()


        if np.abs(angle-angle_offset) <= 45.0:
            # yield from fly2d(dets1,zpssx,-6.5,7,18,zpssy,-5,5.5,14,0.05,return_speed=40)
            # yield from mov_to_image_cen_dsx(-1)

            x_start_real = x_start / np.cos((angle-angle_offset) * np.pi / 180.)/x_scale_factor
            x_end_real = x_end / np.cos((angle-angle_offset) * np.pi / 180.)/x_scale_factor
            y_start_real = y_start
            y_end_real = y_end
            #yield from fly2d([fs, zebra, sclr1, xspress3], zpssy, y_start, y_end, y_num,
            #                 zpssx, x_start_real, x_end_real, x_num, exposure, return_speed=40)
            #RE(fly2d(zpssx, x_start_real, x_end_real, x_num, zpssy,
            #         y_start, y_end, y_num, exposure, return_speed=40))
            yield from fly2dcontpd([fs,eiger2],zpssx,x_start_real, x_end_real, x_num,zpssy,y_start_real,y_end_real,y_num,exposure)

        else:
            # yield from fly2d(dets1,zpssz,-6.5,7,18,zpssy,-5,5.5,14,0.05,return_speed=40)
            # yield from mov_to_image_cen_dsx(-1)

            x_start_real = x_start / np.abs(np.sin((angle-angle_offset) * np.pi / 180.))/z_scale_factor
            x_end_real = x_end / np.abs(np.sin((angle-angle_offset) * np.pi / 180.))/z_scale_factor
            y_start_real = y_start
            y_end_real = y_end
            #yield from fly2d([fs, zebra, sclr1, xspress3],zpssy, y_start, y_end, y_num,
            #                 zpssz, x_start_real, x_end_real, x_num, exposure, return_speed=40)
            #RE(fly2d(zpssz, x_start_real, x_end_real, x_num, zpssy,
            #         y_start, y_end, y_num, exposure, return_speed=40))
            yield from fly2dcontpd([fs,eiger2],zpssz,x_start_real, x_end_real, x_num,zpssy,y_start_real,y_end_real,y_num,exposure)

        #mov_to_image_cen_smar(-1)
        #yield from mov_to_image_cen_dsx(-1)
        #plot2dfly(-1,elem,'sclr1_ch4')
        #insertFig(note='zpsth = {}'.format(check_baseline(-1,'zpsth')))
        #plt.close()
        #merlin2.unstage()
        #xspress3.unstage()
        #yield from bps.sleep(5)
        #insert_xrf_map_to_pdf(-1, elem, title_=['zpsth'])
        flog = open(save_file,'a')
        flog.write('%d %.2f\n'%(db[-1].start['scan_id'],zpsth.position))
        flog.close()
        #if (sclr2_ch2.get() < (0.85*ic_0)):
        #    yield from peak_the_flux()
        #if np.remainder(i+1,5)==0:
        #    yield from peak_bpm_x(-20, 20, 10)
        #    yield from peak_bpm_y(-10, 10, 10)
    save_page()

