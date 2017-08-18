def night_scan():
    print('this is a night scan')
    ang_list = np.arange(73.8, 72.8,-0.04)
    print(ang_list)
    num_ang = np.size(ang_list)
    for ii in range(num_ang):
        mov(zps.zpsth,ang_list[ii])
        #RE(fly1d(zpssx, -.75, .75, 100, 0.1))
        RE(fly1d(zpssx, -2.5, 2.5, 50, 0.1))
        mov_to_line_center(-1,'Ga',threshold=5)
        #RE(fly1d(zpssx, -.3, .3, 30, 0.2))
        #mov_to_line_center(-1,'Ga',threshold=5,moveflag=1,movepiezoflag=1)
        print('wait for 2 sec')
        sleep(2)
        #RE(fly2d(zpssy, -0.9, 2.1, 30, zpssx, -0.9, 1.1, 20, 0.05))
        RE(fly2d(zpssy, -1.1, 1.8, 145, zpssx, -1, 1, 29, 0.2))

        export(-1)

def scan2d_user():
   # x_list = np.arange(-0.65,0.7,.07)
   # print(x_list)
    num_x = 19 #np.size(x_list)
    for ii in range(num_x):
        movr(smarx,0.00007)
        RE(dscan(smary,-.00025,.00025,25,5))
        print('wait for 2 sec')
        sleep(2)
        export(-1)
        plot(-1,'Ga')


def export_experiment():
    exps = np.arange(24369, 26110, 1)
    export(exps)
    #nexps = np.size(exps)
    #for ii in range(nexps):
    #    export(exps[ii])
    #    print('pause 10')
    #    sleep(10)


def night_tomo_scan(angle1, angle2,xs1, xe1, xs2, xe2, x_num, y_start, y_end, y_num, exposure):

    tomo_scan_list_zp(angle1, xs1, xe1, x_num, y_start, y_end, y_num, exposure)
    sleep(10)


    tomo_scan_list_zp(angle2, xs2, xe2, x_num, y_start, y_end, y_num, exposure)
    sleep(10)

def fermat_mosaic():
    bx_start = p_bx.position
    cy_start = p_cy.position
    movr(p_cy,-10)
    for j in range(3):
        mov(p_bx,bx_start)
        movr(p_bx,-0.03)
        for i in range(6):
            RE(fermat(p_ssx,p_ssy, 10, 10, 0.15, 1, 0.5))
            movr(p_bx,0.01)
        movr(p_cy,10)

    mov(p_bx,bx_start)
    mov(p_cy,cy_start)
