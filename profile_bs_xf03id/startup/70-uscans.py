def zscan(z_start, z_end, z_num):
    movr(hz, z_start)
    movr(vz, z_start)
    fly2d(ssx, -.25, 0.25, 50, ssy, -0.25, 0.25, 50, 0.1)
    plot2dfly(-1, elem='Ni', cmap='jet', shift_zeros=False)
    step = (z_end - z_start)/z_num
    for i in range(z_num):
        movr(hz, step)
        movr(vz, step)
        fly2d(ssx, -.25, 0.25, 50, ssy, -0.25, 0.25, 50, 0.1)
        plot2dfly(-1, elem='Ni', cmap='jet', shift_zeros=False)
    movr(hz, -z_end)
    movr(vz, -z_end)

def userscan():
    fly2d(ssx, -3, 3, 200, ssy, -3, 3, 200, 0.12)
    fly2d(ssx, -3, 3, 100, ssy, -3, 3, 100, 0.5)
    mesh ([ssy, -3.2, 3.2], [ssx, -3.2, 3.2], [80, 80], 1)       
    fly2d(ssx, -3.5, 3.5, 140, ssy, -3.5, 3.5, 140, 0.2)

