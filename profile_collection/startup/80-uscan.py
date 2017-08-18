'''
from epics import caput,caget

def fly1d_user(motor,start,end,num_pos,exp):
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"R1")
    RE(fly1d(motor,start,end,num_pos,exp))
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"I1")


def fly2d_user(motor1,start1,end1,num1,motor2,start2,end2,num2,exp):
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"R1")
    RE(fly2d(motor1,start1,end1,num1,motor2,start2,end2,num2,exp))
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"I1")
'''
def theta_fly2d(angle_start,angle_end,num_angle,motor1,start1,end1,num1,motor2,start2,end2,num2,exp):
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
    #RE(fly1d(zps.zpssz,-1,1,50,0.1))
    #mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
    
    for i in range(num_angle+1):
        RE(fly2d(motor1,start1,end1,num1,motor2,start2,end2,num2,exp))
        mov_to_image_cen_zpss(-1,'Ga',1)
        export(-1)
        movr(zps.zpsth,angle_step)
        #RE(fly1d(zps.zpssz,-.25,.25,20,0.1))
        #mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
        print('wait 2 sec ...')
        sleep(2)
        
    mov(zps.zpsth,angle_current)

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

    angle_step = (angle_end - angle_start) / num_angle
    start_angle = angle_current + angle_start
    mov(zps.zpsth,start_angle)
   # RE(fly1d(zps.zpssz,-1,1,50,0.1))
   # mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
    
    for i in range(num_angle+1):
        RE(fly2d(motor1,start1,end1,num1,motor2,start2,end2,num2,exp))
    #    mov_to_image_cen_zpss(-1,'Ga',1)
        export(-1)
        movr(zps.zpsth,angle_step)
        movr(zps.zpssz,-2*angle_step)
    #    RE(fly1d(zps.zpssz,-.25,.25,20,0.1))
    #    mov_to_line_center(-1,elem='Ga',moveflag=1,threshold=0.2,movepiezoflag=1)
        print('wait 2 sec ...')
        sleep(2)
        
   # mov(zps.zpsth,angle_current)

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
