from epics import caput,caget

def fly1d_user(motor,start,end,num_pos,exp):
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"R1")
    RE(fly1d(motor,start,end,num_pos,exp))
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"I1")


def fly2d_user(motor1,start1,end1,num1,motor2,start2,end2,num2,exp):
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"R1")
    RE(fly2d(motor1,start1,end1,num1,motor2,start2,end2,num2,exp))
    caput('XF:03IDC-ES{Flt:2}sendCommand.VAL',"I1")

