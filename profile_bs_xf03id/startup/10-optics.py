from ophyd import (EpicsMotor, Device, Component as Cpt,
                   EpicsSignalRO)
from ophyd.device import FormattedComponent as FCpt


class HxnDCM(Device):
    '''HXN DCM Device'''
    th = Cpt(EpicsMotor, 'XF:03IDA-OP{Mon:1-Ax:Bragg}Mtr')
    x = Cpt(EpicsMotor, 'XF:03IDA-OP{Mon:1-Ax:X}Mtr')
    p = Cpt(EpicsMotor, 'XF:03IDA-OP{Mon:1-Ax:P}Mtr')
    r = Cpt(EpicsMotor, 'XF:03IDA-OP{Mon:1-Ax:R}Mtr')
    pf = Cpt(EpicsMotor, 'XF:03IDA-OP{Mon:1-Ax:PF}Mtr')
    rf = Cpt(EpicsMotor, 'XF:03IDA-OP{Mon:1-Ax:RF}Mtr')


dcm = HxnDCM('', name='dcm')


class HxnMirror1(Device):
    '''HXN Mirror 1 device'''
    x = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:1-Ax:X}Mtr')
    y = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:1-Ax:Y}Mtr')
    p = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:1-Ax:P}Mtr')
    b = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:1-Ax:Bend}Mtr')
    pf = Cpt(EpicsMotor, 'XF:03IDA-OP{HCM:1-Ax:PF}Mtr')


class HxnMirror2(Device):
    '''HXN Mirror 2 device'''
    x = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:2-Ax:X}Mtr')
    y = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:2-Ax:Y}Mtr')
    p = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:2-Ax:P}Mtr')
    b = Cpt(EpicsMotor, 'XF:03IDA-OP{Mir:2-Ax:Bend}Mtr')
    pf = Cpt(EpicsMotor, 'XF:03IDA-OP{HFM:1-Ax:PF}Mtr')


m1 = HxnMirror1('', name='m1')
m2 = HxnMirror2('', name='m2')


class HxnSlitA(Device):
    '''HXN slit device, with top/bottom/inboard/outboard'''
    bot = Cpt(EpicsMotor, '-Ax:Btm}Mtr')
    top = Cpt(EpicsMotor, '-Ax:Top}Mtr')
    inb = Cpt(EpicsMotor, '-Ax:Inb}Mtr')
    outb = Cpt(EpicsMotor, '-Ax:Outb}Mtr')


class HxnSlitA1(HxnSlitA):
    #           ^^^^^^^^ means it includes 'bot, top, inb, outb' too
    # x, y position from the i400 IOC:
    xpos = FCpt(EpicsSignalRO, 'XF:03IDA-BI{{Slt:1}}PosX-I')
    ypos = FCpt(EpicsSignalRO, 'XF:03IDA-BI{{Slt:1}}PosY-I')


s1 = HxnSlitA1('XF:03IDA-OP{Slt:1', name='s1')
s2 = HxnSlitA('XF:03IDA-OP{Slt:2', name='s2')


class HxnI400(Device):
    '''HXN I400 BPM current readout'''
    # raw currents
    i_top = Cpt(EpicsSignalRO, 'I:Raw1-I')
    i_bottom = Cpt(EpicsSignalRO, 'I:Raw2-I')
    i_right = Cpt(EpicsSignalRO, 'I:Raw3-I')
    i_left = Cpt(EpicsSignalRO, 'I:Raw4-I')

    # x/y position
    x = Cpt(EpicsSignalRO, 'PosX-I')
    y = Cpt(EpicsSignalRO, 'PosY-I')


# Slit 1 BPM (drain current from I400)
s1_bpm = HxnI400('XF:03IDA-BI{Slt:1}', name='s1_bpm')


class HxnXYPositioner(Device):
    '''HXN X/Y positioner device'''
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr')


class HxnXYPitchPositioner(Device):
    '''HXN X/Y/Pitch positioner'''
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr')
    p = Cpt(EpicsMotor, '-Ax:P}Mtr')


bpm1 = HxnXYPositioner('XF:03IDA-OP{BPM:1', name='bpm1')
bpm2 = HxnXYPositioner('XF:03IDA-OP{BPM:2', name='bpm2')

bpm3_x = EpicsMotor('XF:03IDA-OP{BPM:3-Ax:X}Mtr', name='bpm3_x')
fs1_y = EpicsMotor('XF:03IDA-OP{FS:1-Ax:Y}Mtr', name='fs1_y')
bpm4_y = EpicsMotor('XF:03IDA-OP{BPM:4-Ax:Y}Mtr', name='bpm4_y')
bpm5_y = EpicsMotor('XF:03IDA-OP{BPM:5-Ax:Y}Mtr', name='bpm5_y')

fl1_y = EpicsMotor('XF:03IDA-OP{Flr:1-Ax:Y}Mtr', name='fl1_y')
fl2_y = EpicsMotor('XF:03IDA-OP{Flr:2-Ax:Y}Mtr', name='fl2_y')

crl = HxnXYPitchPositioner('XF:03IDA-OP{Lens:CRL', name='crl')


# # HCM
# m1pz = EpicsMotor('XF:03IDA-OP{HCM:1-Ax:PF}Mtr', name='m1pz')
#
# # Mon:1
# bragg = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:Bragg}Mtr', name='bragg')
#
# # Diagnostic Manipulators
# f1y = EpicsMotor('XF:03IDA-OP{Flr:1-Ax:Y}Mtr', name='f1y')
# f2y = EpicsMotor('XF:03IDA-OP{Flr:2-Ax:Y}Mtr', name='f2y')
#
# # nanoBPM2@SSA1
# nano2y = EpicsMotor('XF:03IDB-OP{BPM:6-Ax:Y}Mtr', name='nano2y')
#
# # nanoBPM3@SSA2
# nano3y = EpicsMotor('XF:03IDC-OP{BPM:7-Ax:Y}Mtr', name='nano3y')
