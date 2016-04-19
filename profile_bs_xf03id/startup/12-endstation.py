from ophyd import (EpicsMotor, Device, Component as Cpt,
                   EpicsSignalRO, PVPositioner)

class HxnSSAperture(Device):
    hgap = Cpt(EpicsMotor, '-Ax:XAp}Mtr')
    vgap = Cpt(EpicsMotor, '-Ax:YAp}Mtr')
    hcen = Cpt(EpicsMotor, '-Ax:X}Mtr')
    vcen = Cpt(EpicsMotor, '-Ax:Y}Mtr')


ssa1 = HxnSSAperture('XF:03IDB-OP{Slt:SSA1', name='ssa1')
ssa2 = HxnSSAperture('XF:03IDC-OP{Slt:SSA2', name='ssa2')

bpm6_y = EpicsMotor('XF:03IDB-OP{BPM:6-Ax:Y}Mtr', name='bpm6_y')

# idb_m1 = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:6}Mtr', name='idb_m1')
# idb_m2 = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:7}Mtr', name='idb_m2')
# idb_m3 = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:8}Mtr', name='idb_m3')

s3 = HxnSlitA('XF:03IDC-OP{Slt:3', name='s3')


class HxnTurboPmacController(Device):
    m1 = Cpt(EpicsMotor, '-Ax:1}Mtr')
    m2 = Cpt(EpicsMotor, '-Ax:2}Mtr')
    m3 = Cpt(EpicsMotor, '-Ax:3}Mtr')
    m4 = Cpt(EpicsMotor, '-Ax:4}Mtr')
    m5 = Cpt(EpicsMotor, '-Ax:5}Mtr')
    m6 = Cpt(EpicsMotor, '-Ax:6}Mtr')
    m7 = Cpt(EpicsMotor, '-Ax:7}Mtr')
    m8 = Cpt(EpicsMotor, '-Ax:8}Mtr')


# Unpopulated motor controllers:
mc2 = HxnTurboPmacController('XF:03IDC-ES{MC:2', name='mc2')
mc3 = HxnTurboPmacController('XF:03IDC-ES{MC:3', name='mc3')
mc4 = HxnTurboPmacController('XF:03IDC-ES{MC:4', name='mc4')


class HxnSlitB(Device):
    '''HXN slit device, with X/Y/Z/top'''
    vgap = Cpt(EpicsMotor, '-Ax:X}Mtr')
    vcen = Cpt(EpicsMotor, '-Ax:Y}Mtr')
    hgap = Cpt(EpicsMotor, '-Ax:Z}Mtr')
    hcen = Cpt(EpicsMotor, '-Ax:Top}Mtr')


s4 = HxnSlitB('XF:03IDC-ES{Slt:4', name='s4')

# mc6_m5 = EpicsMotor('XF:03IDC-ES{MC:6-Ax:5}Mtr', name='mc6_m5')
# mc6_m6 = EpicsMotor('XF:03IDC-ES{MC:6-Ax:6}Mtr', name='mc6_m6')
# mc6_m7 = EpicsMotor('XF:03IDC-ES{MC:6-Ax:7}Mtr', name='mc6_m7')

bpm7_y = EpicsMotor('XF:03IDC-ES{BPM:7-Ax:Y}Mtr', name='bpm7_y')

mc7 = HxnTurboPmacController('XF:03IDC-ES{MC:7', name='mc7')

questar_f = EpicsMotor('XF:03IDC-ES{MC:8-Ax:1}Mtr', name='questar_f')

mc8 = HxnTurboPmacController('XF:03IDC-ES{MC:8', name='mc8')
# mc9 = HxnTurboPmacController('XF:03IDC-ES{MC:9', name='mc9')


class HxnSlitC(Device):
    '''HXN slit device, with vertical/horizontal gaps/centers'''
    vgap = Cpt(EpicsMotor, '-Ax:Vgap}Mtr')
    vcen = Cpt(EpicsMotor, '-Ax:Vcen}Mtr')
    hgap = Cpt(EpicsMotor, '-Ax:Hgap}Mtr')
    vcen = Cpt(EpicsMotor, '-Ax:Hcen}Mtr')


s5 = HxnSlitC('XF:03IDC-ES{Slt:5', name='s5')
s6 = HxnSlitC('XF:03IDC-ES{Slt:6', name='s6')


# mc10 = HxnTurboPmacController('XF:03IDC-ES{MC:10', name='mc10')


class HxnDetectorPositioner(Device):
    '''HXN X/Y positioner device'''
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr')
    z = Cpt(EpicsMotor, '-Ax:Z}Mtr')


fdet1 = HxnDetectorPositioner('XF:03IDC-ES{Det:Vort', name='fdet1')
fdet2 = HxnDetectorPositioner('XF:03IDC-ES{Det:Bruk', name='fdet2')

bs_x = EpicsMotor('XF:03IDC-ES{MC:12-Ax:4}Mtr', name='bs_x')
bs_y = EpicsMotor('XF:03IDC-ES{MC:12-Ax:5}Mtr', name='bs_y')

mc12 = HxnTurboPmacController('XF:03IDC-ES{MC:12', name='mc12')


class DetectorStation(Device):
    z = Cpt(EpicsMotor, '-Ax:Z}Mtr')
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    y1 = Cpt(EpicsMotor, '-Ax:Y1}Mtr')
    y2 = Cpt(EpicsMotor, '-Ax:Y2}Mtr')
    yaw = Cpt(EpicsMotor, '-Ax:Yaw}Mtr')
    c1 = Cpt(EpicsMotor, '-Ax:C1}Mtr')
    c2 = Cpt(EpicsMotor, '-Ax:C2}Mtr')
    c3 = Cpt(EpicsMotor, '-Ax:C3}Mtr')


diff = DetectorStation('XF:03IDC-ES{Diff', name='diff')
s7 = HxnSlitC('XF:03IDC-ES{Slt:7', name='s7')
