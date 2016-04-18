from ophyd import (EpicsMotor, Device, Component as Cpt,
                   EpicsSignalRO, PVPositioner)

class HxnSSAperture(Device):
    hgap = EpicsMotor('-Ax:XAp}Mtr')
    vgap = EpicsMotor('-Ax:YAp}Mtr')
    hcen = EpicsMotor('-Ax:X}Mtr')
    vcen = EpicsMotor('-Ax:Y}Mtr')


ssa1 = HxnSSAperture('XF:03IDB-OP{Slt:SSA1', name='ssa1')
ssa2 = HxnSSAperture('XF:03IDB-OP{Slt:SSA2', name='ssa2')

bpm6_y = EpicsMotor('XF:03IDB-OP{BPM:6-Ax:Y}Mtr', name='bpm6_y')

# idb_m1 = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:6}Mtr', name='idb_m1')
# idb_m2 = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:7}Mtr', name='idb_m2')
# idb_m3 = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:8}Mtr', name='idb_m3')

s3 = HxnSlitA('XF:03IDC-OP{Slt:3', name='s3')

# mc2_m1 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:1}Mtr', name='mc2_m1')
# mc2_m2 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:2}Mtr', name='mc2_m2')
# mc2_m3 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:3}Mtr', name='mc2_m3')
# mc2_m4 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:4}Mtr', name='mc2_m4')
# mc2_m5 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:5}Mtr', name='mc2_m5')
# mc2_m6 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:6}Mtr', name='mc2_m6')
# mc2_m7 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:7}Mtr', name='mc2_m7')
# mc2_m8 = EpicsMotor('XF:03IDC-ES{MC:2-Ax:8}Mtr', name='mc2_m8')
# mc3_m1 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:1}Mtr', name='mc3_m1')
# mc3_m2 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:2}Mtr', name='mc3_m2')
# mc3_m3 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:3}Mtr', name='mc3_m3')
# mc3_m4 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:4}Mtr', name='mc3_m4')
# mc3_m5 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:5}Mtr', name='mc3_m5')
# mc3_m6 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:6}Mtr', name='mc3_m6')
# mc3_m7 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:7}Mtr', name='mc3_m7')
# mc3_m8 = EpicsMotor('XF:03IDC-ES{MC:3-Ax:8}Mtr', name='mc3_m8')
# mc4_m1 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:1}Mtr', name='mc4_m1')
# mc4_m2 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:2}Mtr', name='mc4_m2')
# mc4_m3 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:3}Mtr', name='mc4_m3')
# mc4_m4 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:4}Mtr', name='mc4_m4')
# mc4_m5 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:5}Mtr', name='mc4_m5')
# mc4_m6 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:6}Mtr', name='mc4_m6')
# mc4_m7 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:7}Mtr', name='mc4_m7')
# mc4_m8 = EpicsMotor('XF:03IDC-ES{MC:4-Ax:8}Mtr', name='mc4_m8')

p_vx = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:1}Mtr', name='p_vx')
p_vy = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:2}Mtr', name='p_vy')
p_vt = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:3}Mtr', name='p_vt')
p_vchi = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:4}Mtr', name='p_vchi')
p_osat = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:5}Mtr', name='p_osat')
p_osay = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:6}Mtr', name='p_osay')
p_osax = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:7}Mtr', name='p_osax')
p_osaz = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:8}Mtr', name='p_osaz')


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

# mc7_m1 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:1}Mtr', name='mc7_m1')
# mc7_m2 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:2}Mtr', name='mc7_m2')
# mc7_m3 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:3}Mtr', name='mc7_m3')
# mc7_m4 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:4}Mtr', name='mc7_m4')
# mc7_m5 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:5}Mtr', name='mc7_m5')
# mc7_m6 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:6}Mtr', name='mc7_m6')
# mc7_m7 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:7}Mtr', name='mc7_m7')
# mc7_m8 = EpicsMotor('XF:03IDC-ES{MC:7-Ax:8}Mtr', name='mc7_m8')

questar_f = EpicsMotor('XF:03IDC-ES{MC:8-Ax:1}Mtr', name='questar_f')

# mc8_m2 = EpicsMotor('XF:03IDC-ES{MC:8-Ax:2}Mtr', name='mc8_m2')
# mc8_m3 = EpicsMotor('XF:03IDC-ES{MC:8-Ax:3}Mtr', name='mc8_m3')
# mc8_m4 = EpicsMotor('XF:03IDC-ES{MC:8-Ax:4}Mtr', name='mc8_m4')
# mc8_m5 = EpicsMotor('XF:03IDC-ES{MC:8-Ax:5}Mtr', name='mc8_m5')
# mc8_m6 = EpicsMotor('XF:03IDC-ES{MC:8-Ax:6}Mtr', name='mc8_m6')
# mc8_m7 = EpicsMotor('XF:03IDC-ES{MC:8-Ax:7}Mtr', name='mc8_m7')
# mc8_m8 = EpicsMotor('XF:03IDC-ES{MC:8-Ax:8}Mtr', name='mc8_m8')
# mc9_m1 = EpicsMotor('XF:03IDC-ES{MC:9-Ax:1}Mtr', name='mc9_m1')
# mc9_m2 = EpicsMotor('XF:03IDC-ES{MC:9-Ax:2}Mtr', name='mc9_m2')
# mc9_m3 = EpicsMotor('XF:03IDC-ES{MC:9-Ax:3}Mtr', name='mc9_m3')
# mc9_m4 = EpicsMotor('XF:03IDC-ES{MC:9-Ax:4}Mtr', name='mc9_m4')

class HxnSlitC(Device):
    '''HXN slit device, with vertical/horizontal gaps/centers'''
    vgap = Cpt(EpicsMotor, '-Ax:Vgap}Mtr')
    vcen = Cpt(EpicsMotor, '-Ax:Vcen}Mtr')
    hgap = Cpt(EpicsMotor, '-Ax:Hgap}Mtr')
    vcen = Cpt(EpicsMotor, '-Ax:Hcen}Mtr')


s5 = HxnSlitC('XF:03IDC-ES{Slt:5', name='s5')
s6 = HxnSlitC('XF:03IDC-ES{Slt:6', name='s6')


# mc10_m5 = EpicsMotor('XF:03IDC-ES{MC:10-Ax:5}Mtr', name='mc10_m5')
# mc10_m6 = EpicsMotor('XF:03IDC-ES{MC:10-Ax:6}Mtr', name='mc10_m6')
# mc10_m7 = EpicsMotor('XF:03IDC-ES{MC:10-Ax:7}Mtr', name='mc10_m7')
# mc10_m8 = EpicsMotor('XF:03IDC-ES{MC:10-Ax:8}Mtr', name='mc10_m8')

class HxnDetectorPositioner(Device):
    '''HXN X/Y positioner device'''
    x = Cpt(EpicsMotor, '-Ax:X}Mtr')
    y = Cpt(EpicsMotor, '-Ax:Y}Mtr')
    z = Cpt(EpicsMotor, '-Ax:Z}Mtr')


fdet1 = HxnDetectorPositioner('XF:03IDC-ES{Det:Vort', name='fdet1')
fdet2 = HxnDetectorPositioner('XF:03IDC-ES{Det:Bruk', name='fdet2')

# mc11_m7 = EpicsMotor('XF:03IDC-ES{MC:11-Ax:7}Mtr', name='mc11_m7')
# mc11_m8 = EpicsMotor('XF:03IDC-ES{MC:11-Ax:8}Mtr', name='mc11_m8')
# mc12_m1 = EpicsMotor('XF:03IDC-ES{MC:12-Ax:1}Mtr', name='mc12_m1')
# mc12_m2 = EpicsMotor('XF:03IDC-ES{MC:12-Ax:2}Mtr', name='mc12_m2')
# mc12_m3 = EpicsMotor('XF:03IDC-ES{MC:12-Ax:3}Mtr', name='mc12_m3')

bs_x = EpicsMotor('XF:03IDC-ES{MC:12-Ax:4}Mtr', name='bs_x')
bs_y = EpicsMotor('XF:03IDC-ES{MC:12-Ax:5}Mtr', name='bs_y')
# mc12_m6 = EpicsMotor('XF:03IDC-ES{MC:12-Ax:6}Mtr', name='mc12_m6')
# mc12_m7 = EpicsMotor('XF:03IDC-ES{MC:12-Ax:7}Mtr', name='mc12_m7')
# mc12_m8 = EpicsMotor('XF:03IDC-ES{MC:12-Ax:8}Mtr', name='mc12_m8')

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

# mc14_m5 = EpicsMotor('XF:03IDC-ES{MC:14-Ax:5}Mtr', name='mc14_m5')
# mc14_m6 = EpicsMotor('XF:03IDC-ES{MC:14-Ax:6}Mtr', name='mc14_m6')
# mc14_m7 = EpicsMotor('XF:03IDC-ES{MC:14-Ax:7}Mtr', name='mc14_m7')
# mc14_m8 = EpicsMotor('XF:03IDC-ES{MC:14-Ax:8}Mtr', name='mc14_m8')
