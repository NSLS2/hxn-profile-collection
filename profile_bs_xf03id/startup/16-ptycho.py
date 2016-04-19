from ophyd import (EpicsMotor, Device, Component as Cpt,
                   PVPositioner)



class HxnPrototypeMicroscope(Device):
    vx = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:1}Mtr')
    vy = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:2}Mtr')
    vt = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:3}Mtr')
    vchi = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:4}Mtr')
    osat = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:5}Mtr')
    osay = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:6}Mtr')
    osax = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:7}Mtr')
    osaz = Cpt(EpicsMotor, 'XF:03IDC-ES{Proto:1-Ax:8}Mtr')

    ssx = Cpt(EpicsMotor, 'XF:03IDC-ES{Ddrive:1-Ax:2}Mtr')
    ssy = Cpt(EpicsMotor, 'XF:03IDC-ES{Ddrive:1-Ax:3}Mtr')
    ssz = Cpt(EpicsMotor, 'XF:03IDC-ES{Ddrive:1-Ax:1}Mtr')

    vz = Cpt(EpicsMotor, 'XF:03IDC-ES{MMC100:1-Ax:1}Mtr')
    cz = Cpt(EpicsMotor, 'XF:03IDC-ES{MMC100:1-Ax:2}Mtr')
    cx = Cpt(EpicsMotor, 'XF:03IDC-ES{MMC100:1-Ax:3}Mtr')
    cy = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:6-Ax:1}Mtr')

    bsx = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:3-Ax:1}Mtr')
    bsy = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:3-Ax:2}Mtr')
    bsz = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:3-Ax:3}Mtr')

    vth = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:4-Ax:1}Mtr')


p = HxnPrototypeMicroscope('', name='p')

p_vx = p.vx
p_vy = p.vy
p_vt = p.vt
p_vchi = p.vchi
p_osat = p.osat
p_osay = p.osay
p_osax = p.osax
p_osaz = p.osaz

p_ssx = p.ssx
p_ssy = p.ssy
p_ssz = p.ssz

p_vz = p.vz
p_cz = p.cz
p_cx = p.cx
p_cy = p.cy


p_bsx = p.bsx
p_bsy = p.bsy
p_bsz = p.bsz

p_vth = p.vth
