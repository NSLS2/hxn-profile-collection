from ophyd.controls import EpicsMotor, PVPositioner, PseudoPositioner

p_ssx = EpicsMotor('XF:03IDC-ES{Ddrive:1-Ax:1}Mtr', name='p_ssx')
p_ssy = EpicsMotor('XF:03IDC-ES{Ddrive:1-Ax:2}Mtr', name='p_ssy')
p_ssz = EpicsMotor('XF:03IDC-ES{Ddrive:1-Ax:3}Mtr', name='p_ssz')

p_vz = EpicsMotor('XF:03IDC-ES{MMC100:1-Ax:1}Mtr', name='p_vz')
p_cz = EpicsMotor('XF:03IDC-ES{MMC100:1-Ax:2}Mtr', name='p_cz')
p_cx = EpicsMotor('XF:03IDC-ES{MMC100:1-Ax:3}Mtr', name='p_cx')

p_bsx = EpicsMotor('XF:03IDC-ES{MCS:3-Ax:1}Mtr', name='p_bsx')
p_bsy = EpicsMotor('XF:03IDC-ES{MCS:3-Ax:2}Mtr', name='p_bsy')
p_bsz = EpicsMotor('XF:03IDC-ES{MCS:3-Ax:3}Mtr', name='p_bsz')

p_vth = EpicsMotor('XF:03IDC-ES{MCS:4-Ax:1}Mtr', name='p_vth')

p_vx = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:1}Mtr', name='p_vx')
p_vy = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:2}Mtr', name='p_vy')
p_vt = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:3}Mtr', name='p_vt')
p_vchi = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:4}Mtr', name='p_vchi')
p_osat = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:5}Mtr', name='p_osat')
p_osay = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:6}Mtr', name='p_osay')
p_osax = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:7}Mtr', name='p_osax')
p_osaz = EpicsMotor('XF:03IDC-ES{Proto:1-Ax:8}Mtr', name='p_osaz')
