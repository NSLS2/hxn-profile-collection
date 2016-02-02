# Zoneplate module fine sample stage axes (closed on cap sensors/interferometer)
zpssx = EpicsMotor('XF:03IDC-ES{Ppmac:1-zpssx}Mtr', name='zpssx')
zpssy = EpicsMotor('XF:03IDC-ES{Ppmac:1-zpssy}Mtr', name='zpssy')
zpssz = EpicsMotor('XF:03IDC-ES{Ppmac:1-zpssz}Mtr', name='zpssz')

# OSA
zposax = EpicsMotor('XF:03IDC-ES{ANC350:5-Ax:0}Mtr', name='zposax')
zposay = EpicsMotor('XF:03IDC-ES{ANC350:5-Ax:1}Mtr', name='zposay')
zposaz = EpicsMotor('XF:03IDC-ES{ANC350:5-Ax:2}Mtr', name='zposaz')

# beamstop
zpbsx = EpicsMotor('XF:03IDC-ES{MCS:3-Ax:1}Mtr', name='zpbsx')
zpbsz = EpicsMotor('XF:03IDC-ES{MCS:3-Ax:2}Mtr', name='zpbsz')
zpbsy = EpicsMotor('XF:03IDC-ES{MCS:3-Ax:3}Mtr', name='zpbsy')

# rotary underneath sample
zpsth = EpicsMotor('XF:03IDC-ES{SC210:1-Ax:1}Mtr', name='zpsth')
# PI controller underneath smarpod
zpsx = EpicsMotor('XF:03IDC-ES{ZpPI:1-zpsx}Mtr', name='zpsx')
zpsz = EpicsMotor('XF:03IDC-ES{ZpPI:1-zpsz}Mtr', name='zpsz')

# TPA stage holding the ZP (underneath long travel range stage)
zpx = EpicsMotor('XF:03IDC-ES{ZpTpa-Ax:X}Mtr', name='zpx')
zpy = EpicsMotor('XF:03IDC-ES{ZpTpa-Ax:Y}Mtr', name='zpy')
zpz = EpicsMotor('XF:03IDC-ES{ZpTpa-Ax:Z}Mtr', name='zpz')

# long travel range z holding the ZP
zpz1 = EpicsMotor('XF:03IDC-ES{MCS:1-Ax:zpz1}Mtr', name='zpz1')
