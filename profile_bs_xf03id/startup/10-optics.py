from ophyd import EpicsMotor

dcm_th = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:Bragg}Mtr', name='dcm_th')
dcm_x = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:X}Mtr', name='dcm_x')
dcm_p = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:P}Mtr', name='dcm_p')
dcm_r = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:R}Mtr', name='dcm_r')
m1_x = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:X}Mtr', name='m1_x')
m1_y = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:Y}Mtr', name='m1_y')
m1_p = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:P}Mtr', name='m1_p')
m1_b = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:Bend}Mtr', name='m1_b')
m2_x = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:X}Mtr', name='m2_x')
m2_y = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:Y}Mtr', name='m2_y')
m2_p = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:P}Mtr', name='m2_p')
m2_b = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:Bend}Mtr', name='m2_b')
s1_bot = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Btm}Mtr', name='s1_bot')
s1_top = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Top}Mtr', name='s1_top')
s1_inb = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Inb}Mtr', name='s1_inb')
s1_outb = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Outb}Mtr', name='s1_outb')
s2_bot = EpicsMotor('XF:03IDA-OP{Slt:2-Ax:Btm}Mtr', name='s2_bot')
s2_top = EpicsMotor('XF:03IDA-OP{Slt:2-Ax:Top}Mtr', name='s2_top')
s2_inb = EpicsMotor('XF:03IDA-OP{Slt:2-Ax:Inb}Mtr', name='s2_inb')
s2_outb = EpicsMotor('XF:03IDA-OP{Slt:2-Ax:Outb}Mtr', name='s2_outb')
bpm1_x = EpicsMotor('XF:03IDA-OP{BPM:1-Ax:X}Mtr', name='bpm1_x')
bpm1_y = EpicsMotor('XF:03IDA-OP{BPM:1-Ax:Y}Mtr', name='bpm1_y')
bpm2_y = EpicsMotor('XF:03IDA-OP{BPM:2-Ax:Y}Mtr', name='bpm2_y')
bpm2_x = EpicsMotor('XF:03IDA-OP{BPM:2-Ax:X}Mtr', name='bpm2_x')
bpm3_x = EpicsMotor('XF:03IDA-OP{BPM:3-Ax:X}Mtr', name='bpm3_x')
fs1_y = EpicsMotor('XF:03IDA-OP{FS:1-Ax:Y}Mtr', name='fs1_y')
bpm4_y = EpicsMotor('XF:03IDA-OP{BPM:4-Ax:Y}Mtr', name='bpm4_y')
bpm5_y = EpicsMotor('XF:03IDA-OP{BPM:5-Ax:Y}Mtr', name='bpm5_y')
fl1_y = EpicsMotor('XF:03IDA-OP{Flr:1-Ax:Y}Mtr', name='fl1_y')
fl2_y = EpicsMotor('XF:03IDA-OP{Flr:2-Ax:Y}Mtr', name='fl2_y')
crl_x = EpicsMotor('XF:03IDA-OP{Lens:CRL-Ax:X}Mtr', name='crl_x')
crl_y = EpicsMotor('XF:03IDA-OP{Lens:CRL-Ax:Y}Mtr', name='crl_y')
crl_p = EpicsMotor('XF:03IDA-OP{Lens:CRL-Ax:P}Mtr', name='crl_p')
m1_pf = EpicsMotor('XF:03IDA-OP{HCM:1-Ax:PF}Mtr', name='m1_pf')
m2_pf = EpicsMotor('XF:03IDA-OP{HFM:1-Ax:PF}Mtr', name='m2_pf')
dcm_pf = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:PF}Mtr', name='dcm_pf')
dcm_rf = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:RF}Mtr', name='dcm_rf')


# # M1
# m1x = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:X}Mtr', name='m1x')
# m1y = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:Y}Mtr', name='m1y')
# m1p = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:P}Mtr', name='m1p')
# m1b = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:Bend}Mtr', name='m1b')
#
# # HCM
# m1pz = EpicsMotor('XF:03IDA-OP{HCM:1-Ax:PF}Mtr', name='m1pz')
#
# # slit:1
# slt = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Top}Mtr', name='slt')
# slb = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Btm}Mtr', name='slb')
# sli = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Inb}Mtr', name='sli')
# slo = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Outb}Mtr', name='slo')
#
# # Mon:1
# bragg = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:Bragg}Mtr', name='bragg')
#
# # HFM
# m2x = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:X}Mtr', name='m2x')
# m2y = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:Y}Mtr', name='m2y')
# m2p = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:P}Mtr', name='m2p')
# m2b = EpicsMotor('XF:03IDA-OP{Mir:2-Ax:Bend}Mtr', name='m2b')
#
#
# # Diagnostic Manipulators
# f1y = EpicsMotor('XF:03IDA-OP{Flr:1-Ax:Y}Mtr', name='f1y')
# f2y = EpicsMotor('XF:03IDA-OP{Flr:2-Ax:Y}Mtr', name='f2y')
#
# # CRL
# crlth = EpicsMotor('XF:03IDA-OP{Lens:CRL-Ax:P}Mtr', name='crlth')
#
# # nanoBPM2@SSA1
# nano2y = EpicsMotor('XF:03IDB-OP{BPM:6-Ax:Y}Mtr', name='nano2y')
#
# # nanoBPM3@SSA2
# nano3y = EpicsMotor('XF:03IDC-OP{BPM:7-Ax:Y}Mtr', name='nano3y')
#
#
