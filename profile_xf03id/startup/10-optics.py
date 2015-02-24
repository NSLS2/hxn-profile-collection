from ophyd.controls import EpicsMotor, PVPositioner

# M1
m1x = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:X}Mtr', name='m1x')
m1y = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:Y}Mtr', name='m1y')
m1p = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:P}Mtr', name='m1p')
m1b = EpicsMotor('XF:03IDA-OP{Mir:1-Ax:Bend}Mtr', name='m1b')

# HCM
m1pz = EpicsMotor('XF:03IDA-OP{HCM:1-Ax:PF}Mtr', name='m1pz')

# slit:1
slt = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Top}Mtr', name='slt')
slb = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Btm}Mtr', name='slb')
sli = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Inb}Mtr', name='sli')
slo = EpicsMotor('XF:03IDA-OP{Slt:1-Ax:Outb}Mtr', name='slo')

# Mon:1
bragg = EpicsMotor('XF:03IDA-OP{Mon:1-Ax:Bragg}Mtr', name='bragg')


# Diagnostic Manipulators
f1y = EpicsMotor('XF:03IDA-OP{Flr:1-Ax:Y}Mtr', name='f1y')
f2y = EpicsMotor('XF:03IDA-OP{Flr:2-Ax:Y}Mtr', name='f2y')
