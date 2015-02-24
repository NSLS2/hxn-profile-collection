from ophyd.controls import EpicsMotor, PVPositioner

# M1A

kwargs = {'act': 'XF:23IDA-OP:1{Mir:1}MOVE_CMD.PROC',
          'act_val': 1,
          'stop': 'XF:23IDA-OP:1{Mir:1}STOP_CMD.PROC',
          'stop_val': 1,
          'done': 'XF:23IDA-OP:1{Mir:1}BUSY_STS',
          'done_val': 0}

m1a_z = PVPositioner('XF:23IDA-OP:1{Mir:1-Ax:Z}Mtr_POS_SP',
                     readback='XF:23IDA-OP:1{Mir:1-Ax:Z}Mtr_MON',
                     name='m1a_z', **kwargs)
m1a_y = PVPositioner('XF:23IDA-OP:1{Mir:1-Ax:Y}Mtr_POS_SP',
                     readback='XF:23IDA-OP:1{Mir:1-Ax:Y}Mtr_MON',
                     name='m1a_y', **kwargs)
m1a_x = PVPositioner('XF:23IDA-OP:1{Mir:1-Ax:X}Mtr_POS_SP',
                     readback='XF:23IDA-OP:1{Mir:1-Ax:X}Mtr_MON',
                     name='m1a_x', **kwargs)
m1a_pit = PVPositioner('XF:23IDA-OP:1{Mir:1-Ax:Pit}Mtr_POS_SP',
                       readback='XF:23IDA-OP:1{Mir:1-Ax:Pit}Mtr_MON',
                       name='m1a_pit', **kwargs)
m1a_yaw = PVPositioner('XF:23IDA-OP:1{Mir:1-Ax:Yaw}Mtr_POS_SP',
                       readback='XF:23IDA-OP:1{Mir:1-Ax:Yaw}Mtr_MON',
                       name='m1a_yaw', **kwargs)
m1a_rol = PVPositioner('XF:23IDA-OP:1{Mir:1-Ax:Rol}Mtr_POS_SP',
                       readback='XF:23IDA-OP:1{Mir:1-Ax:Rol}Mtr_MON',
                       name='m1a_rol', **kwargs)

# VLS-PGM

pgm_energy = PVPositioner('XF:23ID1-OP{Mono}Enrgy-SP',
                          readback='XF:23ID1-OP{Mono}Enrgy-I',
                          stop='XF:23ID1-OP{Mono}Cmd:Stop-Cmd',
                          stop_val=1, put_complete=True,
                          name='pgm_energy')

pgm_mir_pit = EpicsMotor('XF:23ID1-OP{Mono-Ax:MirP}Mtr', name='pgm_mir_pit')
pgm_grt_pit = EpicsMotor('XF:23ID1-OP{Mono-Ax:GrtP}Mtr', name='pgm_grt_pit')
pgm_mir_x = EpicsMotor('XF:23ID1-OP{Mono-Ax:MirX}Mtr', name='pgm_mir_x')
pgm_grt_x = EpicsMotor('XF:23ID1-OP{Mono-Ax:GrtX}Mtr', name='pgm_grt_x')

# M3A Mirror

m3a_x = EpicsMotor('XF:23ID1-OP{Mir:3-Ax:XAvg}Mtr', name='m3a_x')
m3a_pit = EpicsMotor('XF:23ID1-OP{Mir:3-Ax:P}Mtr',   name='m3a_pit')
m3a_bdr = EpicsMotor('XF:23ID1-OP{Mir:3-Ax:Bdr}Mtr',  name='m3a_bdr')

# Fast CCD Shutter

sh_y = EpicsMotor('XF:23ID1-OP{Sh:Fast-Ax:Y}Mtr', name='sh_y')
sh_x = EpicsMotor('XF:23ID1-OP{Sh:Fast-Ax:X}Mtr', name='sh_x')

# Slits

slt1_xg = EpicsMotor('XF:23ID1-OP{Slt:1-Ax:XGap}Mtr', name='slt1_xg')
slt1_xc = EpicsMotor('XF:23ID1-OP{Slt:1-Ax:XCtr}Mtr', name='slt1_xc')
slt1_yg = EpicsMotor('XF:23ID1-OP{Slt:1-Ax:YGap}Mtr', name='slt1_yg')
slt1_yc = EpicsMotor('XF:23ID1-OP{Slt:1-Ax:YCtr}Mtr', name='slt1_yc')

slt2_xg = EpicsMotor('XF:23ID1-OP{Slt:2-Ax:XGap}Mtr', name='slt2_xg')
slt2_xc = EpicsMotor('XF:23ID1-OP{Slt:2-Ax:XCtr}Mtr', name='slt2_xc')
slt2_yg = EpicsMotor('XF:23ID1-OP{Slt:2-Ax:YGap}Mtr', name='slt2_yg')
slt2_yc = EpicsMotor('XF:23ID1-OP{Slt:2-Ax:YCtr}Mtr', name='slt2_yc')

slt3_x = EpicsMotor('XF:23ID1-OP{Slt:3-Ax:X}Mtr', name='slt3_x')
slt3_y = EpicsMotor('XF:23ID1-OP{Slt:3-Ax:Y}Mtr', name='slt3_y')

# Diagnostic Manipulators

diag2_y = EpicsMotor('XF:23ID1-BI{Diag:2-Ax:Y}Mtr', name='diag2_y')
diag3_y = EpicsMotor('XF:23ID1-BI{Diag:3-Ax:Y}Mtr', name='diag3_y')
diag5_y = EpicsMotor('XF:23ID1-BI{Diag:5-Ax:Y}Mtr', name='diag5_y')
diag6_y = EpicsMotor('XF:23ID1-BI{Diag:6-Ax:Y}Mtr', name='diag6_y')
