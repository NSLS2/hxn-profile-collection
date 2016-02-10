from ophyd.controls import (PVPositioner, EpicsMotor)

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


# Smarpod axes - updated to support put completion

smarpod_settings = dict(act='XF:03IDC-ES{SPod:1}Move-Cmd.PROC',
                        act_val=1,
                        done='XF:03IDC-ES{SPod:1}Moving-I',
                        done_val=1,
                        put_complete=True,
                        )


smarx = PVPositioner(setpoint='XF:03IDC-ES{SPod:1-Ax:2}Pos-SP',
                     readback='XF:03IDC-ES{SPod:1-Ax:2}Pos-I',
                     name='smarx',
                     **smarpod_settings
                     )

smary = PVPositioner(setpoint='XF:03IDC-ES{SPod:1-Ax:3}Pos-SP',
                     readback='XF:03IDC-ES{SPod:1-Ax:3}Pos-I',
                     name='smary',
                     **smarpod_settings
                     )

smarz = PVPositioner(setpoint='XF:03IDC-ES{SPod:1-Ax:1}Pos-SP',
                     readback='XF:03IDC-ES{SPod:1-Ax:1}Pos-I',
                     name='smarz',
                     **smarpod_settings
                     )

smarthx = PVPositioner(setpoint='XF:03IDC-ES{SPod:1-Ax:2}Rot-SP',
                       readback='XF:03IDC-ES{SPod:1-Ax:2}Rot-I',
                       name='smarthx',
                       **smarpod_settings
                       )

smarthy = PVPositioner(setpoint='XF:03IDC-ES{SPod:1-Ax:3}Rot-SP',
                       readback='XF:03IDC-ES{SPod:1-Ax:3}Rot-I',
                       name='smarthy',
                       **smarpod_settings
                       )

smarthz = PVPositioner(setpoint='XF:03IDC-ES{SPod:1-Ax:1}Rot-SP',
                       readback='XF:03IDC-ES{SPod:1-Ax:1}Rot-I',
                       name='smarthz',
                       **smarpod_settings
                       )
