from ophyd import (PVPositionerPC, EpicsMotor,
                   EpicsSignal, EpicsSignalRO,
                   Component as Cpt,
                   FormattedComponent as FCpt,
                   )

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


class SmarpodBase(PVPositionerPC):
    actuate = Cpt(EpicsSignal, 'XF:03IDC-ES{SPod:1}Move-Cmd.PROC')
    actuate_value = 1
    done = Cpt(EpicsSignalRO, 'XF:03IDC-ES{SPod:1}Moving-I')
    done_value = 1

    def __init__(self, prefix, axis,  **kwargs):
        self.axis = axis
        super().__init__(prefix, **kwargs)


class SmarpodTranslationAxis(SmarpodBase):
    setpoint = FCpt(EpicsSignal,
                    'XF:03IDC-ES{{SPod:1-Ax:{self.axis}}}Pos-SP')
    readback = FCpt(EpicsSignal,
                    'XF:03IDC-ES{{SPod:1-Ax:{self.axis}}}Pos-I')


class SmarpodRotationAxis(SmarpodBase):
    setpoint = FCpt(EpicsSignal,
                    'XF:03IDC-ES{{SPod:1-Ax:{self.axis}}}Rot-SP')
    readback = FCpt(EpicsSignal,
                    'XF:03IDC-ES{{SPod:1-Ax:{self.axis}}}Rot-I')


smarx = SmarpodTranslationAxis('', axis=2, name='smarx')
smary = SmarpodTranslationAxis('', axis=3, name='smary')
smarz = SmarpodTranslationAxis('', axis=1, name='smary')
smarthx = SmarpodRotationAxis('', axis=2, name='smarthx')
smarthy = SmarpodRotationAxis('', axis=3, name='smarthy')
smarthz = SmarpodRotationAxis('', axis=1, name='smarthy')
