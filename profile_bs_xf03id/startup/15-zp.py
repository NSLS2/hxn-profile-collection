import math

from ophyd import (Device, PVPositionerPC, EpicsMotor, Signal, EpicsSignal,
                   EpicsSignalRO, Component as Cpt, FormattedComponent as FCpt,
                   PseudoSingle, PseudoPositioner,
                   )

from ophyd.pseudopos import (real_position_argument,
                             pseudo_position_argument)


from hxntools.device import NamedDevice


class SmarpodBase(PVPositionerPC):
    actuate = Cpt(EpicsSignal, 'XF:03IDC-ES{SPod:1}Move-Cmd.PROC')
    actuate_value = 1
    done = Cpt(EpicsSignalRO, 'XF:03IDC-ES{SPod:1}Moving-I')
    done_value = 1

    def __init__(self, prefix='', axis=0, **kwargs):
        self.axis = axis
        super().__init__(prefix='', **kwargs)


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


class HxnZPSample(NamedDevice):
    # Zoneplate module fine sample stage axes (closed on cap
    # sensors/interferometer)
    fine_x = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-zpssx}Mtr', name='zpssx')
    fine_y = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-zpssy}Mtr', name='zpssy')
    fine_z = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-zpssz}Mtr', name='zpssz')

    # rotary underneath sample
    theta = Cpt(EpicsMotor, 'XF:03IDC-ES{SC210:1-Ax:1}Mtr', name='zpsth')
    # PI controller underneath smarpod
    coarse_x = Cpt(EpicsMotor, 'XF:03IDC-ES{ZpPI:1-zpsx}Mtr', name='zpsx')
    coarse_z = Cpt(EpicsMotor, 'XF:03IDC-ES{ZpPI:1-zpsz}Mtr', name='zpsz')

    smarx = Cpt(SmarpodTranslationAxis, axis=2, name='smarx')
    smary = Cpt(SmarpodTranslationAxis, axis=3, name='smary')
    smarz = Cpt(SmarpodTranslationAxis, axis=1, name='smary')
    smarthx = Cpt(SmarpodRotationAxis, axis=2, name='smarthx')
    smarthy = Cpt(SmarpodRotationAxis, axis=3, name='smarthy')
    smarthz = Cpt(SmarpodRotationAxis, axis=1, name='smarthy')


zps = HxnZPSample('', name='zps')
zpssx = zps.fine_x
zpssy = zps.fine_y
zpssz = zps.fine_z


class HxnZP_OSA(NamedDevice):
    x = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:5-Ax:0}Mtr', name='zposax')
    y = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:5-Ax:1}Mtr', name='zposay')
    z = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:5-Ax:2}Mtr', name='zposaz')


zposa = HxnZP_OSA('', name='zposa')


class HxnZPBeamStop(NamedDevice):
    x = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:3-Ax:1}Mtr', name='zpbsx')
    y = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:3-Ax:3}Mtr', name='zpbsy')
    z = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:3-Ax:2}Mtr', name='zpbsz')


zpbs = HxnZPBeamStop('', name='zpbs')


class HxnZonePlate(NamedDevice):
    # TPA stage holding the ZP (underneath long travel range stage)
    x = Cpt(EpicsMotor, 'XF:03IDC-ES{ZpTpa-Ax:X}Mtr', name='zpx')
    y = Cpt(EpicsMotor, 'XF:03IDC-ES{ZpTpa-Ax:Y}Mtr', name='zpy')
    z = Cpt(EpicsMotor, 'XF:03IDC-ES{ZpTpa-Ax:Z}Mtr', name='zpz')

    # long travel range z holding the ZP
    long_z = Cpt(EpicsMotor, 'XF:03IDC-ES{MCS:1-Ax:zpz1}Mtr', name='zpz1')


zp = HxnZonePlate('', name='zp')


class FineSampleLabX(PseudoPositioner, NamedDevice):
    '''Pseudo positioner definition for zoneplate fine sample positioner
    with angular correction
    '''
    # pseudo axes
    zpssx_lab = Cpt(PseudoSingle, name='zpssx_lab')
    zpssz_lab = Cpt(PseudoSingle, name='zpssz_lab')

    # real axes
    zpssx = Cpt(EpicsMotor, '{Ppmac:1-zpssx}Mtr', name='zpssx')
    zpssz = Cpt(EpicsMotor, '{Ppmac:1-zpssz}Mtr', name='zpssz')
    theta = Cpt(EpicsMotor, '{SC210:1-Ax:1}Mtr', name='zpsth')

    # configuration settings
    theta0 = Cpt(Signal, value=0.0)

    def __init__(self, prefix, **kwargs):
        super().__init__(prefix, **kwargs)

        # if theta changes, update the pseudo position
        self.theta0.subscribe(self.parameter_updated)

    def parameter_updated(self, value=None, **kwargs):
        self._update_position()

    @property
    def radian_theta(self):
        return math.radians(self.theta.position + self.theta0.get())

    @pseudo_position_argument
    def forward(self, position):
        theta = self.radian_theta
        c = math.cos(theta)
        s = math.sin(theta)

        x = c * position.zpssx_lab + s * position.zpssz_lab
        z = -s * position.zpssx_lab + c * position.zpssz_lab
        return self.RealPosition(zpssx=x, zpssz=z)

    @real_position_argument
    def inverse(self, position):
        theta = self.radian_theta
        c = math.cos(theta)
        s = math.sin(theta)
        x = c * position.zpssx - s * position.zpssz
        z = s * position.zpssx + c * position.zpssz
        return self.PseudoPosition(zpssx_lab=x, zpssz_lab=z)


zplab = FineSampleLabX('XF:03IDC-ES', name='zplab')

zpssx_lab = zplab.zpssx_lab
zpssz_lab = zplab.zpssz_lab
