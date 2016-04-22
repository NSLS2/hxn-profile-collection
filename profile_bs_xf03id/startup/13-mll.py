import math
from ophyd import (Device, EpicsMotor, Signal, Component as Cpt,
                   PseudoSingle, PseudoPositioner,
                   movr)

from ophyd.pseudopos import (real_position_argument,
                             pseudo_position_argument)


def rename_motors(device):
    from ophyd.positioner import PositionerBase
    cls = device.__class__
    for attribute in device.signal_names:
        motor = getattr(device, attribute)
        component = getattr(cls, attribute)
        if isinstance(motor, PositionerBase):
            if 'name' in component.kwargs:
                motor.name = component.kwargs['name']


class HxnMLLSample(Device):
    fine_x = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-ssx}Mtr', name='ssx')
    fine_y = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-ssy}Mtr', name='ssy')
    fine_z = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-ssz}Mtr', name='ssz')
    th = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:1-Ax:0}Mtr', name='th')

    coarse_x = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:4-Ax:5}Mtr', name='sx')
    coarse_y = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:3-Ax:0}Mtr', name='sy')
    coarse_x1 = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:3-Ax:1}Mtr', name='sx1')
    coarse_z = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:3-Ax:2}Mtr', name='sz')
    coarse_z1 = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:3-Ax:3}Mtr', name='sz1')


smll = HxnMLLSample('', name='smll')
ssx = smll.fine_x
ssy = smll.fine_y
ssz = smll.fine_z
sth = smll.th

sx = smll.coarse_x
sy = smll.coarse_y
sx1 = smll.coarse_x1
sz = smll.coarse_z
sz1 = smll.coarse_z1


# NOTE: normally, motors would be named smll_{attribute} - or smll_coarse_x for
#       example. To quickly rename them to what the component line shows, use
#       'rename_motors' on the device:
rename_motors(smll)


class HxnAnc350_3(Device):
    '''3 axis ANC350'''
    ax0 = Cpt(EpicsMotor, '-Ax:0}Mtr')
    ax1 = Cpt(EpicsMotor, '-Ax:1}Mtr')
    ax2 = Cpt(EpicsMotor, '-Ax:2}Mtr')


class HxnAnc350_4(Device):
    '''4 axis ANC350'''
    ax0 = Cpt(EpicsMotor, '-Ax:0}Mtr')
    ax1 = Cpt(EpicsMotor, '-Ax:1}Mtr')
    ax2 = Cpt(EpicsMotor, '-Ax:2}Mtr')
    ax3 = Cpt(EpicsMotor, '-Ax:3}Mtr')


class HxnAnc350_6(Device):
    '''6 axis ANC350'''
    ax0 = Cpt(EpicsMotor, '-Ax:0}Mtr')
    ax1 = Cpt(EpicsMotor, '-Ax:1}Mtr')
    ax2 = Cpt(EpicsMotor, '-Ax:2}Mtr')
    ax3 = Cpt(EpicsMotor, '-Ax:3}Mtr')
    ax4 = Cpt(EpicsMotor, '-Ax:4}Mtr')
    ax5 = Cpt(EpicsMotor, '-Ax:5}Mtr')


# Note that different controllers have different axis counts:
anc350_1 = HxnAnc350_6('XF:03IDC-ES{ANC350:1', name='anc350_1')
anc350_2 = HxnAnc350_6('XF:03IDC-ES{ANC350:2', name='anc350_2')
anc350_3 = HxnAnc350_4('XF:03IDC-ES{ANC350:3', name='anc350_3')
anc350_4 = HxnAnc350_6('XF:03IDC-ES{ANC350:4', name='anc350_4')
anc350_5 = HxnAnc350_6('XF:03IDC-ES{ANC350:5', name='anc350_5')
anc350_6 = HxnAnc350_6('XF:03IDC-ES{ANC350:6', name='anc350_6')
anc350_7 = HxnAnc350_3('XF:03IDC-ES{ANC350:7', name='anc350_7')


class HxnVerticalMLL(Device):
    x = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:2-Ax:0}Mtr', name='vx')
    y = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:2-Ax:1}Mtr', name='vy')
    z = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:2-Ax:2}Mtr', name='vz')
    chi = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:2-Ax:4}Mtr', name='vchi')
    th = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:2-Ax:3}Mtr', name='vth')


vmll = HxnVerticalMLL('', name='vmll')
rename_motors(vmll)


class HxnHorizontalMLL(Device):
    x = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:2-Ax:5}Mtr', name='hx')
    y = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:4-Ax:0}Mtr', name='hy')
    z = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:4-Ax:1}Mtr', name='hz')
    th = EpicsMotor('XF:03IDC-ES{ANC350:1-Ax:2}Mtr', name='hth')


hmll = HxnHorizontalMLL('', name='hmll')
rename_motors(hmll)


class HxnMLL_OSA(Device):
    osax = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:4-Ax:2}Mtr', name='osax')
    osay = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:4-Ax:3}Mtr', name='osay')
    osaz = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:4-Ax:4}Mtr', name='osaz')


mllosa = HxnMLL_OSA('', name='mllosa')
rename_motors(mllosa)


class HxnMLLBeamStop(Device):
    bsx = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:5-Ax:0}Mtr', name='bsx')
    bsy = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:5-Ax:1}Mtr', name='bsy')


mllbs = HxnMLLBeamStop('', name='mllbs')
rename_motors(mllbs)


class PseudoAngleCorrection(PseudoPositioner):
    '''Pseudo positioner definition for MLL coarse and fine sample positioners
    with angular correction
    '''
    # configuration settings
    theta = Cpt(Signal, value=15.0)

    def __init__(self, prefix, **kwargs):
        super().__init__(prefix, **kwargs)

        # if theta changes, update the pseudo position
        self.theta.subscribe(self.parameter_updated)

    def parameter_updated(self, value=None, **kwargs):
        self._update_position()

    @property
    def radian_theta(self):
        return math.radians(self.theta.get())

    @pseudo_position_argument
    def forward(self, position):
        theta = self.radian_theta
        c = math.cos(theta)
        s = math.sin(theta)

        x = c * position.px + s * position.pz
        z = -s * position.px + c * position.pz
        return self.RealPosition(x=x, z=z)

    @real_position_argument
    def inverse(self, position):
        theta = self.radian_theta
        c = math.cos(theta)
        s = math.sin(theta)
        x = c * position.x - s * position.z
        z = s * position.x + c * position.z
        return self.PseudoPosition(px=x, pz=z)


class PseudoMLLFineSample(PseudoAngleCorrection):
    # pseudo axes
    px = Cpt(PseudoSingle, name='pssx')
    pz = Cpt(PseudoSingle, name='pssz')

    # real axes
    x = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-ssx}Mtr', name='ssx')
    z = Cpt(EpicsMotor, 'XF:03IDC-ES{Ppmac:1-ssz}Mtr', name='ssz')


class PseudoMLLCoarseSample(PseudoAngleCorrection):
    # pseudo axes
    px = Cpt(PseudoSingle, name='psx')
    pz = Cpt(PseudoSingle, name='psz')

    # real axes
    x = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:4-Ax:5}Mtr', name='sx')
    z = Cpt(EpicsMotor, 'XF:03IDC-ES{ANC350:3-Ax:2}Mtr', name='sz')


pseudo_mll_fine = PseudoMLLFineSample('', name='pmc')
pssx = pseudo_mll_fine.x
pssz = pseudo_mll_fine.z
rename_motors(pseudo_mll_fine)
# To tweak the angle, set pseudo_mll_fine.theta.put(15.1) for example


pseudo_mll_coarse = PseudoMLLCoarseSample('', name='pmc')
psx = pseudo_mll_coarse.x
psz = pseudo_mll_coarse.z
rename_motors(pseudo_mll_coarse)


def movr_hth(angle):
    radian = angle * math.pi / 180.0
    correction = -1. * math.tan(radian) * 34376.6
    movr(hmll.th, angle)
    movr(hmll.coarse_x, correction)
