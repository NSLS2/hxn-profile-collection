from ophyd import (EpicsMotor, Device, Component as Cpt)


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


_xz_angle = 15. * pi / 180.


def _pssxz_rev(ssx=None, ssz=None):
    if None in [ssx, ssz]:
        return [0.0, 0.0]

    _pssx = ssx * cos(_xz_angle) + ssz * sin(_xz_angle)
    _pssz = -ssx * sin(_xz_angle) + ssz * cos(_xz_angle)
    return [_pssx, _pssz]


def _pssxz_fwd(pssx=None, pssz=None):
    if None in [pssx, pssz]:
        return [0.0, 0.0]

    _ssx = pssx * cos(_xz_angle) - pssz * sin(_xz_angle)
    _ssz = pssx * sin(_xz_angle) + pssz * cos(_xz_angle)
    return [_ssx, _ssz]


# _pssxz = PseudoPositioner('_pssxz', [ssx, ssz], forward=_pssxz_fwd, reverse=_pssxz_rev,
#                           pseudo=['pssx', 'pssz'])
#
# pssx = _pssxz['pssx']
# pssz = _pssxz['pssz']


def _psxz_rev(sx=None, sz=None):
    if None in [sx, sz]:
        return [0.0, 0.0]
    _psx = sx * cos(_xz_angle) + sz * sin(_xz_angle)
    _psz = -sx * sin(_xz_angle) + sz * cos(_xz_angle)
    return [_psx, _psz]


def _psxz_fwd(psx=None, psz=None):
    if None in [psx, psz]:
        return [0.0, 0.0]
    _sx = psx * cos(_xz_angle) - psz * sin(_xz_angle)
    _sz = psx * sin(_xz_angle) + psz * cos(_xz_angle)
    return [_sx, _sz]


# _psxz = PseudoPositioner('_psxz', [sx, sz], forward=_psxz_fwd, reverse=_psxz_rev,
#                          pseudo=['psx', 'psz'])
# # psx, psz = _psxz
# psx = _psxz['psx']
# psz = _psxz['psz']


def movr_hth(angle):
    radian = angle*pi/180.0
    correction = -1.*tan(radian)*34376.6
    movr(hth, angle)
    movr(hx,correction)


def movr_ssx(d):
    dx_n = d * cos(15.*pi/180.)
    dz_n = d * sin(15.*pi/180.)
    movr(ssx, dx_n)
    movr(ssz, dz_n)

def movr_sx(d):
    dx_n = d * cos(15.*pi/180.)
    dz_n = d * sin(15.*pi/180.)
    movr(sx, dx_n)
    movr(sz, dz_n)


def movr_sz(d):
    dx_n = -d * sin(15.*pi/180.)
    dz_n = d * cos(15.*pi/180.)
    movr(sx, dx_n)
    movr(sz, dz_n)


def movr_ssz(d):
    dx_n = -d * sin(15.*pi/180.)
    dz_n = d * cos(15.*pi/180.)
    movr(ssx, dx_n)
    movr(ssz, dz_n)
