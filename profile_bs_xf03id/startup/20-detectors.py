from ophyd import (EpicsSignal, EpicsSignalRO)
from ophyd import (Device, Component as Cpt)

import hxntools.handlers
from hxntools.detectors import (HxnTimepixDetector, HxnMerlinDetector,
                                BeamStatusDetector, HxnMercuryDetector)
from hxntools.detectors.zebra import HxnZebra

# Register all HXN-specific handlers so that filestore can load all detector
# spectra and images directly:
hxntools.handlers.register()

# - 2D pixel array detectors
# -- Timepix 1
timepix1 = HxnTimepixDetector('XF:03IDC-ES{Tpx:1}', name='timepix1',
                              image_name='timepix1',
                              read_attrs=['hdf5', 'cam'])
timepix1.hdf5.read_attrs = []

# -- Timepix 2
timepix2 = HxnTimepixDetector('XF:03IDC-ES{Tpx:2}', name='timepix2',
                              image_name='timepix1',
                              read_attrs=['hdf5', 'cam'])
timepix2.hdf5.read_attrs = []

# -- Merlin 1
merlin1 = HxnMerlinDetector('XF:03IDC-ES{Merlin:1}', name='merlin1',
                            image_name='merlin1',
                            read_attrs=['hdf5', 'cam'])
merlin1.hdf5.read_attrs = []

# - Other detectors and triggering devices
# -- DXP Mercury (1 channel)
mercury1 = HxnMercuryDetector('XF:03IDC-ES{DXP:1}', name='mercury1')
mercury1.read_attrs = ['dxp', 'mca']
mercury1.dxp.read_attrs = []

# -- Quantum Detectors Zebra
zebra = HxnZebra('XF:03IDC-ES{Zeb:1}:', name='zebra')
zebra.read_attrs = []

# ugap scan trigger
ugap_trig = EpicsSignal('SR:C3-ID:G1{IVU20:1-Mtr:2}Sw:Go', name='ugap_trig')


# -- Lakeshores
class HxnLakeShore(Device):
    ch_a = Cpt(EpicsSignalRO, '-Ch:A}C:T-I')
    ch_b = Cpt(EpicsSignalRO, '-Ch:B}C:T-I')
    ch_c = Cpt(EpicsSignalRO, '-Ch:C}C:T-I')
    ch_d = Cpt(EpicsSignalRO, '-Ch:D}C:T-I')

    def set_names(self, cha, chb, chc, chd):
        '''Set names of all channels

        Returns channel signals
        '''
        self.ch_a.name = cha
        self.ch_b.name = chb
        self.ch_c.name = chc
        self.ch_d.name = chd
        return self.ch_a, self.ch_b, self.ch_c, self.ch_d


lakeshore2 = HxnLakeShore('XF:03IDC-ES{LS:2', name='lakeshore2')

# Name the lakeshore channels:
t_hlens, t_vlens, t_sample, t_base = lakeshore2.set_names(
    't_hlens', 't_vlens', 't_sample', 't_base')

# X-ray eye camera sigma X/sigma Y
sigx = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaX_RBV', name='sigx')
sigy = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaY_RBV', name='sigy')


# -- Interferometers

class HxnFPSensor(Device):
    ch0 = Cpt(EpicsSignalRO, '-Chan0}Pos-I')
    ch1 = Cpt(EpicsSignalRO, '-Chan1}Pos-I')
    ch2 = Cpt(EpicsSignalRO, '-Chan2}Pos-I')

    def set_names(self, ch0, ch1, ch2):
        '''Set names of all channels

        Returns channel signals
        '''
        self.ch0.name = ch0
        self.ch1.name = ch1
        self.ch2.name = ch2
        return self.ch0, self.ch1, self.ch2


fpsensor_1 = HxnFPSensor('XF:03IDC-ES{FPS:1', name='fpsensor_1')
fpsensor_2 = HxnFPSensor('XF:03IDC-ES{FPS:2', name='fpsensor_2')
fpsensor_3 = HxnFPSensor('XF:03IDC-ES{FPS:3', name='fpsensor_3')
fpsensor_4 = HxnFPSensor('XF:03IDC-ES{FPS:4', name='fpsensor_4')
fpsensor_5 = HxnFPSensor('XF:03IDC-ES{FPS:5', name='fpsensor_5')
fpsensor_6 = HxnFPSensor('XF:03IDC-ES{FPS:6', name='fpsensor_6')

int_sx, int_sy, int_sz = fpsensor_1.set_names('int_sx',
                                              'int_sy',
                                              'int_sz')
int_hx, int_vy, int_hy = fpsensor_2.set_names('int_hx',
                                              'int_vy',
                                              'int_hy')
int_hz, int_vx, int_vz = fpsensor_3.set_names('int_hz',
                                              'int_vx',
                                              'int_vz')
int_zpssx, int_zpssy, int_zpssz = fpsensor_4.set_names('int_zpssx',
                                                       'int_zpssy',
                                                       'int_zpssz')
int_zpx1, int_zpx2, int_zpy1 = fpsensor_5.set_names('int_zpx1',
                                                    'int_zpx2',
                                                    'int_zpy1')
int_zpy2, int_zpz, int_zpspare1 = fpsensor_6.set_names('int_zpy2',
                                                       'int_zpz',
                                                       'int_zpspare1')

# Front-end Xray BPMs and local bumps
class HxnBpm(Device):
    x = Cpt(EpicsSignalRO, 'Pos:X-I')
    y = Cpt(EpicsSignalRO, 'Pos:Y-I')


xbpm = HxnBpm('SR:C03-BI{XBPM:1}', name='xbpm')

angle_x = EpicsSignalRO('SR:C31-{AI}Aie3:Angle-x-Cal', name='angle_x')
angle_y = EpicsSignalRO('SR:C31-{AI}Aie3:Angle-y-Cal', name='angle_y')

# Diamond Quad BPMs in C hutch
# quad = HxnBpm('SR:C12-BI{XBPM:1}', name='quad')


sr_shutter_status = EpicsSignalRO('SR-EPS{PLC:1}Sts:MstrSh-Sts',
                                  name='sr_shutter_status')
sr_beam_current = EpicsSignalRO('SR:C03-BI{DCCT:1}I:Real-I',
                                name='sr_beam_current')

det_beamstatus = BeamStatusDetector(min_current=100.0)
