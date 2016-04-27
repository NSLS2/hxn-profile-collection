from ophyd import (EpicsSignal, EpicsSignalRO)
from ophyd import (Device, Component as Cpt)
import pandas as pd

from hxntools.detectors import (TimepixDetector, HxnMerlinDetector,
                                BeamStatusDetector)
from hxntools.detectors.zebra import HxnZebra
from hxntools.struck_scaler import (HxnScaler, HxnTriggeringScaler,
                                    StruckScaler)


# Flyscan results are shown using pandas. Maximum rows/columns to use when
# printing the table:
pd.options.display.width = 180
pd.options.display.max_rows = None
pd.options.display.max_columns = 10

# timepix1 = TimepixDetector('XF:03IDC-ES{Tpx:1}', files=['TIFF1:'],
#                            name='timepix1',
#                            file_path='/data', ioc_file_path='/data')
# timepix2 = TimepixDetector('XF:03IDC-ES{Tpx:2}', files=['TIFF1:'],
#                            name='timepix2',
#                            file_path='/data', ioc_file_path='/data')
merlin1 = HxnMerlinDetector('XF:03IDC-ES{Merlin:1}', name='merlin1')
merlin1.tiff1.read_attrs = []

zebra = HxnZebra('XF:03IDC-ES{Zeb:1}:', name='zebra')

# 3IDC RG:C4 VME scalers
# - sclr1 is used for data acquisition. HxnScaler takes care of setting that
#   up:
sclr1 = HxnTriggeringScaler('XF:03IDC-ES{Sclr:1}', name='sclr1')
# let the scans know which detectors sclr1 triggers:
sclr1.scan_type_triggers['step'] = [zebra, merlin1, ]
sclr1.scan_type_triggers['fly'] = []
# - sclr2, on the other hand, is just used as a regular scaler, however the
#   user desires
sclr2 = StruckScaler('XF:03IDC-ES{Sclr:2}', name='sclr2')

n_scaler_mca = 8
sclr1_mca = [sclr1.mca_by_index[i] for i in range(1, n_scaler_mca + 1)]

for mca in sclr1_mca:
    mca.name = 'sclr1_mca{}'.format(mca.index)


# ugap scan trigger
ugap_trig = EpicsSignal('SR:C3-ID:G1{IVU20:1-Mtr:2}Sw:Go', name='ugap_trig')


# Ion chamber
sclr1_ch1 = sclr1.channels.chan1
sclr1_ch2 = sclr1.channels.chan2
sclr1_ch3 = sclr1.channels.chan3
sclr1_ch4 = sclr1.channels.chan4
sclr1_ch5 = sclr1.channels.chan5

sclr2_ch1 = sclr2.channels.chan1
sclr2_ch2 = sclr2.channels.chan2
sclr2_ch3 = sclr2.channels.chan3
sclr2_ch4 = sclr2.channels.chan4
sclr2_ch5 = sclr2.channels.chan5

sclr1_ch4_calc = sclr1.calculations.calc4.value
sclr2_ch4_calc = sclr2.calculations.calc4


def setup_scaler_names(scaler):
    scaler_name = scaler.name
    for i in range(1, 32):
        channel = getattr(scaler.channels, 'chan{}'.format(i))
        channel.name = '{}_ch{}'.format(scaler_name, i)

    for i in range(1, 9):
        c = getattr(scaler.calculations, 'calc{}'.format(i))
        c.name = '{}_ch{}_calc'.format(scaler_name, i)


setup_scaler_names(sclr1)
setup_scaler_names(sclr2)

t_base = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:D}C:T-I', name='t_base')
t_sample = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:C}C:T-I', name='t_sample')
t_vlens = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:B}C:T-I', name='t_vlens')
t_hlens = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:A}C:T-I', name='t_hlens')

# X-ray eye camera sigma X/sigma Y
sigx = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaX_RBV', name='sigx')
sigy = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaY_RBV', name='sigy')

# Interferometers
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
