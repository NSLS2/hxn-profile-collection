from ophyd import (EpicsSignal, EpicsSignalRO)
from ophyd import (Device, Component as Cpt)
import pandas as pd

import hxntools.handlers
from hxntools.detectors import (HxnTimepixDetector, HxnMerlinDetector,
                                BeamStatusDetector)
from hxntools.detectors.zebra import HxnZebra
from hxntools.struck_scaler import (HxnTriggeringScaler, StruckScaler)

# Register all HXN-specific handlers so that filestore can load all detector
# spectra and images directly:
hxntools.handlers.register()

# Flyscan results are shown using pandas. Maximum rows/columns to use when
# printing the table:
pd.options.display.width = 180
pd.options.display.max_rows = None
pd.options.display.max_columns = 10

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

zebra = HxnZebra('XF:03IDC-ES{Zeb:1}:', name='zebra')
zebra.read_attrs = []

# - 3IDC RG:C4 VME scalers
# -- scaler 1 is used for data acquisition. HxnScaler takes care of setting
#    that up:
sclr1 = HxnTriggeringScaler('XF:03IDC-ES{Sclr:1}', name='sclr1')
# let the scans know which detectors sclr1 triggers:
sclr1.scan_type_triggers['step'] = [zebra, merlin1, ]
sclr1.scan_type_triggers['fly'] = []
sclr1.read_attrs = ['channels.chan1', 'channels.chan2', 'channels.chan3',
                    'channels.chan4', 'channels.chan5', 'channels.chan6',
                    'channels.chan7',
                    ]

sclr1_ch1 = sclr1.channels.chan1
sclr1_ch2 = sclr1.channels.chan2
sclr1_ch3 = sclr1.channels.chan3
sclr1_ch4 = sclr1.channels.chan4
sclr1_ch5 = sclr1.channels.chan5

sclr1_ch4_calc = sclr1.calculations.calc4.value

n_scaler_mca = 8
sclr1_mca = [sclr1.mca_by_index[i] for i in range(1, n_scaler_mca + 1)]

for mca in sclr1_mca:
    mca.name = 'sclr1_mca{}'.format(mca.index)

# -- scaler 2, on the other hand, is just used as a regular scaler, however the
#    user desires
sclr2 = StruckScaler('XF:03IDC-ES{Sclr:2}', name='sclr2')
sclr2.read_attrs = sclr1.read_attrs

sclr2_ch1 = sclr2.channels.chan1
sclr2_ch2 = sclr2.channels.chan2
sclr2_ch3 = sclr2.channels.chan3
sclr2_ch4 = sclr2.channels.chan4
sclr2_ch5 = sclr2.channels.chan5

sclr2_ch4_calc = sclr2.calculations.calc4


# Rename all scaler channels according to a common format.
def setup_scaler_names(scaler, channel_format, calc_format):
    for i in range(1, 33):
        channel = getattr(scaler.channels, 'chan{}'.format(i))
        channel.name = channel_format.format(i)

    for i in range(1, 9):
        c = getattr(scaler.calculations, 'calc{}'.format(i))
        c.name = calc_format.format(i)


setup_scaler_names(sclr1, 'sclr1_ch{}', 'sclr1_ch{}_calc')
setup_scaler_names(sclr2, 'sclr2_ch{}', 'sclr2_ch{}_calc')

# ugap scan trigger
ugap_trig = EpicsSignal('SR:C3-ID:G1{IVU20:1-Mtr:2}Sw:Go', name='ugap_trig')


class HxnLakeShore(Device):
    ch_a = Cpt(EpicsSignalRO, '-Ch:A}C:T-I')
    ch_b = Cpt(EpicsSignalRO, '-Ch:B}C:T-I')
    ch_c = Cpt(EpicsSignalRO, '-Ch:C}C:T-I')
    ch_d = Cpt(EpicsSignalRO, '-Ch:D}C:T-I')


lakeshore2 = HxnLakeShore('XF:03IDC-ES{LS:2', name='lakeshore2')

# Name the lakeshore channels:
t_base = lakeshore2.ch_d
t_base.name = 't_base'

t_sample = lakeshore2.ch_c
t_sample.name = 't_sample'

t_vlens = lakeshore2.ch_b
t_vlens.name = 't_vlens'

t_hlens = lakeshore2.ch_a
t_hlens.name = 't_hlens'

# X-ray eye camera sigma X/sigma Y
sigx = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaX_RBV', name='sigx')
sigy = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaY_RBV', name='sigy')



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
