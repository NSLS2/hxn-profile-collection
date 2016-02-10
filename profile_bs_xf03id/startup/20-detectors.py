from ophyd import (EpicsSignal, EpicsSignalRO, EpicsScaler)
import pandas as pd

from hxntools.detectors import (TimepixDetector,
                                HxnMerlinDetector, BeamStatusDetector)
from hxntools.detectors.zebra import HXNZebra


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
# merlin1 = HxnMerlinDetector('XF:31IDA-BI{Cam:Tbl}', name='merlin1')
merlin1.tiff1.read_attrs = []

zebra = HXNZebra('XF:03IDC-ES{Zeb:1}:', name='zebra')

# 3IDC RG:C4 VME scalers
sclr1 = EpicsScaler('XF:03IDC-ES{Sclr:1}', name='sclr1')
sclr2 = EpicsScaler('XF:03IDC-ES{Sclr:2}', name='sclr2')

sclr1_erase_start = EpicsSignal('XF:03IDC-ES{Sclr:1}EraseStart',
                                name='sclr1_erase_start')
sclr1_nuse_all = EpicsSignal('XF:03IDC-ES{Sclr:1}NuseAll',
                             name='sclr1_nuse_all')
sclr1_channel_advance = EpicsSignal('XF:03IDC-ES{Sclr:1}ChannelAdvance',
                                    name='sclr1_channel_advance')
sclr1_input_mode = EpicsSignal('XF:03IDC-ES{Sclr:1}InputMode',
                               name='sclr1_input_mode')

n_scaler_mca = 8
sclr1_mca = [EpicsSignalRO('XF:03IDC-ES{Sclr:1}Mca:%d' % (i, ))
             for i in range(1, n_scaler_mca + 1)]

sclr1_trig = EpicsSignal('XF:03IDC-ES{Sclr:1}.CNT', name='sclr1_trig')
sclr2_trig = EpicsSignal('XF:03IDC-ES{Sclr:2}.CNT', name='sclr2_trig')

# ugap scan trigger
ugap_trig = EpicsSignal('SR:C3-ID:G1{IVU20:1-Mtr:2}Sw:Go', name='ugap_trig')


# Ion chamber
sclr2_ch2 = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_cts1.B', name='sclr2_ch2')
sclr2_ch3 = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_cts1.C', name='sclr2_ch3')
sclr2_ch4_calc = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_calc4.VAL',
                               name='sclr2_ch4_calc')
sclr2_ch4 = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_cts1.D', name='sclr2_ch4')

sclr1_ch2 = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_cts1.B', name='sclr1_ch2')
sclr1_ch3 = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_cts1.C', name='sclr1_ch3')
sclr1_ch4_calc = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_calc4.VAL',
                               name='sclr1_ch4_calc')
sclr1_ch4 = EpicsSignalRO('XF:03IDC-ES{Sclr:2}_cts1.D', name='sclr1_ch4')

t_base = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:D}C:T-I', name='t_base')
t_sample = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:C}C:T-I', name='t_sample')
t_vlens = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:B}C:T-I', name='t_vlens')
t_hlens = EpicsSignalRO('XF:03IDC-ES{LS:2-Ch:A}C:T-I', name='t_hlens')

# X-ray eye camera sigma X/sigma Y
sigx = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaX_RBV', name='sigx')
sigy = EpicsSignalRO('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaY_RBV', name='sigy')

# Interferometers
int_sx = EpicsSignalRO('XF:03IDC-ES{FPS:1-Chan0}Pos-I', name='int_sx')
int_sy = EpicsSignalRO('XF:03IDC-ES{FPS:1-Chan1}Pos-I', name='int_sy')
int_sz = EpicsSignalRO('XF:03IDC-ES{FPS:1-Chan2}Pos-I', name='int_sz')

int_hx = EpicsSignalRO('XF:03IDC-ES{FPS:2-Chan0}Pos-I', name='int_hx')
int_vy = EpicsSignalRO('XF:03IDC-ES{FPS:2-Chan1}Pos-I', name='int_vy')
int_hy = EpicsSignal('XF:03IDC-ES{FPS:2-Chan2}Pos-I', name='int_hy')

int_hz = EpicsSignal('XF:03IDC-ES{FPS:3-Chan0}Pos-I', name='int_hz')
int_vx = EpicsSignal('XF:03IDC-ES{FPS:3-Chan1}Pos-I', name='int_vx')
int_vz = EpicsSignal('XF:03IDC-ES{FPS:3-Chan2}Pos-I', name='int_vz')

int_zpssx = EpicsSignal('XF:03IDC-ES{FPS:4-Chan0}Pos-I', name='int_zpssx')
int_zpssy = EpicsSignal('XF:03IDC-ES{FPS:4-Chan1}Pos-I', name='int_zpssy')
int_zpssz = EpicsSignal('XF:03IDC-ES{FPS:4-Chan2}Pos-I', name='int_zpssz')

int_zpx1 = EpicsSignal('XF:03IDC-ES{FPS:5-Chan0}Pos-I', name='int_zpx1')
int_zpx2 = EpicsSignal('XF:03IDC-ES{FPS:5-Chan1}Pos-I', name='int_zpx2')
int_zpy1 = EpicsSignal('XF:03IDC-ES{FPS:5-Chan2}Pos-I', name='int_zpy1')

int_zpy2 = EpicsSignal('XF:03IDC-ES{FPS:6-Chan0}Pos-I', name='int_zpy2')
int_zpz = EpicsSignal('XF:03IDC-ES{FPS:6-Chan1}Pos-I', name='int_zpz')
int_zpspare1 = EpicsSignal('XF:03IDC-ES{FPS:6-Chan2}Pos-I',
                           name='int_zpspare1')

# Front-end Xray BPMs and local bumps
xbpm_x = EpicsSignalRO('SR:C03-BI{XBPM:1}Pos:X-I', name='xbpm_x')
xbpm_y = EpicsSignalRO('SR:C03-BI{XBPM:1}Pos:Y-I', name='xbpm_y')

angle_x = EpicsSignalRO('SR:C31-{AI}Aie3:Angle-x-Cal', name='angle_x')
angle_y = EpicsSignalRO('SR:C31-{AI}Aie3:Angle-y-Cal', name='angle_y')

# Diamond Quad BPMs in C hutch
quad_x = EpicsSignalRO('SR:C12-BI{XBPM:1}Pos:X-I', name='quad_x')
quad_y = EpicsSignalRO('SR:C12-BI{XBPM:1}Pos:Y-I', name='quad_y')


# Slit 1 BPM (drain current from I400)
slit1_top = EpicsSignalRO('XF:03IDA-BI{Slt:1}I:Raw1-I', name='slit1_top')
slit1_bottom = EpicsSignalRO('XF:03IDA-BI{Slt:1}I:Raw2-I', name='slit1_bottom')
slit1_right = EpicsSignalRO('XF:03IDA-BI{Slt:1}I:Raw3-I', name='slit1_right')
slit1_left = EpicsSignalRO('XF:03IDA-BI{Slt:1}I:Raw4-I', name='slit1_left')

# calculated position (from CALC record)
slit1_xpos = EpicsSignalRO('XF:03IDA-BI{Slt:1}PosX-I', name='slit1_xpos')
slit1_ypos = EpicsSignalRO('XF:03IDA-BI{Slt:1}PosY-I', name='slit1_ypos')

# Mono motors
dcm_th = EpicsSignal('XF:03IDA-OP{Mon:1-Ax:Bragg}Mtr.RBV', name='dcm_th')
dcm_p = EpicsSignal('XF:03IDA-OP{Mon:1-Ax:P}Mtr.RBV', name='dcm_p')

#tpx1_roi = EpicsSignal('XF:03IDC-ES{Tpx:1}Stats1:Total_RBV', name='tpx1_roi')


sr_shutter_status = EpicsSignalRO('SR-EPS{PLC:1}Sts:MstrSh-Sts', name='sr_shutter_status')
sr_beam_current = EpicsSignalRO('SR:C03-BI{DCCT:1}I:Real-I', name='sr_beam_current')

det_beamstatus = BeamStatusDetector(min_current=100.0)
