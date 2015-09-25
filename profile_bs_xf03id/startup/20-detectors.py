from ophyd.controls import EpicsSignal, EpicsScaler
from ophyd.controls import SignalDetector
import pandas as pd

from hxntools.detectors import (TimepixDetector,
                                MerlinDetector, BeamStatusDetector)
from hxntools.detectors.zebra import HXNZebra


# Scaler 1 MCA channels numbers to record with fly scans
# (if above 8, be sure to modify n_scaler_mca below)
FLY_SCALER1_CHANS = [2, 3, 8]

# Flyscan results are shown using pandas. Maximum rows/columns to use when
# printing the table:
pd.options.display.width = 180
pd.options.display.max_rows = None
pd.options.display.max_columns = 10

timepix1 = TimepixDetector('XF:03IDC-ES{Tpx:1}', files=['TIFF1:'],
                           name='timepix1',
                           file_path='/data', ioc_file_path='/data')
timepix2 = TimepixDetector('XF:03IDC-ES{Tpx:2}', files=['TIFF1:'],
                           name='timepix2',
                           file_path='/data', ioc_file_path='/data')
merlin1 = MerlinDetector('XF:03IDC-ES{Merlin:1}', files=['TIFF1:'],
                         name='merlin1',
                         file_path='/data', ioc_file_path='/data')

zebra = HXNZebra('XF:03IDC-ES{Zeb:1}:', name='zebra')

# 3IDC RG:C4 VME scalers
sclr1 = EpicsScaler('XF:03IDC-ES{Sclr:1}', name='sclr1')
sclr2 = EpicsScaler('XF:03IDC-ES{Sclr:2}', name='sclr2')

sclr1_erase_start = EpicsSignal('XF:03IDC-ES{Sclr:1}EraseStart', rw=True,
                                name='sclr1_erase_start')
sclr1_nuse_all = EpicsSignal('XF:03IDC-ES{Sclr:1}NuseAll', rw=True,
                             name='sclr1_nuse_all')
sclr1_channel_advance = EpicsSignal('XF:03IDC-ES{Sclr:1}ChannelAdvance',
                                    rw=True, name='sclr1_channel_advance')
sclr1_input_mode = EpicsSignal('XF:03IDC-ES{Sclr:1}InputMode', rw=True,
                               name='sclr1_input_mode')

n_scaler_mca = 8
sclr1_mca = [EpicsSignal('XF:03IDC-ES{Sclr:1}Mca:%d' % (i, ))
             for i in range(1, n_scaler_mca + 1)]

sclr1_trig = EpicsSignal('XF:03IDC-ES{Sclr:1}.CNT', rw=True, name='sclr1_trig')
sclr2_trig = EpicsSignal('XF:03IDC-ES{Sclr:2}.CNT', rw=True, name='sclr2_trig')

# Ion chamber
sclr2_ch2 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.B', name='sclr2_ch2')
sclr2_ch3 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.C', name='sclr2_ch3')
sclr2_ch4_calc = EpicsSignal('XF:03IDC-ES{Sclr:2}_calc4.VAL',
                             name='sclr2_ch4_calc')
sclr2_ch4 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.D', name='sclr2_ch4')

sclr1_ch2 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.B', name='sclr1_ch2')
sclr1_ch3 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.C', name='sclr1_ch3')
sclr1_ch4_calc = EpicsSignal('XF:03IDC-ES{Sclr:2}_calc4.VAL',
                             name='sclr1_ch4_calc')
sclr1_ch4 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.D', name='sclr1_ch4')

t_base = EpicsSignal('XF:03IDC-ES{LS:2-Ch:D}C:T-I', name='t_base')
t_sample = EpicsSignal('XF:03IDC-ES{LS:2-Ch:C}C:T-I', name='t_sample')
t_vlens = EpicsSignal('XF:03IDC-ES{LS:2-Ch:B}C:T-I', name='t_vlens')
t_hlens = EpicsSignal('XF:03IDC-ES{LS:2-Ch:A}C:T-I', name='t_hlens')

# X-ray eye camera sigma X/sigma Y
sigx = EpicsSignal('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaX_RBV', name='sigx')
sigy = EpicsSignal('XF:03IDB-BI{Xeye-CAM:1}Stats1:SigmaY_RBV', name='sigy')

# Interferometers
int_sx = EpicsSignal('XF:03IDC-ES{FPS:1-Chan0}Pos-I', name='int_sx')
int_sy = EpicsSignal('XF:03IDC-ES{FPS:1-Chan1}Pos-I', name='int_sy')
int_sz = EpicsSignal('XF:03IDC-ES{FPS:1-Chan2}Pos-I', name='int_sz')
int_hx = EpicsSignal('XF:03IDC-ES{FPS:2-Chan0}Pos-I', name='int_hx')
int_vy = EpicsSignal('XF:03IDC-ES{FPS:2-Chan1}Pos-I', name='int_vy')

# Unused interferometer signals
# int_2ch2 = EpicsSignal('XF:03IDC-ES{FPS:2-Chan2}Pos-I', name='int_2ch2')
# int_3ch0 = EpicsSignal('XF:03IDC-ES{FPS:3-Chan0}Pos-I', name='int_3ch0')
# int_3ch1 = EpicsSignal('XF:03IDC-ES{FPS:3-Chan1}Pos-I', name='int_3ch1')
# int_3ch2 = EpicsSignal('XF:03IDC-ES{FPS:3-Chan2}Pos-I', name='int_3ch2')
# int_4ch0 = EpicsSignal('XF:03IDC-ES{FPS:4-Chan0}Pos-I', name='int_4ch0')
# int_4ch1 = EpicsSignal('XF:03IDC-ES{FPS:4-Chan1}Pos-I', name='int_4ch1')
# int_4ch2 = EpicsSignal('XF:03IDC-ES{FPS:4-Chan2}Pos-I', name='int_4ch2')
# int_5ch0 = EpicsSignal('XF:03IDC-ES{FPS:5-Chan0}Pos-I', name='int_5ch0')
# int_5ch1 = EpicsSignal('XF:03IDC-ES{FPS:5-Chan1}Pos-I', name='int_5ch1')
# int_5ch2 = EpicsSignal('XF:03IDC-ES{FPS:5-Chan2}Pos-I', name='int_5ch2')

# Fluorescence detectors
det1_Al = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI10:Value_RBV', name='Det1_Al')
det2_Al = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI10:Value_RBV', name='Det2_Al')
det3_Al = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI10:Value_RBV', name='Det3_Al')

det1_Si = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI16:Value_RBV', name='Det1_Si')
det2_Si = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI16:Value_RBV', name='Det2_Si')
det3_Si = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI16:Value_RBV', name='Det3_Si')

det1_Gd = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI3:Value_RBV', name='Det1_Gd')
det2_Gd = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI3:Value_RBV', name='Det2_Gd')
det3_Gd = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI3:Value_RBV', name='Det3_Gd')

det1_YSZ = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI4:Value_RBV', name='Det1_YSZ')
det2_YSZ = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI4:Value_RBV', name='Det2_YSZ')
det3_YSZ = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI4:Value_RBV', name='Det3_YSZ')

det1_Ca = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI5:Value_RBV', name='Det1_Ca')
det2_Ca = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI5:Value_RBV', name='Det2_Ca')
det3_Ca = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI5:Value_RBV', name='Det3_Ca')

det1_Ti = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI6:Value_RBV', name='Det1_Ti')
det2_Ti = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI6:Value_RBV', name='Det2_Ti')
det3_Ti = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI6:Value_RBV', name='Det3_Ti')

det1_V = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI7:Value_RBV', name='Det1_V')
det2_V = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI7:Value_RBV', name='Det2_V')
det3_V = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI7:Value_RBV', name='Det3_V')

det1_Cr = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI8:Value_RBV', name='Det1_Cr')
det2_Cr = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI8:Value_RBV', name='Det2_Cr')
det3_Cr = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI8:Value_RBV', name='Det3_Cr')

det1_Mn = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI9:Value_RBV', name='Det1_Mn')
det2_Mn = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI9:Value_RBV', name='Det2_Mn')
det3_Mn = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI9:Value_RBV', name='Det3_Mn')

det1_Fe = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI2:Value_RBV', name='Det1_Fe')
det2_Fe = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI2:Value_RBV', name='Det2_Fe')
det3_Fe = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI2:Value_RBV', name='Det3_Fe')

det1_Co = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI11:Value_RBV', name='Det1_Co')
det2_Co = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI11:Value_RBV', name='Det2_Co')
det3_Co = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI11:Value_RBV', name='Det3_Co')

det1_Ni = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI12:Value_RBV', name='Det1_Ni')
det2_Ni = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI12:Value_RBV', name='Det2_Ni')
det3_Ni = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI12:Value_RBV', name='Det3_Ni')

det1_Cu = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI3:Value_RBV', name='Det1_Cu')
det2_Cu = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI3:Value_RBV', name='Det2_Cu')
det3_Cu = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI3:Value_RBV', name='Det3_Cu')

det1_Ce = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI4:Value_RBV', name='Det1_Ce')
det2_Ce = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI4:Value_RBV', name='Det2_Ce')
det3_Ce = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI4:Value_RBV', name='Det3_Ce')

det1_Au = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI15:Value_RBV', name='Det1_Au')
det2_Au = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI15:Value_RBV', name='Det2_Au')
det3_Au = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI15:Value_RBV', name='Det3_Au')

det1_Pt = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI1:Value_RBV', name='Det1_Pt')
det2_Pt = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI1:Value_RBV', name='Det2_Pt')
det3_Pt = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI1:Value_RBV', name='Det3_Pt')

# Fine sample stage readbacks
ssx_rbv = EpicsSignal('XF:03IDC-ES{Ppmac:1-ssx}Mtr.RBV', name='ssx_rbv')
ssy_rbv = EpicsSignal('XF:03IDC-ES{Ppmac:1-ssy}Mtr.RBV', name='ssy_rbv')
ssz_rbv = EpicsSignal('XF:03IDC-ES{Ppmac:1-ssz}Mtr.RBV', name='ssz_rbv')

sr_shutter_status = EpicsSignal('SR-EPS{PLC:1}Sts:MstrSh-Sts', rw=False)
sr_beam_current = EpicsSignal('SR:C03-BI{DCCT:1}I:Real-I', rw=False)

det_beamstatus = BeamStatusDetector(min_current=100.0)
