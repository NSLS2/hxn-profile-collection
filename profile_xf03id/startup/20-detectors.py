from ophyd.controls import EpicsSignal, EpicsScaler
from ophyd.controls import SignalDetector
import pandas as pd

from hxntools.detectors import (Xspress3Detector, TimepixDetector)
from hxntools.detectors import Xspress3ROI
from hxntools.detectors.zebra import HXNZebra

# Fly scan ROIs to display after the scan
FLY_XSPRESS3_ROI = [Xspress3ROI(chan=1, ev_low=9300, ev_high=9600),
                    Xspress3ROI(chan=2, ev_low=9300, ev_high=9600),
                    Xspress3ROI(chan=3, ev_low=9300, ev_high=9600),
                    ]

# Scaler 1 MCA channels numbers to record with fly scans
# (if above 8, be sure to modify n_scaler_mca below)
FLY_SCALER1_CHANS = [2, 3, 8]

# Flyscan results are shown using pandas. Maximum rows/columns to use when printing the table:
pd.options.display.width = 180
pd.options.display.max_rows = None

timepix1 = TimepixDetector('XF:03IDC-ES{Tpx:1}', files=['TIFF1:'], name='timepix1',
                           file_path='/data', ioc_file_path='/data')
timepix2 = TimepixDetector('XF:03IDC-ES{Tpx:2}', files=['TIFF1:'], name='timepix2',
                           file_path='/data', ioc_file_path='/data')
xspress3 = Xspress3Detector('XF:03IDC-ES{Xsp:1}:', cam='', files=['HDF5:'],
                            name='xspress3',
                            file_path='/xspress3_data/',
                            ioc_file_path='/xspress3_data/')

zebra = HXNZebra('XF:03IDC-ES{Zeb:1}:', name='zebra1')

# 3IDC RG:C4 VME scalers
sclr1 = EpicsScaler('XF:03IDC-ES{Sclr:1}', name='sclr1')
sclr2 = EpicsScaler('XF:03IDC-ES{Sclr:2}', name='sclr2')

sclr1_erase_start = EpicsSignal('XF:03IDC-ES{Sclr:1}EraseStart', rw=True, name='sclr1_erase_start')
sclr1_nuse_all = EpicsSignal('XF:03IDC-ES{Sclr:1}NuseAll', rw=True, name='sclr1_nuse_all')
sclr1_channel_advance = EpicsSignal('XF:03IDC-ES{Sclr:1}ChannelAdvance', rw=True, name='sclr1_channel_advance')
sclr1_input_mode = EpicsSignal('XF:03IDC-ES{Sclr:1}InputMode', rw=True, name='sclr1_input_mode')

n_scaler_mca = 8
sclr1_mca = [EpicsSignal('XF:03IDC-ES{Sclr:1}Mca:%d' % (i, ))
             for i in range(1, n_scaler_mca + 1)]

sclr1_trig = EpicsSignal('XF:03IDC-ES{Sclr:1}.CNT', rw=True, name='sclr1_trig')
det_sclr1 = SignalDetector(name='det_sclr1')
det_sclr1.add_acquire_signal(sclr1_trig)

sclr2_trig = EpicsSignal('XF:03IDC-ES{Sclr:2}.CNT', rw=True, name='sclr2_trig')
det_sclr2 = SignalDetector(name='det_sclr2')
det_sclr2.add_acquire_signal(sclr2_trig)

# Ion chamber
ion0 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.B', name='ion0')
ion1 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.C', name='ion1')
ionN = EpicsSignal('XF:03IDC-ES{Sclr:2}_calc4.VAL', name='ionN')
ion3 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.D', name='ion3')

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
Pt_ch1 = EpicsSignal('XF:03IDC-ES{Xsp:1}:C1_ROI1:Value_RBV', name='Pt_ch1')
Pt_ch2 = EpicsSignal('XF:03IDC-ES{Xsp:1}:C2_ROI1:Value_RBV', name='Pt_ch2')
Pt_ch3 = EpicsSignal('XF:03IDC-ES{Xsp:1}:C3_ROI1:Value_RBV', name='Pt_ch3')

# Fine sample stage readbacks
ssx_rbv = EpicsMotor('XF:03IDC-ES{Ppmac:1-ssx}Mtr.RBV', name='ssx_rbv')
ssy_rbv = EpicsMotor('XF:03IDC-ES{Ppmac:1-ssy}Mtr.RBV', name='ssy_rbv')
ssz_rbv = EpicsMotor('XF:03IDC-ES{Ppmac:1-ssz}Mtr.RBV', name='ssz_rbv')
