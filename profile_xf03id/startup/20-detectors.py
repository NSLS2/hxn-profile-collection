from ophyd.controls import ProsilicaDetector, EpicsSignal, EpicsScaler
from ophyd.controls import SignalDetector

# 3IDC RG:C4 VME scalers
sclr1 = EpicsScaler('XF:03IDC-ES{Sclr:1}', name='sclr1')
sclr2 = EpicsScaler('XF:03IDC-ES{Sclr:2}', name='sclr2')

sclr_trig = EpicsSignal('XF:03IDC-ES{Sclr:2}.CNT', rw=True, name='sclr_trig')
det_sclr2 = SignalDetector()
det_sclr2.add_acquire_signal(sclr_trig)

# Ion chamber
ion0 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.B', name='ion0')
ion1 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts1.C', name='ion1')
ionN = EpicsSignal('XF:03IDC-ES{Sclr:2}_calc4.VAL', name='ionN')

ion3 = EpicsSignal('XF:03IDC-ES{Sclr:2}_cts2.D', name='ion3')

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
xrf_erase = EpicsSignal('XF:03IDC-ES{Xsp:1}:ERASE', rw=True, name='xrf_erase')
xrf_trig = EpicsSignal('XF:03IDC-ES{Xsp:1}:Acquire', rw=True, name='xrf_trig')

det_xrf_erase = SignalDetector()
det_xrf_erase.add_acquire_signal(xrf_erase)
det_xrf_trig = SignalDetector()
det_xrf_trig.add_acquire_signal(xrf_trig)
