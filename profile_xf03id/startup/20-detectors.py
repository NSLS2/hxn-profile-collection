from ophyd.controls import ProsilicaDetector, EpicsSignal, EpicsScaler

# CSX-1 Scalar

sclr_trig = EpicsSignal('XF:03IDB-ES{Sclr:1}.CNT', rw=True, name='sclr_trig')
ion0 = EpicsSignal('XF:03IDB-ES{Sclr:1}scaler_cts1.B', name='ion0')
ion1 = EpicsSignal('XF:03IDB-ES{Sclr:1}scaler_cts1.D', name='ion1')
