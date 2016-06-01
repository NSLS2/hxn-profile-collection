from ophyd import (EpicsSignal, EpicsSignalRO)
from ophyd import (Device, Component as Cpt)

import hxntools.handlers
from hxntools.detectors import (HxnTimepixDetector, HxnMerlinDetector,
                                BeamStatusDetector, HxnMercuryDetector)
from hxntools.detectors.zebra import HxnZebra


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
