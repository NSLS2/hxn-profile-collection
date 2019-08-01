from ophyd import (PVPositioner, Component as Cpt, EpicsSignal, EpicsSignalRO,
                   Signal, EpicsMotor)
from ophyd.utils import ReadOnlyError
import time as ttime



class InsertionDevice(Device):
    gap = Cpt(EpicsMotor, '-Ax:Gap}-Mtr',
              kind='hinted', name='', settle_time=.1)
    brake = Cpt(EpicsSignal, '}BrakesDisengaged-Sts',
                write_pv='}BrakesDisengaged-SP',
                kind='omitted', add_prefix=('read_pv', 'write_pv', 'suffix'))


    def set(self, *args, **kwargs):
        from ophyd.utils import set_and_wait
        set_and_wait(self.brake, 1)
        return self.gap.set(*args, **kwargs)

    def stop(self, *, success=False):
        return self.gap.stop(success=success)

    @property
    def position(self):
        return self.gap.position

ugap = InsertionDevice('SR:C3-ID:G1{IVU20:1', name='ugap')
ugap.gap.user_readback.name = 'ugap_readback'
ugap.gap.user_setpoint.name = 'ugap_setpoint'

#ugap = UgapPositioner(prefix='', settle_time=3., name='ugap')


# Front End Slits (Primary Slits)
fe_tb = EpicsMotor('FE:C03A-OP{Slt:1-Ax:T}Mtr', name='fe_tb')
fe_bb = EpicsMotor('FE:C03A-OP{Slt:2-Ax:B}Mtr', name='fe_bb')
fe_ib = EpicsMotor('FE:C03A-OP{Slt:2-Ax:I}Mtr', name='fe_ib')
fe_ob = EpicsMotor('FE:C03A-OP{Slt:1-Ax:O}Mtr', name='fe_ob')

# Diamond XBPM motor
# xbpmb_x = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:8}Mtr', name='xbpmb_x')
# xbpmb_y = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:7}Mtr', name='xbpmb_y')

xbpmc_y = EpicsMotor('XF:03IDC-ES{BPM:7-Ax:Y}Mtr', name='xbpmc_y')

xbpmb_y = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:7}Mtr', name='xbpmb_y')
xbpmb_x = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:8}Mtr', name='xbpmb_x')


# Shutter operation
shutter_open = EpicsSignal('XF:03IDB-PPS{PSh}Cmd:Opn-Cmd', name='shutter_open')
shutter_close = EpicsSignal('XF:03IDB-PPS{PSh}Cmd:Cls-Cmd',
                            name='shutter_close')
