from ophyd import (PVPositioner, Component as Cpt,
                   EpicsSignal, EpicsSignalRO)

# Undulator

class UgapPositioner(PVPositioner):
    setpoint = Cpt(EpicsSignal, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Inp:Pos')
    readback = Cpt(EpicsSignalRO, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.RBV')
    actuate = Cpt(EpicsSignal, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Sw:Go')

    actuate_value = 1
    stop_signal = Cpt(EpicsSignal, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.STOP')

    stop_value = 1
    done = Cpt(EpicsSignalRO, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.DMOV')


ugap = UgapPositioner(prefix='', settle_time=1.5, name='ugap')

# Front End Slits (Primary Slits)
fe_tb = EpicsMotor('FE:C03A-OP{Slt:1-Ax:T}Mtr.VAL', name='fe_tb')
fe_bb = EpicsMotor('FE:C03A-OP{Slt:2-Ax:B}Mtr.VAL', name='fe_bb')
fe_ib = EpicsMotor('FE:C03A-OP{Slt:2-Ax:I}Mtr.VAL', name='fe_ib')
fe_ob = EpicsMotor('FE:C03A-OP{Slt:1-Ax:O}Mtr.VAL', name='fe_ob')
