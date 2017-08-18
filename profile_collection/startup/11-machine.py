from ophyd import (PVPositioner, Component as Cpt, EpicsSignal, EpicsSignalRO,
                   Signal)
from ophyd.utils import ReadOnlyError
import time as ttime


class UVDoneMOVN(Signal):
    # borrowed from srx's startup config
    """Signal for use as done signal for use in individual mode undulator
    motors This is a soft-signal that monitors several real PVs to sort out
    when the positioner is done moving.

    If the positioner looks like it has stopped (ex, moving is 0) but the
    readback is not close enough to the target, then re-actuate the motor.

    Parameters
    ----------
    parent : Device
         This comes from Cpt magic
    moving : str
        Name of the 'moving' signal on parent
    readback : str
        Name of the 'readback' signal on the parent
    actuate : str
        Name of the 'actuate' signal on the parent
    stop : str
        Name of the stop signal on the parent
    kwargs : ??
        All passed through to base Signal

    Attributes
    ----------
    target : float or None
        Where the positioner is going.  If ``None``, the callbacks
        short-circuit
    """
    def __init__(self, parent, moving, readback, actuate='actuate',
                 stop='stop_signal',
                 **kwargs):
        super().__init__(parent=parent, value=1, **kwargs)
        self._rbv = readback
        self._brake = moving
        self._act = actuate
        self._stp = stop
        self.target = None
        self._next_reactuate_time = 0

    def put(self, *arg, **kwargs):
        raise ReadOnlyError("You cannot tell an undulator motor it is done")

    def _put(self, *args, **kwargs):
        return super().put(*args, **kwargs)

    def _watcher(self, obj=None, **kwargs):
        '''The callback to install on readback and moving signals

        This callback watches if the position has gotten close enough to it's
        target _and_ has stopped moving, and then flips this signal to 1 (which
        in turn flips the Status object)
        '''
        target = self.target
        if target is None:
            return

        rbv = getattr(self.parent, self._rbv)
        moving = getattr(self.parent, self._brake)

        cur_value = rbv.get()
        not_moving = not moving.get()

        # come back and check this threshold value
        # this is 1 microns
        if not_moving and abs(target - cur_value) < 0.001:
            self._put(1)
            self._remove_cbs()
            return

        # if it is not moving, but we are not where we want to be,
        # poke it again
        if not_moving:
            cur_time = ttime.time()
            if cur_time > self._next_reactuate_time:
                actuate = getattr(self.parent, self._act)
                print('re actuated', self.parent.name)
                actuate.put(1)
                self._next_reactuate_time = cur_time + 1

    def _stop_watcher(self, *arg, **kwargs):
        '''Call back to be installed on the stop signal

        If this gets flipped, clear all of the other callbacks and tell the
        status object that it is done.

        TODO: mark status object as failed
        TODO: only trigger this on 0 -> 1 transposition
        '''
        # print('STOPPED')
        # set the target to None and remove all callbacks
        self.reset(None)
        # flip this signal to 1 to signal it is done
        self._put(1)
        # push stop again 'just to be safe'
        # this is paranoia related to the re-kicking the motor is the
        # other callback
        stop = getattr(self.parent, self._stp)
        stop.put(1)

    def reset(self, target):
        self.target = target
        self._put(0)
        self._remove_cbs()

    def _remove_cbs(self):
        rbv = getattr(self.parent, self._rbv)
        stop = getattr(self.parent, self._stp)
        moving = getattr(self.parent, self._brake)

        rbv.clear_sub(self._watcher)
        moving.clear_sub(self._watcher)
        stop.clear_sub(self._stop_watcher)


class UgapReadbackRO(EpicsSignalRO):
    '''Undulator gap readback scaling factor'''
    factor = 1e6
    precision = 4

    def subscribe(self, _fcn, event_type=None, **kwargs):
        # TODO this should be a DerivedSignal and have this taken care of there
        fcn = _fcn
        if event_type in (None, 'value'):
            def wrapped_fcn(*args, value=None, **kwargs):
                if value is not None:
                    value /= self.factor
                return _fcn(*args, value=value, **kwargs)

            fcn = wrapped_fcn

        return super().subscribe(fcn, event_type=event_type, **kwargs)

    def get(self, *args, **kwargs):
        return super().get(*args, **kwargs) / self.factor


class FixedPrecisionEpicsSignal(EpicsSignal):
    precision = 4


class UgapPositioner(PVPositioner):
    readback = Cpt(UgapReadbackRO, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.RBV')
    setpoint = Cpt(FixedPrecisionEpicsSignal,
                   'SR:C3-ID:G1{IVU20:1-Mtr:2}Inp:Pos')
    actuate = Cpt(EpicsSignal, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Sw:Go')
    actuate_value = 1

    stop_signal = Cpt(EpicsSignal, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.STOP')
    stop_value = 1

    # done = Cpt(EpicsSignalRO, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.DMOV')
    moving = Cpt(EpicsSignalRO, 'SR:C3-ID:G1{IVU20:1-Mtr:2}Pos.MOVN')
    done = Cpt(UVDoneMOVN, moving='moving', readback='readback')

    def move(self, position, *args, **kwargs):
        self.done.reset(position)
        self.done._next_reactuate_time = ttime.time() + 2
        ret = super().move(position, *args, **kwargs)
        self.moving.subscribe(self.done._watcher,
                              event_type=self.moving.SUB_VALUE)
        self.readback.subscribe(self.done._watcher,
                                event_type=self.readback.SUB_VALUE)

        self.stop_signal.subscribe(self.done._stop_watcher,
                                   event_type=self.stop_signal.SUB_VALUE,
                                   run=False)

        return ret

    def stop(self):
        self.done.reset(None)
        super().stop()


ugap = UgapPositioner(prefix='', settle_time=3., name='ugap')
ugap.read_attrs = ['setpoint', 'readback']

# Front End Slits (Primary Slits)
fe_tb = EpicsMotor('FE:C03A-OP{Slt:1-Ax:T}Mtr', name='fe_tb')
fe_bb = EpicsMotor('FE:C03A-OP{Slt:2-Ax:B}Mtr', name='fe_bb')
fe_ib = EpicsMotor('FE:C03A-OP{Slt:2-Ax:I}Mtr', name='fe_ib')
fe_ob = EpicsMotor('FE:C03A-OP{Slt:1-Ax:O}Mtr', name='fe_ob')

# Diamond XBPM motor
# xbpmb_x = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:8}Mtr', name='xbpmb_x')
# xbpmb_y = EpicsMotor('XF:03IDB-OP{Slt:SSA1-Ax:7}Mtr', name='xbpmb_y')

# Shutter operation
shutter_open = EpicsSignal('XF:03IDB-PPS{PSh}Cmd:Opn-Cmd', name = 'shutter_open')
shutter_close = EpicsSignal('XF:03IDB-PPS{PSh}Cmd:Cls-Cmd', name = 'shutter_close')
