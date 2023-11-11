print(f"Loading {__file__!r} ...")

class HXN_FuncGen(Device):
    freq = Cpt(EpicsSignal, '{FG:1}OUTPUT1:FREQ:SP') #freq.set('#.##'){.wait()}
    freq_readout = Cpt(EpicsSignal, '{FG:1}OUTPUT1:FREQ') #freq.set('#.##'){.wait()}
    volt = Cpt(EpicsSignal, '{FG:1}OUTPUT1:VOLT:SP')
    volt_readout = Cpt(EpicsSignal, '{FG:1}OUTPUT1:VOLT')
    offset = Cpt(EpicsSignal, '{FG:1}OUTPUT1:VOLT:OFFSET:SP')
    sym = Cpt(EpicsSignal, '{FG:1}OUTPUT1:RAMP:SYMM:SP')
    func = Cpt(EpicsSignal, '{FG:1}OUTPUT1:FUNC:SP')
    burst_count = Cpt(EpicsSignal, '{FG:1}OUTPUT1:BURST_NCYCLES:SP')
    burst = Cpt(EpicsSignal, '{FG:1}OUTPUT1:BURST_STATUS:SP')
    output = Cpt(EpicsSignal, '{FG:1}OUTPUT1:STATUS:SP') #output.set('ON'){.wait()}
    output_readout = Cpt(EpicsSignal, '{FG:1}OUTPUT1:STATUS') #output.set('ON'){.wait()}
    #slt_hcen = Cpt(EpicsMotor, '{Slt:4-Ax:Top}Mtr')

    def on(self):
        yield from abs_set(self.output,"ON",wait=True)
        yield from bps.sleep(0.2)

    def off(self):
        yield from abs_set(self.output,"OFF",wait=True)
        yield from bps.sleep(0.2)

pt_fg = HXN_FuncGen('XF:03IDC-ES', name='pt_fg')
# pt_fg.output._put_complete = True
# pt_fg.burst._put_complete = True
