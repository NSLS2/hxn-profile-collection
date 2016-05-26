from hxntools.anc350 import dc_toggle


def dc_on(*, frequency=100):
    '''Enable duty cycling for all ANC350 controllers'''
    dc_toggle(True, freq=frequency)


def dc_off(*, frequency=1000):
    '''Disable duty cycling for all ANC350 controllers'''
    dc_toggle(False, freq=frequency)
