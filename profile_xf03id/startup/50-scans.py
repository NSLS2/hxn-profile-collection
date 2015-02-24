from ophyd.userapi.scan_api import Scan, AScan, DScan, Count

scan = Scan()
ascan = AScan()
ascan.default_triggers = [sclr_trig]
ascan.default_detectors = [ion0, ion1]
dscan = DScan()

# Use ct as a count which is a single scan.

ct = Count()
