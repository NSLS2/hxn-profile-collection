
from fly1d_zp import test_fly1d
from fly2d_zp import test_fly2d

from dscan_zp import test_dscan
from d2scan_zp import test_d2scan

from mesh_zp import test_mesh

print("=====================================================================================")
print("                              Testing 'test_fly1d' ...                               ")
test_fly1d()

print("=====================================================================================")
print("                              Testing 'test_fly2d' ...                               ")
test_fly2d()

print("=====================================================================================")
print("                              Testing 'test_dscan' ...                               ")
test_dscan()

print("=====================================================================================")
print("                              Testing 'test_d2scan' ...                              ")
test_d2scan()

print("=====================================================================================")
print("                              Testing 'test_mesh' ...                                ")
test_mesh()

