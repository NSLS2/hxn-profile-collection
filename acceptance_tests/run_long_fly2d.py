# Running the test from IPython:
# %run -i ~/.ipython/profile_collection/acceptance_tests/run_long_fly2d.py

npts_fast_axis = 200
npts_slow_axis = 180
# npts_fast_axis = 20
# npts_slow_axis = 18


def test_long_fly2d():
    """
    Test ``fly2d`` scan with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(fly2d([sclr1,zebra,merlin1,xspress3],zpssx,-5,5,npts_fast_axis,zpssy,-4.5,4.5,npts_slow_axis,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")

    print("Creating data file using PyXRF ('make_hdf') ...")
    from pyxrf.api import make_hdf
    dir_name = "/tmp/acceptance_tests"
    os.makedirs(dir_name, exist_ok=True)
    make_hdf(uid, wd=dir_name, create_each_det=True, file_overwrite_existing=True)

    scan_id = db[uid].start["scan_id"]
    fln = f"scan2D_{scan_id}.h5"
    path_fln = os.path.join(dir_name, fln)
    if not os.path.isfile(path_fln):
        assert False, f"PyXRF failed to create data file {path_fln!r} ..."
    print(f"Data file {path_fln!r} was successfully created")


print("=================================================================================")
print(f"  The following test functions were loaded in the environment:")
print(f"      test_long_fly2d()")
print(f"  Run those functions manually to complete the test.")
print("=================================================================================")

# print("=====================================================================================")
# print("                         Testing 'test_long_fly2d' ...                               ")
# test_long_fly2d()

