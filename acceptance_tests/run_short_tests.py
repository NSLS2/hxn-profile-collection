# Running the test from IPython:
# %run -i ~/.ipython/profile_collection/acceptance_tests/run_short_tests.py

def test_fly1d():
    """
    Test ``fly1d`` scan with ZP motor.
    """
    print("Running scan ..")
    uid, = RE(fly1d([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")


def test_fly2d():
    """
    Test ``fly2d`` scan with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(fly2d([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,zpssy,-1,1,10,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")

    print("Creating data file using PyXRF ('make_hdf') ...")
    from pyxrf.api import make_hdf
    dir_name = "/tmp/acceptance_tests"
    os.makedirs(dir_name, exist_ok=True)
    make_hdf(uid, wd=dir_name, create_each_det=True, file_overwrite_existing=True)

    scan_id = db[uid].start.scan_id
    fln = f"scan2D_{scan_id}.h5"
    path_fln = os.path.join(dir_name, fln)
    if not os.path.isfile(path_fln):
        assert False, f"PyXRF failed to create data file {path_fln!r} ..."
    print(f"Data file {path_fln!r} was successfully created")


def test_dscan():
    """
    Test ``dscan`` scan (1D step scan) with ZP motor.
    """
    print("Running scan ..")
    uid, = RE(dscan([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")



def test_d2scan():
    """
    Test ``d2scan`` scan (1D step scan) along two axes with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(d2scan([sclr1,zebra,merlin1,xspress3],10,zpssx,-1,1,zpssy,-1,1,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")



def test_mesh():
    """
    Test ``mesh`` scan (2D step scan) with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(mesh([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,zpssy,-1,1,10,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")

print("=================================================================================")
print(f"  The following test functions were loaded in the environment:")
print(f"      test_fly1d(), test_fly2d(), test_dscan(), test_d2scan(), test_mesh()")
print(f"  Run those functions manually to complete the test.")
print("=================================================================================")

# print("=====================================================================================")
# print("                              Testing 'test_fly1d' ...                               ")
# test_fly1d()

# print("=====================================================================================")
# print("                              Testing 'test_fly2d' ...                               ")
# test_fly2d()

# print("=====================================================================================")
# print("                              Testing 'test_dscan' ...                               ")
# test_dscan()

# print("=====================================================================================")
# print("                              Testing 'test_d2scan' ...                              ")
# test_d2scan()

# print("=====================================================================================")
# print("                              Testing 'test_mesh' ...                                ")
# test_mesh()

