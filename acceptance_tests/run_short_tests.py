# Running the test from IPython:
# %run -i ~/.ipython/profile_collection/acceptance_tests/run_short_tests.py

def test_fly1dpd():
    """
    Test ``fly1d`` scan with ZP motor.
    """
    print("Running scan ..")
    uid, = RE(fly1dpd(dets_fast,zpssx,-1,1,10,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")


def test_fly2dpd():
    """
    Test ``fly2d`` scan with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(fly2dpd(dets_fast,zpssx,-1,1,50,zpssy,-1,1,50,0.005))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")

def test_fly1d():
    """
    Test ``fly1d`` scan with ZP motor.
    """
    print("Running scan ..")
    uid, = RE(fly1d([sclr1,zebra,merlin1,xspress3,eiger1],zpssx,-1,1,10,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")


def test_fly2d():
    """
    Test ``fly2d`` scan with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(fly2d([sclr1,zebra,merlin1,xspress3,eiger1],zpssx,-1,1,10,zpssy,-1,1,10,0.03))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")

def test_dscan():
    """
    Test ``dscan`` scan (1D step scan) with ZP motor.
    """
    print("Running scan ..")
    uid, = RE(dscan([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,0.1))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")



def test_d2scan():
    """
    Test ``d2scan`` scan (1D step scan) along two axes with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(d2scan([sclr1,zebra,merlin1,xspress3],10,zpssx,-1,1,zpssy,-1,1,0.1))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")



def test_mesh():
    """
    Test ``mesh`` scan (2D step scan) with ZP motors.
    """
    print("Running scan ..")
    uid, = RE(mesh([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,5,zpssy,-1,1,5,0.1))
    print("Scan is completed")
    print("Filling the table ...")
    _ = db[uid].table(fill=True)
    print("Table is filled")

print("="*90)
print("="*90)
print(f"  The following test functions were loaded in the environment:")
print(f"      test_fly1d(), test_fly2d(), test_dscan(), test_d2scan(), test_mesh()")
print(f"  Run those functions manually to complete the test.")
print("="*90)
print("="*90)


print("="*90)
print("                              Testing 'test_fly1dpd' ...                               ")
print("="*90)
test_fly1dpd()

print("="*90)
print("                              Testing 'test_fly2dpd' ...                               ")
print("="*90)
test_fly2dpd()


print("="*90)
print("                              Testing 'test_fly1d' ...                               ")
print("="*90)
test_fly1d()

print("="*90)
print("                              Testing 'test_fly2d' ...                               ")
print("="*90)
test_fly2d()

print("="*90)
print("                              Testing 'test_dscan' ...                               ")
print("="*90)
test_dscan()

print("="*90)
print("                              Testing 'test_d2scan' ...                              ")
print("="*90)
test_d2scan()

print("="*90)
print("                              Testing 'test_mesh' ...                                ")
print("="*90)
test_mesh()
