def test_dscan():
    """
    Test ``dscan`` scan (1D step scan) with ZP motor.
    """
    RE(dscan([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,0.03))

