def test_d2scan():
    """
    Test ``d2scan`` scan (2D step scan) with ZP motors.
    """
    RE(d2scan([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,zpssy,-1,1,10,0.03))

