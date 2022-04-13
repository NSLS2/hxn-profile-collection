def test_mesh():
    """
    Test ``mesh`` scan (2D step scan) with ZP motors.
    """
    RE(mesh([sclr1,zebra,merlin1,xspress3],zpssx,-1,1,10,zpssy,-1,1,10,0.03))

