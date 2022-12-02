import getopt
import os
import sys

import RPi.GPIO as gp

gp.setwarnings(False)
gp.setmode(gp.BOARD)

gp.setup(7, gp.OUT)
gp.setup(11, gp.OUT)
gp.setup(12, gp.OUT)


def input(argv):
    itime = 1000
    gain = 8
    try:
        opts, args = getopt.getopt(argv, "h:t:g:", ["help=", "itime=", "gain="])
    except getopt.GetoptError:
        print("captureA.py -t <itime> -g <gain>")
        sys.exit(2)
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            print("captureA.py -t <itime> -g <gain>")
            sys.exit()
        elif opt in ("-t", "--itime"):
            itime = arg
        elif opt in ("-g", "--gain"):
            gain = arg
        return itime, gain


itime, gain = input(sys.argv[1:])


def main():
    print("Start testing the camera A")
    i2c = "i2cset -y 1 0x70 0x00 0x04"
    os.system(i2c)
    gp.output(7, False)
    gp.output(11, False)
    gp.output(12, True)
    print("selected integration V : ", itime)
    capture(1)


def capture(cam):
    # cmd = "libcamera-hello -t 0"
    cmd = (
        "libcamera-still --analoggain gain --shutter itime --denoise off --rawfull "
        "--raw --awbgains 1,1 --nopreview -o capture_%d.jpg" % cam
    )
    os.system(cmd)


if __name__ == "__main__":
    main()

    gp.output(7, False)
    gp.output(11, False)
    gp.output(12, True)
