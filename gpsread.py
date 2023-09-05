from datetime import datetime
import time
import os
import threading
import serial
import pynmea2


# Get gps position


def getPositionData():
    global end
    global lat
    global lon
    global alt
    global nbSats
    global times

    SERIAL_PORT = "/dev/ttyACM0"
    try:
        gps = serial.Serial(SERIAL_PORT, baudrate=9600, timeout=0.5)
    except serial.SerialException:
        print("No GPS module found...")

    while end == 0:
        try:
            data = gps.readline()
            message = str(data[0:6])
            message = message[2:8]
            message_all = str(data)
            message_all = message_all[2:-1]

            if (message == "$GPGGA"):
                # GPGGA = Global Positioning System Fix Data
                # Reading the GPS fix data is an alternative approach that also works
                today = str(datetime.now())
                hour = int(today[11:13])
                minute = int(today[14:16])
                second = float(today[17:])
                times[0] = times[1]
                times[1] = hour*3600 + minute*60 + second
                data = str(data)
                data = data[2:-5]

                parts = pynmea2.parse(data)
                print(parts)

                if int(parts.gps_qual) == 0:
                    # Equal to 0 = No gps fix...
                    print("No gps fix")
                    lat[0] = 0
                    lon[0] = 0
                    alt[0] = 0
                    lat[1] = 0
                    lon[1] = 0
                    alt[1] = 0
                    nbSats = 0
                else:
                    # Get the position data that was transmitted with the GPGGA message
                    lat[0] = lat[1]
                    lon[0] = lon[1]
                    alt[0] = alt[1]

                    lat[1] = float("{:.6f}".format(parts.latitude))
                    lon[1] = float("{:.6f}".format(parts.longitude))
                    alt[1] = float("{:.6f}".format(parts.altitude))
                    nbSats = int(parts.num_sats)

            else:
                # Handle other NMEA messages and unsupported strings
                pass
        except KeyboardInterrupt:
            gps.close()
            print("Application closed!")
        except:
            lat[0] = 0
            lon[0] = 0
            alt[0] = 0
            lat[1] = 0
            lon[1] = 0
            alt[1] = 0
            nbSats = 0

            time.sleep(0.9)
            try:
                if SERIAL_PORT == "/dev/ttyACM0":
                    SERIAL_PORT = "/dev/ttyACM1"
                elif SERIAL_PORT == "/dev/ttyACM1":
                    SERIAL_PORT = "/dev/ttyACM0"

                gps = serial.Serial(SERIAL_PORT, baudrate=9600, timeout=0.5)

            except serial.SerialException:
                print("No GPS module found...")


# initialisation

today = str(datetime.now())
year = today[0:4]
month = today[5:7]
day = today[8:10]
lat = [0, 0]
lon = [0, 0]
alt = [0, 0]
nbSats = 0
times = [0, 0]

# Gps thread initialisation
tGps = threading.Thread(target=getPositionData, name="Gps thread")
tGps.start()

print("Coords=",lat[1],lon[1])
