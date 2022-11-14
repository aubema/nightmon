#!/usr/bin/python3
# move the camera assemply using half step mode on a bipolar stepper
#
# usage: move_cams.py steps slow_level
# steps = number of steps
# slow_level 1=fastest n=max_speed/n
import RPi.GPIO as GPIO
import time
import sys

# Variables
reverse=1
steps = int(sys.argv[1])
if steps<0:
   steps=-1*steps
   reverse=0
   # set a maximum number of steps
   if steps>1500:
   	steps=1500
   	print("Too much steps")
steps=steps-1
delay = float(sys.argv[2]) * 0.0075

GPIO.setmode(GPIO.BCM)
GPIO.setwarnings(False)

# set the sensor pin
limit_gpio=5
GPIO.setup(limit_gpio, GPIO.IN)

# Enable pins for IN1-4 to control step sequence

coil_A_1_pin = 16
coil_A_2_pin = 12
coil_B_1_pin = 20
coil_B_2_pin = 21

# Set pin states

GPIO.setup(coil_A_1_pin, GPIO.OUT)
GPIO.setup(coil_A_2_pin, GPIO.OUT)
GPIO.setup(coil_B_1_pin, GPIO.OUT)
GPIO.setup(coil_B_2_pin, GPIO.OUT)

# Function for step sequence
def setStep(w1, w2, w3, w4):
  time.sleep(delay)
  GPIO.output(coil_A_1_pin, w1)
  GPIO.output(coil_A_2_pin, w2)
  GPIO.output(coil_B_1_pin, w3)
  GPIO.output(coil_B_2_pin, w4)

# loop through step sequence based on number of steps
j=0
if reverse==0:
	for i in range(0, steps):
		j=j+1
		if j==1:
			setStep(1,0,1,0)
		if j==2:
			setStep(1,0,0,0)
		if j==3:
			setStep(1,0,0,1)
		if j==4:
			setStep(0,0,0,1)
		if j==5:
			setStep(0,1,0,1)
		if j==6:
			setStep(0,1,0,0)
		if j==7:
			setStep(0,1,1,0)
		if j==8:
			setStep(0,0,1,0)
			j=0


# Reverse previous step sequence to reverse motor direction
else:
	for i in range(0, steps):
		if GPIO.input(limit_gpio)==0:
			print("Limit switch activated")
			setStep(0,0,0,0)
			break
		j=j+1
		if j==1:
			setStep(0,0,1,0)
		if j==2:
			setStep(0,1,1,0)
		if j==3:
			setStep(0,1,0,0)
		if j==4:
			setStep(0,1,0,1)
		if j==5:
			setStep(0,0,0,1)
		if j==6:
			setStep(1,0,0,1)
		if j==7:
			setStep(1,0,0,0)
		if j==8:
			setStep(1,0,1,0)
			j=0
setStep(0,0,0,0)
