#!/usr/bin/env python
"""
Script to control the fan speed of an NVidia gpu using a custom fan speed/temperature curve.
Copyright (C) 2012  Luke Frisken

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see http://www.gnu.org/licenses.
"""

from subprocess import *
import time
import os
import sys
import signal
import threading

"""
SETUP:
In order to get this program to work, my xorg.conf was generated
and saved using the x server display configuration tab in nvidia-settings.
These lines:

	Section "Device"
		Identifier     "Device0"
		Driver         "nvidia"
		VendorName     "NVIDIA Corporation"
		BoardName      "GeForce GTX 570"
	EndSection

in /etc/X11/xorg.conf
were then edited to look like this:

	Section "Device"
		Identifier     "Device0"
		Driver         "nvidia"
		VendorName     "NVIDIA Corporation"
		BoardName      "GeForce GTX 570"
		Option "Coolbits" "4"
	EndSection

then the compter is restarted.


RUNNING:
There are a few bugs which I haven't been able to solve, so observe carefully here:
 - try to run from an already open terminal.
 - try to end program using ctrl-C
 - ensure that when finished, the box in nvidia-settings for the manual fan control is unticked!!!!!
"""
class Vector():
	def __init__(self, x, y):
		self.x = x
		self.y = y

	def gradient(self):
		return float(y)/float(x)

	def __str__(self):
		return "Vector: {0}, {1}".format(self.x, self.y)


class Curve():
	def __init__(self, curve_point_array):
		#convert curve to use vectors for points
		"""self.vector_curve = []
		for v in curve:
			self.vector_curve.append(Vector(v[0], v[1]))"""

		self.cpa = curve_point_array

	def evaluate(self, x):
		point_i = 0
		while(point_i < len(self.cpa) - 1):

			if(self.cpa[point_i][0] <= x and self.cpa[point_i + 1][0] > x):
				point_1 = self.cpa[point_i]
				point_2 = self.cpa[point_i + 1]
				delta_x = point_2[0] - point_1[0]
				delta_y = point_2[1] - point_1[1]
				#print("delta_x: {0}, delta_y: {1}".format(delta_x, delta_y))
				gradient = float(delta_y)/float(delta_x)
				#print("x: {0} point_1: [{1}, {2}] gradient: {3}".format(x, point_1[0], point_1[1], gradient))
				x_bit = x - point_1[0]
				y_bit = int(float(x_bit) * gradient)
				y = point_1[1] + y_bit
				return y

			point_i += 1




def clearScreen():
	os.system("clear")

class StoppableThread(threading.Thread):
    """Thread class with a stop() method. The thread itself has to check
    regularly for the stopped() condition."""

    def __init__(self):
        super(StoppableThread, self).__init__()
        self._stop = threading.Event()

    def stop(self):
        self._stop.set()

    def stopped(self):
        return self._stop.isSet()

class FanController(StoppableThread):
	"""Room here for arguments to implement multigpu fan controll"""
	def __init__(self):
		super(FanController, self).__init__()
		self.daemon = True
		
	def run(self):
		self.customFanSpeed()
		
	def stop(self):
		super(FanController, self).stop()
		time.sleep(1.5)
		self.resetFan()
	
	def resetFan(self):
		print("\nReset to Auto Fan")
		process = Popen("nvidia-settings -a [gpu:0]/GPUFanControlState=0", shell=True, stdin=PIPE, stdout=PIPE)
	
	def getTemp(self):
		process = Popen("nvidia-settings -q gpucoretemp", shell=True, stdin=PIPE, stdout=PIPE)
		line_array = process.stdout.readlines()
		tmp_line = line_array[1]
		#grab number from end of line
		return int(tmp_line[-4:-2])

	def setFanSpeed(self, speed):
        #watch -n0 aticonfig --adapter=0 --od-gettemperature
        #that's the line to get ati fan speed
		process = Popen("nvidia-settings -a [gpu:0]/GPUFanControlState=1 -a [fan:0]/GPUCurrentFanSpeed={0}".format(speed), shell=True, stdin=PIPE, stdout=PIPE)
		return
	
	def customFanSpeed(self):
		"""custom fan speed curve example:
			[[temp, speed],
			 [temp, speed],
			 [temp, speed]]
			always start at low temp/speed and head towards high temp/speed
			ensure that first point is always lower temp than possible
			ensure that gradient is always positive and less than infinity
		"""
		curve_point_array = [[10, 5],
							 [30, 20],
							 [40, 30],
							 [50, 50],
							 [60, 70],
							 [66, 80],
							 [70, 99],
							 [100, 100]]

		curve = Curve(curve_point_array)


		while(not self.stopped()):
			current_temp = self.getTemp()
			new_fan_speed = curve.evaluate(current_temp)
			clearScreen()
			print("CurrTemp: {0} FanSpd: {1}".format(current_temp,new_fan_speed))
			self.setFanSpeed(new_fan_speed)

			time.sleep(1.0)
		
		#finished and ready to exit
		return


class MainThread(threading.Thread):
	def __init__(self):
		super(MainThread, self).__init__()
		self.fan_controller = FanController()
		signal.signal(signal.SIGINT, self.exit_signal_handler)
		
	def exit_signal_handler(self, signal, frame):
		self.fan_controller.stop()
        #sys.exit(0)
        
	def start(self):
		super(MainThread, self).start()
		signal.pause()
		
	def run(self):
		self.fan_controller.start()


t = MainThread()
t.start()



