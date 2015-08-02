#test script for grabbing images for the optical pointing of the antenna
# Revision 1 : 30/3/2015
# Charles Copley


import cv2 #for ubuntu 12.04 install see: http://karytech.blogspot.com/2012/05/opencv-24-on-ubuntu-1204.html
import cv
import numpy as np
import time
import datetime
import sys, optparse
import pdb
import matplotlib.pylab as plt
import os

#which device (i.e. /dev/video0) should we use?
CAMERA_DEVICE_INDEX=0


cam0t=cv2.VideoCapture(CAMERA_DEVICE_INDEX)
#set up v4l2 to the correct input on the frame grabber
#os.system("v4l2-ctl -d /dev/video0 -i 0 --set-fmt-video=width=768,height=480")
#time.sleep(1)
os.system("v4l2-ctl -d /dev/video0 -i 1 --set-fmt-video=width=768,height=480")
#time.sleep(1)
#and set up the correct width
width =768
height = 480
cam0t.set(3,width)
cam0t.set(4,height)
while(1):
	#read in the image from the frame grabber
	filename=datetime.datetime.utcnow().strftime("%Y-%m-%d-%H%M%S")
	ret,frame = cam0t.read()
	for i in range(1,10):
		time.sleep(0.1)
		ret,frame2 = cam0t.read()
		frame = cv2.addWeighted(frame,0.5,frame2,0.5,0)
	cv2.imwrite(filename+'.png',frame)
	#v2.imshow("Test",frame)
	print filename
	#v2.waitKey(0)
	time.sleep(30)
cam0t.release()
