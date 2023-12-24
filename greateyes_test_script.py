# -*- coding: utf-8 -*-
"""
Created on Wed Mar 18 11:22:47 2020

@author: marcel
"""

import numpy as np
from PIL import Image
import sys
sys.path.append(r'C:\Users\idoi\Dropbox\DeepSpec\Detectors\sdk\python\22.5')
import greateyesSDK as ge
import time
from matplotlib.pyplot import imshow, show, colorbar
import time
# --------------------------------------------------------------------------------------------------------
# 1 Options for this example script
# --------------------------------------------------------------------------------------------------------

connectionType = ge.connectionType_Ethernet  # or ge.ConnectionType_USB instead
#connectionType = ge.connectionType_USB  # or ge.ConnectionType_USB instead

CameraIP = '192.168.1.234'  # only needed for Ethernet cameras

CameraIP_1 = '192.168.1.231'
CameraIP_2 = '192.168.1.232'
CameraIP_3 = '192.168.1.233'
CameraIP_4 = '192.168.1.234'

CameraIP  = [CameraIP_1,CameraIP_2]
CameraIP  = [CameraIP_3,CameraIP_4]


showImages = True  # displays all new images using matplotlib
saveSampleImage = True
CameraSensorCooling = True  # before setting this True, please ensure that it is safe to cool the CCD sensor.
Temperature_SetPoint = -70
binning = [2,2]
# I.e. if you are using a camera with vacuum flange or inVacuum type, make sure the chip is evacuated
# and the pressure is at least 1e-4 mbars or lower.
t_exp = 1  # ms
add_test = 0
# --------------------------------------------------------------------------------------------------------

# 2 function declarations
# --------------------------------------------------------------------------------------------------------

# This function starts a measurement and once it is finished, calls the image data from the greateyes SDK.
# It uses the non-blocking functions, so python is operating in the meantime
def Camera_PerformMeasurement(prints=True, addr = 0):
	if (ge.StartMeasurement_DynBitDepth(addr = addr) == True):
		if prints:
			print('   measurement started')
			print('   waiting, while DLL is busy')
		t_meas = 0
		global t_exp
		t_crit = t_exp / 1000 + 10  # seconds for measurement timeout
		dt = 0.01
		while ge.DllIsBusy(addr = addr):
			time.sleep(dt)
			t_meas = t_meas + dt
			if (t_meas >= t_crit):  # if measurement takes took long
				if ge.StopMeasurement(addr = addr) == True:
					print('measurement stopped, measurement took too long')
		if prints:
			print('   ...finished\n')
		imageArray = ge.GetMeasurementData_DynBitDepth(addr = addr)
		return imageArray

def display_image(imageArray):
	imshow(imageArray)
	colorbar()
	show()

def camera_connect_single():
	connectionSetupWorked = None
	image = None
	image_blocking = None
	NumberOfOutputModes = None
	# Check greateyes DLL is working
	print('DLL Version:')
	print(ge.GetDLLVersion())
	ge.DisconnectCamera(addr = 0)
	ge.DisconnectCameraServer(addr = 0)

	# setup camera connection
	print('\n-----------------------------------------------------\n')
	print('setting up connection')
	if connectionType == ge.connectionType_USB:
		print('   using USB connection')
		connectionSetupWorked = ge.SetupCameraInterface(connectionType, addr = 0)
	elif connectionType == ge.connectionType_Ethernet:
		print('using TCP connection on ip address:', CameraIP)
		connectionSetupWorked = ge.SetupCameraInterface(connectionType, ipAddress=CameraIP, addr = 0)
		if connectionSetupWorked:
			connectionSetupWorked = ge.ConnectToSingleCameraServer(addr = 0)
	else:
		print('unexpected connectionType index')
	if connectionSetupWorked == True:
		print('   ok')
	else:
		print('   failed')
		print('   Status:', ge.StatusMSG)
		sys.exit()
	print('\n-----------------------------------------------------\n')
	print('attempting to connect to camera')
	CameraConnected = False
	N_Cams = ge.GetNumberOfConnectedCams()
	print('   ' + str(N_Cams), 'camera(s) detected')
	if (N_Cams == 1):
		print('\n-----------------------------------------------------\n')
		print('connecting')
		CameraModel = []
		CameraConnected = ge.ConnectCamera(model=CameraModel, addr=0)
		if (CameraConnected == True):
			CameraModelID = CameraModel[0]
			CameraModelStr = CameraModel[1]
			print('   connected to camera ' + CameraModelStr)
			print('\n-----------------------------------------------------\n')
			print('initializing camera')
			if (ge.InitCamera(addr=0) == True):
				print('Status:' + ge.StatusMSG)
				print('   ok')
				Cooling_limits = ge.TemperatureControl_Init()
				ge.SetBitDepth(4)
				print('Switching off LED: ' + str(ge.SetLEDStatus(False)))
				ge.OpenShutter(1)
				time.sleep(1)
				ge.OpenShutter(0)
			else:
				print('   failed')
				print('   Status:', ge.StatusMSG)
				ge.DisconnectCamera()
				sys.exit()
		else:
			print('   failed')
			print('   Status:', ge.StatusMSG)
			ge.DisconnectCamera()
			if connectionType == ge.connectionType_Ethernet:
				ge.DisconnectCameraServer()
			sys.exit()
		return  CameraConnected

def camera_connect_multi(IP_list):
	connectionSetupWorked = None
	image = None
	image_blocking = None
	NumberOfOutputModes = None
	# Check greateyes DLL is working
	print('DLL Version:')
	print(ge.GetDLLVersion())
	ge.DisconnectCamera()
	ge.DisconnectCameraServer()

	# setup camera connection
	print('\n-----------------------------------------------------\n')
	print('setting up connection')
	for i,CamIP in enumerate(IP_list):
		if connectionType == ge.connectionType_Ethernet:
			print('using TCP connection on ip address:', CameraIP[i])
			connectionSetupWorked = ge.SetupCameraInterface(connectionType, ipAddress=CamIP, addr = i)
			if connectionSetupWorked:
				connectionSetupWorked = ge.ConnectToSingleCameraServer(addr = i)
		else:
			print('unexpected connectionType index')
		if connectionSetupWorked == True:
			print('   ok')
		else:
			print('   failed')
			print('   Status:', ge.StatusMSG)
			#sys.exit()
	print('\n-----------------------------------------------------\n')
	print('attempting to connect to cameras')
	CameraConnected = False
	CameraConnectedList = []
	N_Cams = len(IP_list)
	for add in range(N_Cams):
		CameraConnected = False
		#print('   failed. Number of cameras is not 1.')
		print('\n-----------------------------------------------------\n')
		print('connecting to cam {0}'.format(add))
		CameraModel = []
		CameraConnected = ge.ConnectCamera(model=CameraModel, addr=add)
		if (CameraConnected == True):
			CameraModelID = CameraModel[0]
			CameraModelStr = CameraModel[1]
			print('   connected to camera ' + CameraModelStr)
			print('\n-----------------------------------------------------\n')
			print('initializing camera')
			if (ge.InitCamera(addr=add) == True):
				print('Status:' + ge.StatusMSG)
				print('   ok')
				Cooling_limits = ge.TemperatureControl_Init(addr = add)
				ge.SetBitDepth(4)
				print('Switching off LED: ' + str(ge.SetLEDStatus(False)))
			else:
				print('   failed')
				print('   Status:', ge.StatusMSG)
				ge.DisconnectCamera(addr = add)
				#sys.exit()
		else:
			print('   failed')
			print('   Status:', ge.StatusMSG)
			ge.DisconnectCamera(addr = add)
			if connectionType == ge.connectionType_Ethernet:
				ge.DisconnectCameraServer(addr = add)
			#sys.exit()	
		CameraConnectedList.append(CameraConnected)
	return CameraConnectedList




def get_camera_info(addr = 0):
	# Get camera information
	if CameraConnectedList[addr]:
		print('\n-----------------------------------------------------\n')
		print('Gathering Camera information:')
		print('   Firmware Version:', ge.GetFirmwareVersion(addr = addr))
		print('   Image Size:', ge.GetImageSize(addr = addr)[0:2])
		print('   Digital Resolution:', ge.GetImageSize(addr = addr)[2] * 8, 'bit')
		print('   Pixel Size: ', ge.GetSizeOfPixel(addr = addr), 'um')
		print('   Camera is busy:', ge.DllIsBusy(addr = addr))
		print('   max. Exposure time:', ge.GetMaxExposureTime(addr = addr), 'ms')
		print('   max. binning x:', ge.GetMaxBinningX(addr = addr), 'y:', ge.GetMaxBinningY(addr = addr))
		print('   camera supports capacity mode:', ge.SupportedSensorFeature(ge.sensorFeature_capacityMode,addr = addr))
		print('   camera supports horizontal hardware binning:', ge.SupportedSensorFeature(ge.sensorFeature_binningX))
		print('   camera supports horizontal hardware cropping:', ge.SupportedSensorFeature(ge.sensorFeature_cropX))
		print('   camera provides the following output mode(s):')
		NumberOfOutputModes = ge.GetNumberOfSensorOutputModes(addr = addr)
		for i in range(NumberOfOutputModes):
			print('      mode ' + str(i) + ':', ge.GetSensorOutputModeStrings(i,addr = addr))
		ge.SetupSensorOutputMode(NumberOfOutputModes-1,addr = addr)
		
def cool_down(T_goal = -70, n_max = 80,addr = 0, verbose = True):
	# temperature management
	if verbose:
		print('\n-----------------------------------------------------\n')
		print('initializing sensor cooling')
	Cooling_limits = ge.TemperatureControl_Init(addr = addr)
	if verbose:
		print('   Temperature values shall be in this range:')
		print('      lowest possible setpoint =', Cooling_limits[0], '°C')
		print('      highest possible setpoint =', Cooling_limits[1], '°C')
		print('   Actual Temperatures are:')
		print('      CCD (TEC frontside):', ge.TemperatureControl_GetTemperature(0,addr = addr), '°C')
		print('      TEC backside:', ge.TemperatureControl_GetTemperature(1,addr = addr), '°C')
		print('   Setting {0} °C, monitoring for 10 minutes'.format(T_goal))
	if ge.TemperatureControl_SetTemperature(T_goal) == True:
		x = 0
		n = 0
		if verbose:
			print('warming up')
		while x == 0:
			Temp = ge.TemperatureControl_GetTemperature(0,addr = addr)
			dT = np.abs(Temp - T_goal)
			print('  ', 10*n, 's   T_CCD:', Temp, '°C', '   T_back:',
			      ge.TemperatureControl_GetTemperature(1,addr = addr), '°C')
			n+=1
			if dT<1:
				x = 1
				if verbose:
					print('reached goal Temperature. Control is active.')
			elif n>=n_max:
				x = 1
				if verbose:
					print('reached max attempts. Control is active.')
			if x==0:
				time.sleep(10)
		return True
	else: 
		if verbose:
			print('Temperature control failed')	
		return False

def warm_up(T_goal = 0,n_max = 10,addr = 0):
	if ge.TemperatureControl_SetTemperature(T_goal,addr = addr) == True:
		x = 0
		n = 0
		print('warming up')
		while x == 0:
			Temp = ge.TemperatureControl_GetTemperature(0,addr = addr)
			dT = np.abs(Temp - T_goal)
			print('  ', 10*n, 's   T_CCD:', Temp, '°C', '   T_back:',
			      ge.TemperatureControl_GetTemperature(1,addr = addr), '°C')
			n+=1
			if dT<1:
				x = 1
			elif n>=n_max:
				x = 1
			if x==0:
				time.sleep(10)
	off = ge.TemperatureControl_SwitchOff(addr = addr)
	if ge.TemperatureControl_SwitchOff(addr = addr) == True:
		print('Sensor cooling switched off')
	else: 
		print('Sensor cooling switch off failed')
	return off

def set_params(t_exp = t_exp, addr = 0 ):
	# Set Measurement Parameters
	if CameraConnectedList[addr]:
		print('\n-----------------------------------------------------\n')
		print('setting measurement parameters:')
		if (ge.SetExposure(t_exp,addr = addr) == True):
			print('   exposure time set to', t_exp, 'ms')
		print('\n-----------------------------------------------------\n')
			# binning mode
		print('\n   Setting Binning at {0} x {1}'.format(int(binning[0]), int(binning[1])))
		if ge.SetBinningMode(binning[0], binning[1],addr = addr) == True:
			print('   Successfully set binning to {0} x {1}'.format(int(binning[0]), int(binning[1])))
		else:
			print('   Error while setting binning {0} x {1}'.format(int(binning[0]), int(binning[1])))

		# readout frequency
		print('\n  Setting to 50 kHz readout')
		if ge.SetReadOutSpeed(ge.readoutSpeed_50_kHz,addr = addr) == True:
			print('   Success')
		else:
			print('   Error while setting 50 kHz')
			
def disconnect(addr = 0):
	print('\n-----------------------------------------------------\n')
	print('disconnecting')
	if (ge.DisconnectCamera(addr = addr) == True):
		print('   done')
	if connectionType == ge.connectionType_Ethernet:
		if ge.DisconnectCameraServer(addr = addr) == True:
			print('   CameraServer connection closed')
	else:
		print('   failed')
		print('   Status:', ge.StatusMSG)
		sys.exit()
		
def set_up_cam(addr = 0):
	CameraConnected = camera_connect(addr = addr)
	# Get camera information
	get_camera_info(addr = addr)
	# Set Measurement Parameters
	set_params(addr = addr)
	# cool down sensor
	cool_down(addr = addr)
	
def take_image(t_exp = 1000, showImages = False, addr = 0):
	set_params(t_exp = t_exp)
	print('\n-----------------------------------------------------\n')
	print('taking single shot non-blocking:')
	image = Camera_PerformMeasurement(addr = addr)
	print('   measurement time:   ', '{:6.3f}'.format(ge.GetLastMeasTimeNeeded()), 'seconds')
	print('   mean intensity:     ', '{:6.3f}'.format(np.mean(image)), 'ADU')
	print('   standard deviation: ', '{:6.3f}'.format(np.std(image)), 'ADU')
	if showImages:
		display_image(image)
	return image
def shut_down_cam(addr = 0):
	warm_up(addr = addr)
	disconnect(addr = addr)
	
import threading
import time
import pandas as pd






def log_T_func(addr = 0,path_fold = 'C:\\Users\\idoi\\Dropbox\\MAST\\Deep_Spec\\detectors\\Temperature_stability'):
	global monitor
	columns = ['current_time', 'detector_temperature', 'backside_temperature']
	data = pd.DataFrame(columns=columns)
	# Open file initially to write headers
	#file = open(path_fold + "\\temperature_log_{0}.txt".format(addr), "a")
	#file.write("current_time, detector_temperature, backside_temperature \n")
	n = 1
	while monitor:
		temperature = ge.TemperatureControl_GetTemperature(0, addr=addr)
		temperature_back = ge.TemperatureControl_GetTemperature(1, addr=addr)
		# Open, write to, and close the file
		current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
		#file.write(f"{current_time}, {temperature}, {temperature_back}\n")
		new_data = pd.DataFrame([[current_time, temperature, temperature_back]], columns=columns)
		data = pd.concat([data, new_data], ignore_index=True)
		time.sleep(10)  # Wait for 10 seconds
		if np.mod(n,10)==0:
			try:
				data.to_csv(path_fold + "\\temperature_log_{0}.csv".format(addr), index=False)
			except:
				pass
		n+=1
	#file.close()
	data.to_csv(path_fold + "\\temperature_log_{0}.csv".format(addr), index=False)
	pass


def stop_logging():
	global monitor
	monitor = False
	pass
def cool_and_log_temperature(T_goal = -70, addr = 0,path_fold = 'C:\\Users\\idoi\\Dropbox\\MAST\\Deep_Spec\\detectors\\Temperature_stability'):
	# temperature management
	print('\n-----------------------------------------------------\n')
	print('initializing sensor cooling')
	Cooling_limits = ge.TemperatureControl_Init(addr = addr)
	print('   Temperature values shall be in this range:')
	print('      lowest possible setpoint =', Cooling_limits[0], '°C')
	print('      highest possible setpoint =', Cooling_limits[1], '°C')
	print('   Actual Temperatures are:')
	print('      CCD (TEC frontside):', ge.TemperatureControl_GetTemperature(0,addr = addr), '°C')
	print('      TEC backside:', ge.TemperatureControl_GetTemperature(1,addr = addr), '°C')	
	if ge.TemperatureControl_SetTemperature(T_goal, addr = addr) == True:
		global monitor
		monitor = True
		temp_thread = threading.Thread(target=log_T_func, args = [addr,path_fold])	
		temp_thread.start()
	else: 
		print('failed')
	return temp_thread

# Main thread can continue to perform other tasks
# For example, sending commands to cameras


CameraConnected = camera_connect_multi(CameraIP)
if not isinstance(CameraConnected,list):
	CameraConnectedList = [CameraConnected]
else:
	CameraConnectedList = CameraConnected
N_cam = len(CameraConnectedList)	

for add in range(N_cam):
	set_params(addr = add)
	#cool_down(addr = add)

Thread_list = []
for cam in range(N_cam):
	cool_down(T_goal = -70, n_max = 1,addr = cam, verbose = False)
	temp_thread = cool_and_log_temperature(T_goal = -70, addr = cam)
	Thread_list.append(temp_thread)
for cam in range(N_cam):
	warm_up(addr = cam, n_max = 5)

for cam in range(N_cam):
	disconnect(addr = cam)

#
#
#N = 10
#image_bias = take_image(t_exp = 1, addr = add_test)
#image_bias = np.zeros((image_bias.shape[0],image_bias.shape[1],N))
#for i in range(N):
#	img = take_image(t_exp = 1, addr = add_test)
#	image_bias[:,:,i] = img
#
#image_bias_median = np.zeros_like(img)
#image_bias_std = np.zeros_like(img)
#
#for i in range(image_bias.shape[0]):
#	for j in range(image_bias.shape[1]):
#		image_bias_median[i,j] = np.median(image_bias[i,j,:])
#		p75 = np.percentile(image_bias[i,j,:] - image_bias_median[i,j], 75)
#		p25 = np.percentile(image_bias[i,j,:] - image_bias_median[i,j], 25)
#
#		image_bias_std[i,j] = p75 - p25
#
#display_image(image_bias_median)
#display_image(image_bias_std)
#
#plt.imshow(image_bias[:,:,1] - image_bias_median)