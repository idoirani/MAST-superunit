
# -*- coding: utf-8 -*-


"""
Created on Wed Dec 24 2023

@author: Ido Irani
"""

import numpy as np
from DS_config import *
import sys
sys.path.append(path_sdk)
import greateyesSDK as ge
import time
import threading
import pandas as pd
import astropy.io.fits as fits
import astropy.time as atime
import os
import logging 


class CamAPI(object):
    '''
    This class is a wrapper for the SDK for controling a single DeepSpec Camera. It provides a set of functions to control the camera and to take measurements.
    last update: 24-12-2023
    author: Ido Irani
    '''
    def __init__(self, id, IP):
        self.id = id                        # index of camera in SDK
        self.IP = IP                        # IP address of camera
        self.band = dic_band[IP]            # DeepSpec Band this Camera is assigned to
        #self.status = 'Null'                # status of camera. NOT YET IMPLEMENTED
        self.connected = False              # connection status
        self.monitor_temperature = False    # flag for monitoring temperature
        self.binning = [1,1]                # binning mode
        self.t_exp = -1                     # exposure time
        self.readout_speed = -1             # pixel readout speed
        self.safe_Temp = True               # flag for detector backside temperature safety
        self.T_goal = -80                   # goal temperature for cooling
        self.cooled = False                 # flag for cooling status
        self.initialized = False            # flag for initialization status
        self.bitpix = 0                     # number of bits per pixel
        self.logger = self.setup_logger(f'camera_{id}')
        self.connectionType = ge.connectionType_Ethernet # or ge.ConnectionType_USB instead`
        self.amp_id = int(2)
    def setup_logger(self, name, log_file=-1, level=logging.INFO):
        '''
        Sets up a logger for each instance of CamAPI.
        IN:       name                name of logger
        IN:       log_file            name of log file
        IN:       level               logging level
        OUT:      logger              logger object
        last update: 24-12-2023
        author: Ido Irani
        '''
        if log_file == -1:
            log_file =  path_logging + '\camera_log_{0}.log'.format(self.band)
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler = logging.FileHandler(log_file)        
        handler.setFormatter(formatter)

        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)

        return logger
    def DLLisbusy(self):
        '''
        check if DLL is busy
        OUT:      busy                boolean flag for busy status
        last update: 24-12-2023
        author: Ido Irani
        '''
        
        busy = ge.DllIsBusy(addr = self.id)
        return busy
        
    def camera_info(self,verbose = True):
        '''
        call some basic info on firmware, resolution, pixelsize, max exposure time, binning, capacity and output modes provided by this model of camera
        IN:       verbose             flag to print output
        RESULT:    sizex,sizey, pixsize and bitpix attributes are set
        last update: 24-12-2023
        author: Ido Irani
        '''
        size = ge.GetImageSize(addr = self.id)[0:3]
        self.sizex= size[0]
        self.sizey= size[1]
        self.bitpix = size[2]
        self.pix_size= ge.GetSizeOfPixel(addr = self.id)
        if verbose:
            self.logger.info('\n-----------------------------------------------------\n')
            self.logger.info('Gathering Camera information:')
            self.logger.info('   Firmware Version: %s', ge.GetFirmwareVersion(addr = self.id))
            self.logger.info('   Image Size: %s',size[0:2] )
            self.logger.info('   Digital Resolution: %s bit', ge.GetImageSize(addr = self.id)[2] * 8)
            self.logger.info('   Pixel Size: %s um', ge.GetSizeOfPixel(addr = self.id))
            self.logger.info('   Camera is busy: %s', ge.DllIsBusy(addr = self.id))
            self.logger.info('   max. Exposure time: %s ms', ge.GetMaxExposureTime(addr = self.id))
            self.logger.info('   max. binning x: %s y: %s', ge.GetMaxBinningX(addr = self.id), ge.GetMaxBinningY(addr = self.id))
            self.logger.info('   camera supports capacity mode: %s', ge.SupportedSensorFeature(ge.sensorFeature_capacityMode,addr = self.id))
            self.logger.info('   camera supports horizontal hardware binning: %s', ge.SupportedSensorFeature(ge.sensorFeature_binningX))
            self.logger.info('   camera supports horizontal hardware cropping: %s', ge.SupportedSensorFeature(ge.sensorFeature_cropX))
            self.logger.info('   camera provides the following output mode(s): ')
            NumberOfOutputModes = ge.GetNumberOfSensorOutputModes(addr = self.id)
            for i in range(NumberOfOutputModes):
                self.logger.info('      mode ' + str(i) + ': %s', ge.GetSensorOutputModeStrings(i,addr = self.id))
                #ge.SetupSensorOutputMode(NumberOfOutputModes-1,addr = self.id)
        pass
    

    def initialize(self):
        '''
        Initialize camera. Will try to initialize the camera at self.IP. If successful, it will set the binning mode to [1,1].
        OUT:      initialized         boolean flag for success/failure of initialization
        last update: 24-12-2023
        author: Ido Irani
        '''
        self.logger.info('initializing camera')
        self.camera_info(verbose = False)
        if (ge.InitCamera(addr=self.id) == True):
            self.initialized = True
            self.logger.info('Status:' + ge.StatusMSG)
            self.logger.info('   ok')
            _ = ge.TemperatureControl_Init(addr = self.id)
            ge.SetBitDepth(4)
            self.logger.info('Switching off LED: ' + str(self.LEDchange(False)))
            return True
        else:
            self.logger.error('   Initialization failed. Status: {0}'.format(ge.StatusMSG))
            return False
        

    def disconnect(self):
        '''
        Disconnect from camera. Will try to disconnect from the camera at self.IP. If successful, it will disconnect the camera.
        OUT:      dis                 boolean flag for success/failure of disconnection
        last update: 24-12-2023
        author: Ido Irani
        '''
        dis = ge.DisconnectCamera(addr = self.id)
        return dis
        
    def connect(self):
        '''
        Connect to camera. Will try to connect to the camera at self.IP. If successful, it will initialize the camera and set the binning mode to [1,1].
        OUT:      connected           boolean flag for success/failure of connection
        last update: 24-12-2023
        author: Ido Irani
        '''
        # setup camera connection
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('setting up connection')
        connectionSetupWorked = ge.SetupCameraInterface(self.connectionType, ipAddress=self.IP, addr =  self.id)
        if connectionSetupWorked:
            connectionSetupWorked = ge.ConnectToSingleCameraServer(addr = self.id)
            self.logger.info('using TCP connection on ip address: %s', self.IP)
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('connecting to camera')
        CameraConnected = self.connected
        CameraModel = []
        CameraConnected = ge.ConnectCamera(model=CameraModel, addr=self.id)
        if (CameraConnected == True):
            CameraModelID = CameraModel[0]
            CameraModelStr = CameraModel[1]
            self.logger.info('   connected to camera %s', CameraModelStr)
            self.connected = True
        self.logger.info('initializing camera')
        init = self.initialize()
        if not init:
            self.disconnect()
            self.connected = False
        return self.connected


    def cool_down(self, T_goal = -80, n_max = 0, verbose = True):
        '''
        Cool down to goal temperature. Will monitor until T_goal is reached and report backside and detector temperature every 10 seconds.
        n_max is the numer of reports given before minotioring terminates. Control will continute either way. 
        IN:       T_goal              goal temperature for cooling
        IN:       n_max               number of reports given before minotioring terminates. Control will continute either way.
        IN:       verbose             flag to print output
        OUT:      cooled              boolean flag for success/failure of cooling
        last update: 24-12-2023
        author: Ido Irani
        '''
        self.T_goal = T_goal
        if verbose:
            self.logger.info('\n-----------------------------------------------------\n')
            self.logger.info('initializing sensor cooling')
        Cooling_limits = ge.TemperatureControl_Init(addr = self.id)
        if verbose:
            self.logger.info('   Temperature values shall be in this range:')
            self.logger.info('      lowest possible setpoint = %s °C', Cooling_limits[0])
            self.logger.info('      highest possible setpoint = %s °C', Cooling_limits[1])
            self.logger.info('      highest possible setpoint = %s °C', Cooling_limits[1])
            self.logger.info('   Actual Temperatures are:')
            self.logger.info('      CCD (TEC frontside): = %s °C', ge.TemperatureControl_GetTemperature(0,addr = self.id))
            self.logger.info('      TEC backside: = %s °C', ge.TemperatureControl_GetTemperature(1,addr = self.id))
            self.logger.info('   Setting %s °C',self.T_goal)
        if ge.TemperatureControl_SetTemperature(self.T_goal,addr = self.id) == True:
            x = 0
            n = 0
            if verbose:
                self.logger.info('monitoring until goal temperature is reached or until max attempts')
            while x == 0:
                Temp = ge.TemperatureControl_GetTemperature(0,addr = self.id)
                backTemp = ge.TemperatureControl_GetTemperature(1,addr = self.id)
                dT = np.abs(Temp - self.T_goal)
                self.logger.info('  %s s   T_CCD: %s °C \n   T_back: %s °C', str(10*n), Temp,backTemp) 
                n+=1
                if dT<1:
                    x = 1
                    if verbose:
                        self.logger.info('reached goal Temperature. Control is active.')
                        self.cooled = True
                elif n>=n_max:
                    x = 1
                    if verbose:
                        self.logger.info('reached max attempts. Control is active.')
                if x==0:
                    time.sleep(10)
            return True
        else: 
            if verbose:
                self.logger.error('Temperature control failed')    
            return False    
            
    def warm_up(self, T_goal = 0,n_max = 0):
        '''
        warms up to T_goal and then shuts down temperature control
        IN:       T_goal              goal temperature for warming up
        IN:       n_max               number of reports given before minotioring terminates. Control will continute either way.
        last update: 24-12-2023
        author: Ido Irani
        '''
        self.T_goal = T_goal
        if ge.TemperatureControl_SetTemperature(self.T_goal,addr = self.id) == True:
            x = 0
            n = 0
            self.logger.info('warming up')
            self.cooled = False
            while x == 0:
                Temp = ge.TemperatureControl_GetTemperature(0,addr = self.id)
                backTemp = ge.TemperatureControl_GetTemperature(1,addr = self.id)
                dT = np.abs(Temp - self.T_goal)
                self.logger.info('  %s s   T_CCD: %s °C \n   T_back: %s °C', str(10*n), Temp,backTemp) 

                if dT<1:
                    x = 1
                elif n>=n_max:
                    x = 1
                if x==0:
                    time.sleep(10)
        off = ge.TemperatureControl_SwitchOff(addr = self.id)
        if ge.TemperatureControl_SwitchOff(addr = self.id) == True:
            self.logger.info('Sensor cooling switched off')
            self.cooled = False
        else: 
            self.logger.error('Sensor cooling switch off failed')
        return off
    
    def set_readout_amps(self,amp_id = -1):
        '''
        setup readout amplifiers. rotate between left and right amplifiers for readout or both. 
        # IN:       sensorOutputMode (int)        0 = 1 Amp (OSR); 1 = 1 Amp (OSL); 2 = 2 Amp (OSR & OSL)
        # In:       addr                    index of connected devices; begins at addr = 0 for first device
        # Result:   Bool                    success true/false
        '''
        if amp_id == -1:
            amp_id = self.amp_id
        self.amp_id = int(amp_id)
        ge.SetupSensorOutputMode(self.amp_id, addr = self.id)
        pass
        
    def set_binning(self, binning):
        '''
        Set binning mode. 
        # IN:       binning             list of 2 elements [binX, binY], where binX is the number of points to be binned in the X direction
        # OUT:      statusMSG           updates index and string of status message
        # In:       addr                index of connected devices; begins at addr = 0 for first device
        # Result:   Bool                success true/false
        last update: 24-12-2023
        author: Ido Irani
        '''
        res = ge.SetBinningMode(binning[0], binning[1],addr = self.id) == True
        self.logger.info('\n   Setting Binning at %s x %s',int(binning[0]), int(binning[1]))
        if res:
            self.logger.info('   Successfully set binning to  %s x %s',int(binning[0]), int(binning[1]))
            self.binning = binning
        else:
            self.logger.error('   Error while setting binning  %s x %s',int(binning[0]), int(binning[1]))
        return res 
    
    def set_exposure_time(self,t_exp, **kwargs):
        '''
        # sets the exposure time for measurements
        # IN:       t_exp               exposure time [s]
        # OUT:      statusMSG           updates index and string of status message
        # Result:   Bool                success true/false
        last update: 24-12-2023
        author: Ido Irani
        '''
        t_exp_ms = int(t_exp*1000)
        exp = ge.SetExposure(t_exp_ms,addr = self.id,**kwargs) == True
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('setting measurement parameters:')
        if (exp):
            self.logger.info('   exposure time set to %s s', t_exp)
        self.logger.info('\n-----------------------------------------------------\n')        
        self.t_exp = t_exp
        return exp
    

    def set_readout_speed(self, readoutSpeed = 7):
        '''
        # sets the readout frequency (pixels/sec) for measurements
        # IN:       readoutSpeed        
        #                                   0 -> 1 MHz
        #                                   3 -> 3 MHz
        #                                   4 -> 500 kHz
        #                                   5 -> 250 kHz
        #                                   6 -> 100 kHz
        #                                   7 -> 50 kHz
        # IN:       addr                index of connected devices; begins at addr = 0 for first device
        # OUT:      statusMSG           updates index and string of status message
        # Result:   Bool                success true/false
        last update: 24-12-2023
        author: Ido Irani
        '''
        speed = {0: '1 MHz', 3: '3 MHz', 4:'500 kHz', 5:'250 kHz',6:'100 kHz',7:'50 kHz'}
        rdspeed =  {0: ge.readoutSpeed_1_MHz,
                    3: ge.readoutSpeed_3_MHz,
                    4:ge.readoutSpeed_500_kHz,
                    5: ge.readoutSpeed_250_kHz,
                    6:ge.readoutSpeed_100_kHz, 
                    7:ge.readoutSpeed_50_kHz}
        if readoutSpeed not in speed.keys():
            raise ValueError('not in the list of keys for permitted readout modes. These are: {0: 1 MHz, 3: 3 MHz, 5:500 kHz,6:50 kHz}')
        self.logger.info('\n  Setting pixel readout frequency to %s',speed[readoutSpeed])
        res = ge.SetReadOutSpeed(rdspeed[readoutSpeed],addr = self.id)
        if  res:
            self.logger.info('   Success')
        else:
            self.logger.error('   Error while setting %s',speed[readoutSpeed])
        self.readout_speed = readoutSpeed
        return res
    def Detector_Temperature(self):
        '''
        # gives the current detector and backside temperature
        # 
        # OUT:      Temperature         a tuple containing detector and backside temepratures 
        last update: 24-12-2023
        author: Ido Irani
        '''
        temperature = ge.TemperatureControl_GetTemperature(0, addr = self.id)
        temperature_back = ge.TemperatureControl_GetTemperature(1, addr = self.id)
        return temperature, temperature_back
        
    def log_T_func(self,path_fold, out_name = -1):
        '''
        A function to monitor the temperature and logs it to a file every 100 seconds in 10 second intervals
        Will stop if self.monitor_temperature is set to False.
        IN:      path_fold           path to folder where temperature log is saved
        IN:      out_name            name of output file. Default is  "temperature_log_#.txt" where # is the band name
        last update: 24-12-2023
        author: Ido Irani
        '''
        if out_name == -1:
            out_name = 'temperature_log_' + self.band + '.txt'
        columns = ['current_time', 'detector_temperature', 'backside_temperature']
        data = pd.DataFrame(columns=columns)
        n = 1
        if path_fold[-1]== "\\":
            out_path = path_fold + out_name
        else:
            out_path = path_fold + "\\"+ out_name    
        
        while self.monitor_temperature:
            temperature = ge.TemperatureControl_GetTemperature(0, addr = self.id)
            temperature_back = ge.TemperatureControl_GetTemperature(1, addr = self.id)
                
            # Open, write to, and close the file
            current_time = time.strftime("%Y-%m-%d %H:%M:%S", time.localtime())
            #file.write(f"{current_time}, {temperature}, {temperature_back}\n")
            new_data = pd.DataFrame([[current_time, temperature, temperature_back]], columns=columns)
            data = pd.concat([data, new_data], ignore_index=True)
            
            if temperature - self.T_goal<=2:
                self.cooled = True
            else: 
                self.cooled = False
            time.sleep(10)  # Wait for 10 seconds
            if np.mod(n,10)==0:
                try:
                    data.to_csv(out_path, index=False)
                except:
                    pass
                
            if temperature_back >= 55:
                self.safe_Temp = False
                self.T_goal += 5
                set_T = ge.TemperatureControl_SetTemperature(self.T_goal, addr = self.id)
                self.logger.warning('Backside temperature exceed 55 C. Raising set temperature by 5 degrees')
                time.sleep(20)
                if not set_T:
                    try:
                        data.to_csv(path_fold + out_name, index=False)
                    except:
                        pass    
                    ge.StopMeasurement(addr = self.id)
                    self.disconnect() #replace with hard shut sown
                    self.logger.error('Temperature control failed, shutting down camera')
                if self.T_goal>-60:
                    self.logger.error('After several attempts, set point exceeding -60. Teminating')
                    ge.StopMeasurement(addr = self.id)
                    self.disconnect() #replace with hard shut sown

            n+=1
        #file.close()
        data.to_csv(path_fold + out_name, index=False)
        pass

    
    def stop_logging(self):
        '''
        stops the temperature logging thread
        RESULT:    monitor_temperature  is set to False
        last update: 24-12-2023
        author: Ido Irani
        '''
        self.monitor_temperature = False
        pass
    
    def cool_and_log_temperature(self, T_goal,path_fold,out_name = -1):
        '''
        starts a thread that monitors the temperature and logs it to a file every 100 seconds in 10 second intervals
        Will stop if self.monitor_temperature is set to False.
        
        IN:       T_goal              goal temperature for cooling
        IN:       path_fold           path to folder where temperature log is saved
        IN:       out_name            name of output file. Default is  "temperature_log_#.txt" where # is the band name
        OUT:      temp_thread         thread that monitors the temperature
        last update: 24-12-2023
        author: Ido Irani
        '''
        if not os.path.exists(path_fold):
                msgstr = "Directory {0} does not exist.".format(path_fold)
                self.logger.error(msgstr)
                return msgstr
        self.T_goal = T_goal
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('initializing sensor cooling')
        Cooling_limits = ge.TemperatureControl_Init(addr = self.id)
        self.logger.info('   Temperature values shall be in this range:')
        self.logger.info('      lowest possible setpoint = %s °C', Cooling_limits[0])
        self.logger.info('      highest possible setpoint = %s °C', Cooling_limits[1])
        self.logger.info('   Actual Temperatures are:')
        T,T_back = self.Detector_Temperature()
        self.logger.info('      CCD (TEC frontside): %s °C', T)
        self.logger.info('      TEC backside: %s °C', T_back)  
        TconSet = ge.TemperatureControl_SetTemperature(self.T_goal, addr = self.id)
        if TconSet:
            self.monitor_temperature = True
            temp_thread = threading.Thread(target=self.log_T_func, args = [path_fold, out_name])    
            temp_thread.start()
        else: 
            self.logger.error('failed setting temperature')
        return temp_thread
    def LEDchange(self, status):
        '''
        # Initiates a measurement with the current settings and returns the image data and header. 
        # IN:       t_exp            exposure time [s]. Default is the current setting self.t_exp
        # OUT:      hdul             FITS file with image data and header
        last update: 24-12-2023
        author: Ido Irani
        '''
        success = ge.SetLEDStatus(status, addr = self.id)
        if success:
            self.LED = status
        return success

    def expose(self,t_exp = -1, verbose=True):
        '''
        # Initiates a measurement with the current settings and returns the image data and header. 
        # IN:       t_exp            exposure time [s]. Default is the current setting self.t_exp
        # IN:       verbose          flag to print output
        # OUT:      hdul             FITS file with image data and header
        last update: 24-12-2023
        author: Ido Irani
        '''
        if not self.initialized: 
            self.initialize()
        if not self.connected:
            self.logger.error("Error: camera not connected.")
            return None
        if t_exp == -1:
            t_exp = self.t_exp
        if t_exp != self.t_exp:    
            _ = self.set_exposure_time(t_exp)
        if not self.cooled: 
            self.logger.warning('Warning: detector is not cooled. Signal to noise may be severly affected')
        start_time = time.localtime() 
        timer = time.time()
        self.LEDchange(False)
        time.sleep(0.5)
        measure = ge.StartMeasurement_DynBitDepth(addr = self.id)
        if measure:
            if verbose:
                self.logger.info('   measurement started')
                self.logger.info('   waiting, while DLL is busy')
            t_meas = 0
            t_crit = t_exp + 25  # seconds for measurement timeout
            dt = 0.1
            time.sleep(1)
            while self.DLLisbusy():
                time.sleep(dt)
                t_meas = time.time() - timer
                if (t_meas >= t_crit):  # if measurement takes took long
                    if ge.StopMeasurement(addr = self.id) == True:
                        import ipdb; ipdb.set_trace()
                        self.logger.error('measurement stopped, measurement took too long')
                        #return False
            end_time = time.localtime()
            utc_start_time = time.gmtime(time.mktime(start_time))
            utc_end_time = time.gmtime(time.mktime(end_time))

            mid_timestamp = (time.mktime(start_time) + time.mktime(end_time))/2
            utc_mid_time = time.gmtime(mid_timestamp)
            mid_time = time.localtime(mid_timestamp) 
            
            start_time = time.strftime("%Y-%m-%dT%H:%M:%S", start_time)
            mid_time = time.strftime("%Y-%m-%dT%H:%M:%S", mid_time)
            end_time = time.strftime("%Y-%m-%dT%H:%M:%S", end_time)
            utc_start_time = time.strftime("%Y-%m-%dT%H:%M:%S", utc_start_time)
            utc_end_time = time.strftime("%Y-%m-%dT%H:%M:%S", utc_end_time)
            utc_mid_time = time.strftime("%Y-%m-%dT%H:%M:%S", utc_mid_time)
            #ge.SetLEDStatus(True, addr = self.id)

            if verbose:
                self.logger.info('   ...finished\n')
            imageArray = ge.GetMeasurementData_DynBitDepth(addr = self.id)
            hdr = {}
            for key in FITS_std_head.keys():
                hdr[key] = FITS_std_head[key]
            hdr['BAND'] = 'DeepSpec_' + self.band
            hdr['CAMERA_IP'] = self.IP
            hdr['TYPE'] = 'RAW'
            hdr['LOCAL_T_START'] = start_time
            hdr['LOCAL_T_MID']   = mid_time
            hdr['LOCAL_T_END']   = end_time
            hdr['T_START']       = utc_start_time
            hdr['T_MID']         = utc_mid_time
            hdr['T_END']         = utc_end_time
            hdr['T_EXP']         = t_exp
            hdr['TEMP_GOAL'] = self.T_goal
            hdr['TEMP_SAFE_FLAG'] = self.safe_Temp
            hdr['DATE-OBS'] = mid_time
            hdr['MJD-OBS'] = atime.Time(utc_mid_time).mjd
            hdr['READOUT_SPEED'] = speed[self.readout_speed] 
            hdr['CDELT1'] = self.binning[0]
            hdr['CDELT2'] = self.binning[1]
            hdr['NAXIS'    ] = 2
            hdr['NAXIS1'    ] = self.sizex
            hdr['NAXIS2'    ] = self.sizey
            hdr['PIXEL_SIZE'] = self.pix_size
            hdr['BITPIX'] = self.bitpix
            for key in hdr.keys():
                hdr[key] = (hdr[key], FITS_HEADER_comment[key])
            HDR = fits.Header()
            for key in hdr.keys():
                HDR[key] = hdr[key]
            hdu = fits.PrimaryHDU(imageArray, header = HDR)
            hdul = fits.HDUList([hdu])
            return hdul
 
        

class DeepSpecAPI(object):
    '''
    This class is a wrapper for the SDK for controling the DeepSpec assembly. It provides a set of functions to control the instrument and to take measurements simultaneously with all cameras.
    last update: 24-12-2023
    author: Ido Irani
    '''
    def __init__(self, IP = IPS):
        self.IP = IP
        self.n_cams = len(IP)
        self.connected = False
        self.initialized = False
        self.monitor_T = False
        self.logger = self.setup_logger('DeepSpec')
        self.connectionType = ge.connectionType_Ethernet # or ge.ConnectionType_USB instead`
        self.cameras = []
        self.cooled = False
        for i in range(self.n_cams):
            self.cameras.append(CamAPI(i, IP[i]))
        self.bands = [self.cameras[i].band for i in range(self.n_cams)]
        self.T_goal = -80
    
    def setup_logger(self, name, log_file=-1, level=logging.INFO):
        '''
        Sets up a logger for each instance of CamAPI.
        IN:       name                name of logger
        IN:       log_file            name of log file
        IN:       level               logging level
        OUT:      logger              logger object
        last update: 24-12-2023
        author: Ido Irani
        '''
        if log_file == -1:
            log_file =  path_logging + '\DeepSpec_log.log'
        formatter = logging.Formatter('%(asctime)s - %(levelname)s - %(message)s')
        handler = logging.FileHandler(log_file)
        handler.setFormatter(formatter)

        logger = logging.getLogger(name)
        logger.setLevel(level)
        logger.addHandler(handler)
        return logger
    
    def initialize(self):
        '''
        Initialize camera. Will try to initialize the camera at self.IP. If successful, it will set the binning mode to [1,1].
        OUT:      initialized         boolean flag for success/failure of initialization
        Result:   Bool                change self.initialized to True/False
        last update: 24-12-2023)     
        author: Ido Irani
        '''
        self.logger.info('initializing DeepSpec')
        success = np.array([False]*self.n_cams)
        for i in  range(self.n_cams):
            self.logger.info('initializing camera %s', self.bands[i])
            init = self.cameras[i].initialize()
            success[i] = init
            if not init:
                self.logger.error('   Initialization failed. Status: {0}'.format(ge.StatusMSG))
        if success.all():
            self.initialized =  True
        else:
            self.initialized =  False  
        return self.initialized
    
    def disconnect(self):
        '''
        Disconnect from DeepSpec. Will try to disconnect from the camera at self.IP. If successful, it will disconnect the camera.
        OUT:      dis                 boolean flag for success/failure of disconnection
        last update: 24-12-2023
        author: Ido Irani
        '''
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('disconnecting DeepSpec')
        dis = np.array([False]*self.n_cams)
        for i in  range(self.n_cams):
            self.logger.info('disconnecting camera %s', self.bands[i])
            dis[i] = self.cameras[i].disconnect()
            if dis[i]:
                self.logger.info('   ok')
            else:
                self.logger.error('   failed')
        if dis.all():
            self.logger.info('   DeepSpec disconnected')
            self.connected =  False
        else:
            self.logger.error('   DeepSpec disconnect failed')
            self.connected =  True
        return dis
    
    def connect(self):
        '''
        Connect to DeepSpec. Will try to connect to the camera at self.IP. If successful, it will initialize the camera and set the binning mode to [1,1].
        OUT:      connected           boolean flag for success/failure of connection
        last update: 24-12-2023
        author: Ido Irani
        '''
        # setup camera connection
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('setting up connection')
        connected = np.array([False]*self.n_cams)
        for i in range(self.n_cams):
            self.cameras[i].connect()
            connected[i] = self.cameras[i].connected
            if connected[i]:
                self.logger.info('   ok')
            else:
                self.logger.error('   failed')
        if connected.all():
            self.logger.info('   DeepSpec connected')
            self.connected =  True
        else:
            self.logger.error('   DeepSpec connection failed')
            self.connected =  False
        return self.connected
    def get_temperature(self):
        T_out = []
        for i in range(self.n_cams):
            T_out.append(self.cameras[i].Detector_Temperature())
        return T_out
    def set_T_goal(self, T_goal):
        res_all = []
        for cam in self.cameras:
            res = cam.cool_down(T_goal , 0, False)
            if res: 
                cam.T_goal = T_goal
            res_all.append(res)
        if np.array(res_all).all():
            self.T_goal = T_goal
            return True
        else: 
            return False
    def cool(self, T_goal = -80, n_max = 0, verbose = True, path_fold = -1):
        '''
        Cool down to goal temperature. Will monitor until T_goal is reached and report backside and detector temperature every 10 seconds.
        n_max is the numer of reports given before minotioring terminates. Control will continute either way. 
        IN:       T_goal              goal temperature for cooling
        IN:       n_max               number of reports given before minotioring terminates. Control will continute either way.
        IN:       verbose             flag to print output
        OUT:      cooled              boolean flag for success/failure of cooling
        last update: 24-12-2023
        author: Ido Irani
        '''
        self.T_goal = T_goal
        if path_fold == -1:
            #add to path date folder
            path_fold = path_T_logging + '\\' + time.strftime("%Y%m%d", time.localtime())
            
        if not os.path.exists(path_fold):
            os.makedirs(path_fold)
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('cooling DeepSpec cameras')
        cooled = np.array([False]*self.n_cams)
        for i in range(self.n_cams):
            out_name = 'temperature_log_' + self.cameras[i].band + '.txt'
            self.cameras[i].cool_down(self.T_goal,0,False)
            cooled[i] = self.cameras[i].cool_and_log_temperature(T_goal,path_fold)
        if cooled.all():
            self.monitor_T = True
        return cooled
        
    def slit_pos(self, pos):
        '''
        PLACEHOLDER for slit control
        '''
        return True
    def start_up(self, t_exp_s = 1e-3, binning = [1,1], readoutSpeed = 7):
        '''
        Startup DeepSpec. 
        OUT:      ready           boolean flag for readeness of DeepSpec
        last update: 24-12-2023
        author: Ido Irani
        '''
        self.logger.info('\n-----------------------------------------------------\n')
        self.logger.info('starting up DeepSpec')
        if not self.connected:
            self.logger.info('connecting cameras')
            _   = self.disconnect()
            out = self.connect()
            if not out:
                return False
        if not self.cooled:
            self.logger.info('cooling...')
            _ = self.cool()
            self.cooled = True
        pos = self.slit_pos(0)
        if not pos:
            return False
        self.logger.info('setting exposure parameters')
        self.binning = binning
        self.t_exp = t_exp_s
        self.readoutSpeed = readoutSpeed
        for i in range(self.n_cams):
            self.cameras[i].set_readout_amps()
            self.cameras[i].set_binning(binning)
            self.cameras[i].set_exposure_time(t_exp_s)
            self.cameras[i].set_readout_speed(readoutSpeed)
        return True

     
    def expose(self, t_exp = -1, verbose = True):
        '''
        # Initiates a measurement with the current settings and returns the image data and header. 
        # IN:       t_exp            exposure time [s]. Default is the current setting self.t_exp
        # IN:       verbose          flag to print output
        # OUT:      hdul             FITS file with 4 extentions with image data and headers
        last update: 24-12-2023)
        author: Ido Irani
        '''
        hdul = [None] * self.n_cams

        def thread_function(camera, hdul, index):
            hdu = camera.expose(t_exp, verbose)
            hdul[index] = hdu[0]
        # Create a list to hold threads
        threads = []
        # Start a new thread for each camera
        if self.n_cams>1:
            for i in range(self.n_cams):
                thread = threading.Thread(target=thread_function, args=(self.cameras[i], hdul, i))
                threads.append(thread)
                thread.start()
            # Wait for all threads to complete
            for thread in threads:
                thread.join()

        else:
            camera = self.cameras[0]
            hdu = camera.expose(t_exp, verbose)
            hdul[0] = hdu[0]
        hdull = fits.HDUList(hdus=hdul)

        return hdull
    
