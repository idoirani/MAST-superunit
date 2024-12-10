# -*- coding: utf-8 -*-
"""
Configuration file for DeepSpec.py class definitions
"""

path_sdk = r'C:\Users\idoi\Dropbox\MAST\Deep_Spec\detectors\DeepSpec_control'
path_logging =  r'C:\Users\idoi\Dropbox\MAST\Deep_Spec\detectors\DeepSpec_control\log'
path_T_logging = r'C:\Users\idoi\Dropbox\MAST\Deep_Spec\detectors\DeepSpec_control\log'

CameraIP_1 = '192.168.1.231'
CameraIP_2 = '192.168.1.232'
CameraIP_3 = '192.168.1.233'
CameraIP_4 = '192.168.1.234'
IPS = [CameraIP_1, CameraIP_2, CameraIP_4, CameraIP_3]#, CameraIP_1]
#IP/band mapping
dic_band =  {'192.168.1.231' : 'U'
			,'192.168.1.232' : 'G'
			,'192.168.1.234' : 'R'
			,'192.168.1.233' : 'I'}
#serial_numer/band mapping
serial_numer_band =  {'192.168.1.231' :'ELSEi100365' 
					 ,'192.168.1.232' :'ELSEi100366' 
					 ,'192.168.1.233' :'ELSEi100367' 
					 ,'192.168.1.231' :'ELSEi100368'}


#standard header key/value paris
FITS_std_head = {'NAXIS':'2',
				 'TELESCOPE':'WAO-MAST',
				 'INSTRUMENT':'DEEPSPEC'
				 ,'DETECTOR':'GE 1024 1024 BI DD'}

#readout speed keyword dictionary
speed = {-1: 'NULL', 0: '1 MHz', 3: '3 MHz', 4:'500 kHz',5:'250 kHz', 6:'100 kHz',7:'50 kHz'}


FITS_HEADER_comment =   {'TELESCOPE': 'TELESCOPE NAME'
                        ,'INSTRUMENT': 'INSTRUMENT NAME'
                        ,'DETECTOR': 'DETECTOR MODEL'
						,'BAND': 'DEEPSPEC BAND' 
						,'CAMERA_IP': 'CAMERA IP' 
						,'TYPE': 'EXPOSURE TYPE' 
						,'LOCAL_T_START':  'EXPOSURE START TIME [local]' 
						,'LOCAL_T_MID':    'EXPOSURE MID TIME [local]'
						,'LOCAL_T_END':    'EXPOSURE END TIME [local]'
						,'T_START': 'EXPOSURE START TIME [UTC]' 
						,'T_MID':   'EXPOSURE MID TIME [UTC]'
						,'T_END':   'EXPOSURE END TIME [UTC]'
						,'T_EXP': 'TOTAL INTEGRATION TIME' 
						,'TEMP_GOAL': 'GOAL DETECTOR TEMPERATURE' 
						,'TEMP_SAFE_FLAG': 'DETECTOR BACKSIDE TEMPERATURE SAFETY FLAG' 
						,'DATE-OBS': 'OBSERVATION DATE' 
						,'MJD-OBS': 'MJD OF OBSERVATION MIDPOINT' 
						,'READOUT_SPEED': 'PIXEL READOUT FREQUENCY'
						,'CDELT1': 'BINNING IN THE X DIRECTION' 
						,'CDELT2': 'BINNING IN THE Y DIRECTION' 
						,'NAXIS':'NUMBER OF AXES IN FRAME'     
						,'NAXIS1':'NUMBER OF PIXELS IN THE X DIRECTION'     
						,'NAXIS2':'NUMBER OF PIXELS IN THE Y DIRECTION'    
						,'PIXEL_SIZE': 'PIXEL SIZE IN MICRONS'
                        ,'BITPIX':'# of bits storing pix values'}
	
if __name__ == '__main__':
    path_sdk = r'C:\Users\idoi\Dropbox\MAST\Deep_Spec\detectors\DeepSpec_control'
    import sys 
    sys.path.append(path_sdk)
    from DeepSpec import *
    