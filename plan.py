# -*- coding: utf-8 -*-
"""
Module to store Plans and routines for DeepSpec
@author: idoi
"""


from DS_config import *
import sys
sys.path.append(path_sdk)
from DeepSpec import *
import astropy.io.fits as fits 
import time
import os
import numpy as np 
import ctypes
import pyautogui


# Define necessary constants
SC_MONITORPOWER = 0xF170
MONITOR_OFF = 2
MONITOR_ON = -1
HWND_BROADCAST = 0xFFFF
WM_SYSCOMMAND = 0x0112


# Load the user32 DLL
user32 = ctypes.WinDLL('user32')

def turn_off_display():
    # Sends the system command to turn off the display
    user32.SendMessageW(HWND_BROADCAST, WM_SYSCOMMAND, SC_MONITORPOWER, MONITOR_OFF)
    pass 
def turn_on_display():
    # Simulate mouse movement to turn on the display
    pyautogui.moveRel(0, 1)
    pyautogui.moveRel(0, -1)
    pass






class ExposeAction(object):
    def __init__(self, DS_obj,out_fold, outfile, t_exp, verbose, type = ''):
        self.DS_obj = DS_obj
        self.t_exp = t_exp
        self.verbose = verbose
        self.type = type
        self.out_fold = out_fold
        #add backslash to out_fold if needed according to OS
        if out_fold.endswith(os.sep):
            self.out_fold = out_fold
        else:
            self.out_fold = out_fold + os.sep
        if outfile.startswith(os.sep):
            self.outfile = self.out_fold + outfile.lstrip(os.sep)
        else:
            self.outfile = self.out_fold + outfile
    def execute(self, write = True):
        hdul = self.DS_obj.expose(self.t_exp, self.verbose)
        if self.type != '':
            for hdu in hdul:
                hdu.header['TYPE'] = self.type
        
        if not os.path.exists(self.out_fold):
            os.mkdir(self.out_fold)
        # save hdul fits file to disk
        if write:
            try:
                hdul.writeto(self.outfile, overwrite=True)
            except Exception as e:
                print('Write to file failed. Trying to rewrite with a new file name')
                new_outfile = self.outfile
                new_outfile = new_outfile.split('.fits')[0] + '_(1).fits'
                hdul.writeto(new_outfile, overwrite=True)
        return hdul


class BiasAction(object):
    def __init__(self, DS_obj, out_fold, outfile, n_exp = 10, t_exp = 1, verbose = True):
        self.DS_obj = DS_obj
        self.t_exp = t_exp
        self.verbose = verbose
        self.n_exp = n_exp
        self.bands = DS_obj.bands
        if out_fold.endswith(os.sep):
            self.out_fold = out_fold
        else:
            self.out_fold = out_fold + os.sep
        if outfile.startswith(os.sep):
            self.outfile = self.out_fold + outfile.lstrip(os.sep)
        else:
            self.outfile = self.out_fold + outfile
            
    def execute(self, write = True):
        hdul = [[]]*len(self.bands)
        # expose n_exp bias frames
        for i in range(self.n_exp):
            hdull = self.DS_obj.expose(self.t_exp, self.verbose)
            for j in range(len(hdull)):
                hdull[j].header['TYPE'] = 'BIAS'
                hdul[j] = hdul[j] + [hdull[j]]
        # create a master frame from the median of all frames
        master_hdul = fits.HDUList(hdus=[])
        for i in range(len(hdul)):
            all_exps = np.array([hdul[i][j].data for j in range(len(hdul[i]))])
            med_exp = np.median(all_exps, axis=0)
            hdr = hdul[i][len(hdul)//2].header
            out_hdul = fits.HDUList([fits.PrimaryHDU(med_exp, header=hdr)])
            master_hdul.append(out_hdul[0])  
        if not os.path.exists(self.out_fold):
            os.mkdir(self.out_fold)
        # save hdul fits file to disk
        master_hdul.writeto(self.outfile, overwrite=True)
        return master_hdul
            
class FullBiasAction(object):
    def __init__(self, DS_obj, out_fold, outfile, n_exp = 10, t_exp = 1, verbose = True):
        self.DS_obj = DS_obj
        self.t_exp = t_exp
        self.verbose = verbose
        self.n_exp = n_exp
        self.bands = DS_obj.bands
        if out_fold.endswith(os.sep):
            self.out_fold = out_fold
        else:
            self.out_fold = out_fold + os.sep
        if outfile.startswith(os.sep):
            self.outfile = self.out_fold + outfile.lstrip(os.sep)
        else:
            self.outfile = self.out_fold + outfile
            
    def execute(self, write = True):
        hdul = [[]]*len(self.bands)
        # expose n_exp bias frames
        for i in range(self.n_exp):
            hdull = self.DS_obj.expose(self.t_exp, self.verbose)
            if i==0:
                time.sleep(0.5)
                for m in range(DS_obj.n_cams):
                    print('DLL is busy? ({0} band): '.format(DS_obj.bands[m])+str(DS_obj.cameras[m].DLLisbusy()))
            for j in range(len(hdull)):
                hdull[j].header['TYPE'] = 'BIAS'
                hdul[j] = hdul[j] + [hdull[j]]
                # create a master frame from the median of all frames
        HDU_dic = {}
        for i,cam in enumerate(self.bands):
            HDU_dic[cam] = hdul[i]
        if not os.path.exists(self.out_fold):
            os.mkdir(self.out_fold)
        # save hdul fits file to disk
        for i,cam in enumerate(self.DS_obj.bands):
            try:
                for j in range(self.n_exp):
                    hdu_list = fits.HDUList(hdus=HDU_dic[cam][j])
                    hdu_list.writeto(self.outfile.split('.fits')[0] + '_'+cam+'_exp'+f'{j}'+'.fits', overwrite=True)
            except:
                import ipdb; ipdb.set_trace()
        return HDU_dic

class IdleAction(object):
    def __init__(self, duration):
        self.duration = duration

    def execute(self):
        time.sleep(self.duration)

class CoolToGoalAction(object):
    def __init__(self, DS_obj, interval = 30):
        self.DS_obj = DS_obj
        self.interval = interval
    def execute(self):
        self.T_goal = self.DS_obj.T_goal
        x = 1 
        while x==1:
            Temps = self.DS_obj.get_temperature()
            Temps_front = np.array([Temps[i][0] for i in range(len(Temps))])
            if self.T_goal<-70:
                diff_T= np.sum(np.abs(Temps_front - self.T_goal))
                max_diff = np.max(np.abs(Temps_front - self.T_goal))
            else: 
                diff_T= np.sum(np.abs(Temps_front - self.T_goal))
                max_diff = np.max(np.abs(Temps_front - self.T_goal))

            if (diff_T>5)|(max_diff>3): 
                print('Detector Temmperature:')
                for i,cam in enumerate(DS_obj.bands):
                    print('{0} = {1} °C'.format(cam,Temps[i][0]))
                time.sleep(self.interval)
            else:
                x = 0
                print('Detector Temmperature:')
                for i,cam in enumerate(DS_obj.bands):
                    print('{0} = {1} °C'.format(cam,Temps[i][0]))
                #print('Detector Temmperature: U = {0} °C, G = {1} °C, R = {2} °C, I = {3} °C'.format(Temps[0][0],Temps[1][0],Temps[2][0],Temps[3][0]))
                print('goal temperature reached to within threshold')
        pass         
            
class StartUpAction(object):
    def __init__(self):
        self.DS_obj = DeepSpecAPI()
    def execute(self):
        res = self.DS_obj.start_up()
        if not res: 
            print('failed to start cameras. Try reseting power')
            self.DS_obj.disconnect()
            sys.exit(1)
        else:
            return self.DS_obj


class CoolAction(object):
    def __init__(self, DS_obj, T_goal):
        self.DS_obj = DS_obj
        self.T_goal = T_goal
    def execute(self):
        if not self.DS_obj.monitor_T:
            self.DS_obj.cool(self.T_goal)
        else:
            self.DS_obj.set_T_goal(self.T_goal)
 
class DisconnectAction(object):
    def __init__(self, DS_obj):
        self.DS_obj = DS_obj
    def execute(self):
        self.DS_obj.disconnect()
    pass





class Plan(object):
    def __init__(self,DS_obj, filename,path_fold,write = False, screen_off = True):
        self.actions = []
        self.action_type_dic = {'idle': 'misc.', 'start_up':'misc.', 'cool':'misc.', 'cool_to_goal':'misc.', 'expose':'expose', 'bias':'expose', 'fullbias':'expose'}
        self.action_type = []
        self.parse_plan(filename,path_fold)
        self.action_number = 0
        self.write = write
        self.DS_obj = DS_obj
    def parse_plan(self, filename, path_fold):
        with open(filename, 'r') as file:
            for line in file:
                parts = line.strip().split(', ')
                if not parts or parts[0].startswith('#'):  # Skip empty lines or comments
                    continue
                action = parts[0]
                self.action_type.append(action)
                if action == 'cool':
                    self.actions.append(CoolAction(DS_obj, int(parts[1])))
                if action == 'cool_to_goal':
                    self.actions.append(CoolToGoalAction(DS_obj, int(parts[1])))
                elif action == 'idle':
                    self.actions.append(IdleAction(float(parts[1])))
                elif action == 'bias':
                    self.actions.append(BiasAction(DS_obj, parts[1], parts[2], int(parts[3]), float(parts[4])))
                elif action == 'fullbias':
                    self.actions.append(FullBiasAction(DS_obj, parts[1], parts[2], int(parts[3]), float(parts[4])))
                elif action == 'expose':
                    self.actions.append(ExposeAction(DS_obj, parts[1], parts[2], float(parts[3]), parts[4]))
                elif action not in ['cool','cool_to_goal', 'idle', 'bias','expose' ]:
                    if action!=['']:
                        print(f"Unknown action type: {action}")
            self.actions.append(DisconnectAction(DS_obj))
        pass
    def execute(self):
        print(self.actions)
        for action,atype in zip(self.actions,self.action_type):
            if atype == 'cool':
                T_goal = action.T_goal
                print(f'Setting cool goal to {T_goal} C')
            if atype == 'cool_to_goal':
                print(f'Cooling until goal temperature is reached')
            elif atype == 'idle':
                dur = action.duration
                print(f'pausing for {dur} s')
            elif atype == 'bias':
                n_exp = action.n_exp
                t_exp = action.t_exp
                print(f'bias with {n_exp} exposures of {t_exp} s')
            elif atype == 'fullbias':
                n_exp = action.n_exp
                t_exp = action.t_exp
                print(f'Full bias with {n_exp} exposures of {t_exp} s')
            elif atype == 'expose':
                t_exp = action.t_exp
                print(f'exposure {t_exp} s')
            self.action_number += 1
            if self.action_type_dic[atype] == 'misc.':
                action.execute()
            elif self.action_type_dic[atype] == 'expose':
                #turn_off_display()
                action.execute(write = self.write)
                #turn_on_display()
        pass
    
if __name__ == "__main__":
    if len(sys.argv) != 4:
        print("Usage: python plan.py <filename> <path_fold_save> <write? True/False>")
        sys.exit(1)
    filename = sys.argv[1]
    path_fold = sys.argv[2]
    write = sys.argv[3]
    DS_obj = DeepSpecAPI()
    print('starting-up....')
    DS_obj.start_up()
    plan = Plan(DS_obj,filename,path_fold, write)
    try:
        plan.execute()
        print('Finished plan successfully')
    except Exception as e:
        import ipdb; ipdb.set_trace()
        plan.DS_obj.disconnect()
        print(f"An error occurred: {e}")
        sys.exit(1)

            
            