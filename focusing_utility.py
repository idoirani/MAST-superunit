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

import tkinter as tk
from tkinter import simpledialog, messagebox, filedialog
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
import matplotlib.pyplot as plt
import threading
from datetime import datetime


### GUI
class RealTimeGUI:
    def __init__(self, master, Type = 'lFWHM'):
        self.master = master
        self.master.title("Real-Time Exposure Viewer")
        self.running = False
        self.t_exp = 0.2  # Default exposure time
        self.ind = 0  # Default camera index
        self.row_range = (0, None)  # Row range for cropping (start, end)
        self.col_range = (0, None)  # Column range for cropping (start, end)
        self.thresh_signi = 100 
        self.bkg = np.zeros((1056,1027))
        self.bkg_hdul = None
        self.img_hdul = None
        self.img_sub_hdul = None
        self.connected = False
        self.cooled = False
        self.T_goal = -70
        self.save_loc = ''
        self.grid_ind = {'U':[0,0],'G':[1,0],'R':[0,1],'I':[1,1]}
        # Control Frame
        top_panel  = tk.Frame(master)
        top_panel.pack(side = 'top', pady=5)
        # Create left panel frame
        left_panel = tk.Frame(master)
        left_panel.pack(side="left", fill="y", padx=5, pady=5)
        if Type.upper() == 'ALL':
            self.ALL = True
            self.FWHM = False
            self.lFWHM = False

        elif Type.upper() == 'FWHM':
            self.FWHM = True
            self.ALL = False
            self.lFWHM = False
        elif Type.upper() == 'LFWHM':
            self.lFWHM = True
            self.FWHM = False
            self.ALL = False
        elif Type.upper() == 'SINGLE':
            self.lFWHM = False
            self.FWHM = False
            self.ALL = False
        else: 
            messagebox.showerror("Error", f"Unknown type. Showing single camera")
            self.lFWHM = False
            self.FWHM = False
            self.ALL = False




        self.connect_button = tk.Button(top_panel, text="Connect", command=self.connect)
        self.connect_button.pack(side=tk.LEFT, padx=5)

        self.disconnect_button = tk.Button(top_panel, text="Disconnect", command=self.disconnect)
        self.disconnect_button.pack(side=tk.LEFT, padx=5)

        self.cool_button = tk.Button(top_panel, text="Cool to:", command=self.cool)
        self.cool_button.pack(side=tk.LEFT, padx=5)

        self.start_button = tk.Button(top_panel, text="Start", command=self.start_sequence)
        self.start_button.pack(side=tk.LEFT, padx=5)
        
        self.stop_button = tk.Button(top_panel, text="Stop", command=self.stop_sequence, state=tk.DISABLED)
        self.stop_button.pack(side=tk.LEFT, padx=5)

        self.take_bkg = tk.Button(top_panel, text="Take Background", command=self.take_bkg)
        self.take_bkg.pack(side="left", padx=5)
        
        self.save_img_as = tk.Button(top_panel, text="Save As", command=self.save_as)
        self.save_img_as.pack(side="left", padx=5)
        self.save_img_p = tk.Button(top_panel, text="Save", command=self.save)
        self.save_img_p.pack(side="left", padx=5)
        if self.FWHM:
            self.thresh_signi_button = tk.Button(left_panel, text="Set trace significance", command=self.set_thresh_signi)
            self.thresh_signi_button.pack(side="top", padx=5)

        self.exp_time_button = tk.Button(left_panel, text="Set Exposure Time", command=self.set_exposure_time)
        self.exp_time_button.pack(side="top", padx=5)
        if not self.ALL:        
            self.camera_index_button = tk.Button(left_panel, text="Set Camera Index", command=self.set_camera_index)
            self.camera_index_button.pack(side="top", padx=5)
        self.row_range_button = tk.Button(left_panel, text="Set Row Range", command=self.set_row_range)
        self.row_range_button.pack(side="top", padx=5)
        
        self.col_range_button = tk.Button(left_panel, text="Set Column Range", command=self.set_col_range)
        self.col_range_button.pack(side="top", padx=5)

        # Temperature Display
        self.temperature_label = tk.Label(master, text="Current Temperature: N/A", font=("Helvetica", 12))
        self.temperature_label.pack(pady=5)

        # Image Display Figure
        plot_panel = tk.Frame(master)
        plot_panel.pack(side="right", fill="both", expand=True)
        
        if self.FWHM:
            self.fig, (self.ax_image, self.ax_plot) = plt.subplots(2, 1, figsize=(6, 8))
            self.ax_plot.set_title("FWHM Plot")
        elif self.lFWHM:
            self.fig, (self.ax_image, self.ax_plot) = plt.subplots(2, 1, figsize=(6, 8))
            self.ax_plot.set_title("Line FWHM Plot")
        elif self.ALL: 
            self.fig, self.ax_grid = plt.subplots(2, 2, figsize=(8, 8))
        else: 
            self.fig, (self.ax_image) = plt.subplots(1, 1, figsize=(8, 8))
        if not self.ALL:
            self.ax_image.set_title("Real-Time Exposure")

        self.canvas = FigureCanvasTkAgg(self.fig, plot_panel)
        self.canvas_widget = self.canvas.get_tk_widget()
        self.canvas_widget.pack(fill="both", expand=True)
    
    def start_sequence(self):
        self.running = True
        self.start_button.config(state=tk.DISABLED)
        self.stop_button.config(state=tk.NORMAL)
        threading.Thread(target=self.capture_sequence, daemon=True).start()
    
    def stop_sequence(self):
        self.running = False
        self.start_button.config(state=tk.NORMAL)
        self.stop_button.config(state=tk.DISABLED)

    def connect(self):
        try:
            messagebox.showinfo('Connect','Connecting. This might take a moment')
            connect  = DS_obj.connect()
            if not connect:
                raise Exception
            else:
                self.connected  = True
                messagebox.showinfo('Success!','Connected to all cameras')

        except Exception as e:
            messagebox.showerror("Error", f"Failed to connect to all cameras. Try disconnecting and connecting again. If this doesn't work, restart the cameras and GUI")
    
    def disconnect(self):
        try:
            disconnect  = DS_obj.disconnect()
            if not disconnect.all():
                raise Exception
            else:
                self.connected  = False
                messagebox.showinfo('Success!','Diconnected from all cameras')

        except Exception as e:
            messagebox.showerror("Error", f"Failed to disconnect from all cameras.")
    def cool(self):
        try:
            new_T_goal = simpledialog.askinteger("T_goal", "Enter target detector temperature:", initialvalue=self.T_goal)
            if (new_T_goal is not None) and ((new_T_goal <40)&(new_T_goal>-100)):
                self.T_goal = new_T_goal
                out = DS_obj.cool(T_goal = new_T_goal)
                if out.all():
                    self.cooled = True
                    messagebox.showinfo('Success!',f'set T_goal to {new_T_goal} on all cameras')
                    self.update_temperature()
                else: 
                    raise Exception
            else: 
                raise Exception
        except Exception as e:
            messagebox.showerror("Error", f"failed to cool to T_goal. T should be between 40 and -100 C: {e}")
    
    def update_temperature(self):
       try:
           temp = DS_obj.get_temperature()
           string = 'Back/Front Temperature: '
           for i in range(len(DS_obj.bands)):
               band = DS_obj.bands[i]
               T_back = temp[i][0]
               T_front = temp[i][1]
               T_str = f' {band}: {T_back}C/{T_front}C |'
               string+=T_str
           self.temperature_label.config(text=string)
       except Exception as e:
           self.temperature_label.config(text="Current Temperature: Error")
       
    def save_img(self, type = 'AS'):
        try:
            if type != 'AS':
                new_save_path = self.save_loc
            else:
                new_save_path = filedialog.asksaveasfile()
                new_save_path = new_save_path.name
                self.save_loc = new_save_path
            
            ihdul = self.img_hdul
            bhdul = self.bkg_hdul
            isub_hdul = self.img_sub_hdul
            hdul = []
            phdu = fits.PrimaryHDU(None, header = None)
            hdr = {}
            hdul.append(phdu)
            ext = 1
            for hdu in ihdul: 
                hdul.append(hdu)
                hdr['EXT_'+str(ext)] = 'Raw Exsposure, {0}'.format(hdu.header['BAND'])
                ext+=1
            for hdu in bhdul: 
                hdul.append(hdu)
                hdr['EXT_'+str(ext)] = 'Bkg Exsposure, {0}'.format(hdu.header['BAND'])
                ext+=1            
            for hdu in isub_hdul: 
                hdul.append(hdu)
                hdr['EXT_'+str(ext)] = 'Bkg subtracted Exsposure, {0}'.format(hdu.header['BAND'])
                ext+=1 
            HDR = fits.Header()
            HDR['simple'] = 'T'
            HDR['BITPIX'] = 32
            HDR['NAXIS'] = 0

            for key in hdr.keys():
                HDR[key] = hdr[key]
            hdul[0].header = HDR
            hdul = fits.HDUList(hdus=hdul)
            hdul.writeto(new_save_path, overwrite=True)
        except Exception as e:
            messagebox.showerror("Error", f"Save file failed: {e}")
            import ipdb; ipdb.set_trace()
    def save(self):
        if self.save_loc == '':
            self.save_img(type = 'AS')
        else:
            self.save_img(type = 'Prev')
        pass 
    def save_as(self):
        self.save_img(type = 'AS')
        pass 

        messagebox.showinfo('Info','Saved file successfully')
        pass
    def set_exposure_time(self):
        try:
            new_exp = simpledialog.askfloat("Exposure Time", "Enter exposure time (s):", initialvalue=self.t_exp)
            if new_exp is not None and new_exp > 0:
                self.t_exp = new_exp
        except Exception as e:
            messagebox.showerror("Error", f"Invalid exposure time: {e}")
    
    def set_thresh_signi(self):
        try:
            new_thresh = simpledialog.askfloat("Significance Threshold", "Enter significance (fold over bkg):", initialvalue=self.thresh_signi)
            if new_thresh is not None and new_thresh > 0:
                self.new_thresh = new_thresh
        except Exception as e:
            messagebox.showerror("Error", f"Invalid significance: {e}")
    
    def take_bkg(self):
        try:
            bkg = DS_obj.expose(t_exp=self.t_exp)  # Capture the bkg image
            bhdul = []
            for i in range(len(DS_obj.bands)):
                phdu = bkg[i]
                bhdu = fits.ImageHDU(data=phdu.data.astype(float), header=phdu.header)
                bhdu.header['TYPE'] = 'BKG'
                bhdul.append(bhdu)
            bhdul = fits.HDUList(hdus=bhdul)
            self.bkg_hdul = bhdul
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            timestamp = 'bkg: '+timestamp
            if not self.ALL:
                data = bkg[self.ind].data.T
                data = data.astype(float)
                self.display_image(data, timestamp, ax = self.ax_image)
                self.update_temperature()  # Update temperature after exposure

            elif self.ALL:
                for i in range(len(DS_obj.bands)): 
                    ind = self.grid_ind[DS_obj.bands[i]]
                    ax = self.ax_grid[ind[0],ind[1]]
                    self.ax_grid[ind[0],ind[1]].set_title(DS_obj.bands[i])
                    data = bkg[i].data.T
                    data = data.astype(float)
                    self.display_image(data, timestamp,ax = ax)
                    self.update_temperature()  # Update temperature after exposure
                    band = DS_obj.bands[i]
                    ax.set_xlabel(f"{band}: {timestamp}")
        except Exception as e:
            messagebox.showerror("Error", f"failed to take bkg. Error msg: {e}")
            import ipdb; ipdb.set_trace()  
        pass
    def set_camera_index(self):
        try:
            new_index = simpledialog.askinteger("Camera Index", "Enter camera index (0-3):", initialvalue=self.ind)
            if new_index is not None and 0 <= new_index <= 3:
                self.ind = new_index
            else:
                messagebox.showerror("Error", "Camera index must be between 0 and 3.")
        except Exception as e:
            messagebox.showerror("Error", f"Invalid camera index: {e}")
        try: 
            if self.img_hdul is not None: 
                self.refresh_plot()
        except Exception as e: 
            messagebox.showerror("Error", f"Failed to refresh image: {e}")
        pass
    def set_row_range(self):
        try:
            start = simpledialog.askinteger("Row Range Start", "Enter start row:", initialvalue=self.row_range[0])
            end = simpledialog.askinteger("Row Range End", "Enter end row (None for no limit):", initialvalue=self.row_range[1])
            if start is not None and (end is None or end > start):
                self.row_range = (start, end)
            else:
                messagebox.showerror("Error", "Invalid row range.")
        except Exception as e:
            messagebox.showerror("Error", f"Invalid row range: {e}")
        try: 
            if self.img_hdul is not None: 
                self.refresh_plot()
        except Exception as e: 
            messagebox.showerror("Error", f"Failed to refresh image: {e}")
        pass

    def set_col_range(self):
        try:
            start = simpledialog.askinteger("Column Range Start", "Enter start column:", initialvalue=self.col_range[0])
            end = simpledialog.askinteger("Column Range End", "Enter end column (None for no limit):", initialvalue=self.col_range[1])
            if start is not None and (end is None or end > start):
                self.col_range = (start, end)
            else:
                messagebox.showerror("Error", "Invalid column range.")
        except Exception as e:
            messagebox.showerror("Error", f"Invalid column range: {e}")
        try: 
            if self.img_hdul is not None: 
                self.refresh_plot()
        except Exception as e: 
            messagebox.showerror("Error", f"Failed to refresh image: {e}")
        pass

    def capture_sequence(self):
        while self.running:
            try:
                img = DS_obj.expose(t_exp=self.t_exp)  # Capture the image
                ihdul = []
                for i in range(len(img)):
                    phdu = img[i]
                    ihdu = fits.ImageHDU(data=phdu.data.astype(float), header=phdu.header)
                    ihdu.header['TYPE'] = 'EXP'
                    ihdul.append(ihdu)
                ihdul = fits.HDUList(hdus=ihdul)
                self.img_hdul = ihdul
                if self.bkg_hdul is not None: 
                   isub_hdul = []
                   for i in range(len(img)):
                       ihdu = ihdul[i]
                       bhdu = self.bkg_hdul[i]
                       isub_hdu = fits.ImageHDU(data=ihdu.data - bhdu.data, header=phdu.header)
                       isub_hdu.header['TYPE'] = 'BKG SUB'
                       isub_hdul.append(isub_hdu)
                   isub_hdul = fits.HDUList(hdus=isub_hdul)
                   self.img_sub_hdul = isub_hdul
                self.refresh_plot()
                self.update_temperature()  # Update temperature after exposure

            except Exception as e:
                messagebox.showerror("Error", f"Failed to capture image: {e}")
                self.stop_sequence()
    def refresh_plot(self):
        if self.ind < len(self.img_hdul):
            timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
            if self.ALL:
                for i in range(len(DS_obj.bands)):
                    ind = self.grid_ind[DS_obj.bands[i]]
                    ax = self.ax_grid[ind[0],ind[1]]
                    data = self.img_hdul[i].data.T  # Transpose the image for FWHM computation
                    data = data.astype(float)
                    if self.bkg_hdul is not None: 
                        background = self.bkg_hdul[i].data.T
                    else: 
                        background = np.zeros_like(data)
                    self.display_image(data-background, timestamp, ax = ax)  # Display the selected camera's image
                    band = DS_obj.bands[i]
                    ax.set_xlabel(f"{band}: {timestamp}")


                    ind = self.grid_ind[DS_obj.bands[i]]
                    self.ax_grid[ind[0],ind[1]].set_title(DS_obj.bands[i])
            else:
                data = self.img_hdul[self.ind].data.T  # Transpose the image for FWHM computation
                data = data.astype(float)
                if self.bkg_hdul is not None: 
                    background = self.bkg_hdul[self.ind].data.T
                else: 
                    background = np.zeros_like(data)
                self.display_image(data-background, timestamp, ax = self.ax_image)  # Display the selected camera's image
        else:
            if not self.ALL:
                self.ax_image.clear()
                self.ax_image.set_title(f"Camera {self.ind} not available")
                self.canvas.draw()

    def display_image(self, image_data, timestamp, ax):
        cropped_data = image_data[
            self.row_range[0] : self.row_range[1],
            self.col_range[0] : self.col_range[1],
        ]
        ax.imshow(cropped_data)
        if not self.ALL:
            ax.set_title(f"Real-Time Exposure (Camera {self.ind})")
            ax.set_xlabel(f"Timestamp: {timestamp}")

        # Compute and plot FWHM
        mask = np.ones_like(cropped_data)  # Replace with actual mask if needed
        mask = mask.astype(bool)

        if self.FWHM:
            try:
                r = get_cont_FWHM(cropped_data, mask, thresh_signi=self.thresh_signi)
            except: 
                r = get_cont_FWHM(cropped_data.T, mask.T, thresh_signi=self.thresh_signi)
            self.ax_plot.clear()
            try:
                self.ax_plot.plot(r.T)
            except:
                print('failed to get FWHM')
            self.ax_plot.set_ylabel('FWHM [microns]')
        elif self.lFWHM:
            try: 
                self.ax_plot.clear()
                sp = get_1d_spec(cropped_data, N_traces=2, buffer = 5)
                for i in range(2):
                    s = sp[:,i]
                    try:                   
                        fwhm_valid, peaks_valid = get_res_fwhm(s, np.zeros_like(s) , half_width = 5,pix_size = 13)
                        self.ax_plot.plot(peaks_valid, fwhm_valid)
                    except:
                        print('failed to get FWHM')
                self.ax_plot.set_ylabel('Line FWHM [microns]')
            except Exception as e:
                messagebox.showerror("Error", f"Failed to plot lFWHM: {e}")                
        self.canvas.draw()




class WrapperGUI:
    def __init__(self, master):
        self.master = master
        self.master.title("Launch DeepSpec Exposure Viewer/Focusing Assitant")
        
        tk.Label(master, text="Select Viewer Type:", font=("Helvetica", 14)).pack(pady=10)

        # Buttons for each type
        self.fwhm_button = tk.Button(master, text="Continuum source analysis", command=lambda: self.launch_gui("FWHM"))
        self.fwhm_button.pack(pady=5)

        self.lfwhm_button = tk.Button(master, text="Line source analysis", command=lambda: self.launch_gui("LFWHM"))
        self.lfwhm_button.pack(pady=5)

        self.single_button = tk.Button(master, text="Single Camera Viewer", command=lambda: self.launch_gui("SINGLE"))
        self.single_button.pack(pady=5)

        self.all_button = tk.Button(master, text="All Cameras Viewer", command=lambda: self.launch_gui("ALL"))
        self.all_button.pack(pady=5)

    def launch_gui(self, gui_type):
        try:
            new_window = tk.Toplevel(self.master)
            RealTimeGUI(new_window, Type=gui_type)
        except Exception as e:
            messagebox.showerror("Error", f"Failed to launch GUI with type {gui_type}: {e}")



# Main Application
if __name__ == "__main__":
    DS_obj = DeepSpecAPI()
    print('starting-up....')
    DS_obj.start_up()
    DS_obj.cool(T_goal = -70)
    for cam in DS_obj.cameras:
        cam.set_readout_speed(3) #3 = 3MHz, 7 = 50Khz
    root = tk.Tk()
    WrapperGUI(root)
    root.mainloop()

