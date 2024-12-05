import scipy.signal as sig
from scipy.optimize import curve_fit
import time
import os
import numpy as np 




### analysis
def get_1d_spec(cimg, N_traces=2, buffer = 5):
    x = np.sum(cimg,axis = 1)
    x = x - np.median(x)
    peaks = sig.find_peaks(x)
    traces = peaks[0][np.argsort(x[peaks[0]])][-N_traces:]
    d1_spec = np.zeros([np.shape(cimg)[1],N_traces])
    for i,trace in enumerate(traces): 
        ind = list(range(int(trace - buffer), int(trace + buffer)))
        spec_2d = cimg[ind,:]
        s1d = spec_2d.sum(axis = 0)
        d1_spec[:,i] = s1d
    return d1_spec


def find_traces(img, N_traces=2):
    img = remove_bkg_2d(img)
    spat_sum = img.sum(axis = 1)
    spat_peaks,_ = sig.find_peaks(spat_sum)
    prominences = sig.peak_prominences(spat_sum, spat_peaks)[0]
    spat_peaks = spat_peaks[np.argsort(prominences)[-N_traces:]]
    return spat_peaks


def remove_bkg_2d(img):
    iqr = np.nanpercentile(img,75) - np.nanpercentile(img,25)+1
    med = np.nanmedian(img)
    med_bkg = np.nanmedian(img[img<med+3*iqr])
    img = img - med_bkg
    return img 


def get_cont_FWHM(img_in, mask, buffer = 30, N_traces = 2, thresh_signi = 15, pix_size_um=13, width_um = 100, res =10, M = 2.35/1.5 ):
    img = img_in.astype(float)
    spat_peaks = find_traces(img, N_traces=N_traces)
    L = np.shape(img)[1]
    FWHM = np.zeros((N_traces,L))
    FWHM_mask = np.zeros((N_traces,L))
    N_width = width_um//pix_size_um
    N_width = int(N_width//2)
    bkg_img = img.copy()
    for peak in spat_peaks:
        bkg_img[peak-buffer:peak+buffer,:] = np.nan
    bkg_mean = np.nanmean(bkg_img)
    bkg_std = np.nanstd(bkg_img)
    img_m = img.copy()
    img_m[mask] = np.nan
    for j,peak in enumerate(spat_peaks):
        cutout_2d = img[peak-buffer:peak+buffer,:]
        cutout_2d_mask = img_m[peak-buffer:peak+buffer,:]
        cutout_1d_mask = np.isnan(img_m[peak-buffer:peak+buffer,:].sum(axis = 0) )
        FWHM_mask[j,:] = cutout_1d_mask
        for i in range(L):
            slice_1d = cutout_2d[:,i-N_width:i+N_width].sum(axis = 1)
            m = np.max(slice_1d)
            if m>bkg_mean + thresh_signi*bkg_std:
                slice_1d = slice_1d/m
                grid = np.arange(len(slice_1d)*res)/res
                inter = np.interp(grid,np.arange(len(slice_1d)),slice_1d ) 
                fwhm = np.where(inter>0.5)[0][-1] - np.where(inter>0.5)[0][0]+1 
            else: 
                fwhm = np.nan
            FWHM[j,i] = fwhm
    FWHM = FWHM * pix_size_um/res/M
    return FWHM

def remove_bkg(spec_1d):
    iqr = np.percentile(spec_1d,75) - np.percentile(spec_1d,25)
    if iqr == 0: 
        iqr = iqr+0.1
    med = np.median(spec_1d)
    spec_1d_smo = spec_1d.copy()
    x, y = (np.arange(len(spec_1d_smo))[spec_1d_smo<med+1*iqr], spec_1d_smo[spec_1d_smo<med+1*iqr])
    pol, _ = polynomial_smooth(x,y,y_err=None,n=3)
    spec_1d_smo = spec_1d_smo.copy() - pol.eval(np.arange(len(spec_1d_smo)))
    iqr = np.percentile(spec_1d_smo,75) - np.percentile(spec_1d_smo,25)
    if iqr == 0: 
        iqr = iqr+0.1
    med = np.median(spec_1d_smo)
    x, y = (np.arange(len(spec_1d_smo))[spec_1d_smo<med+1*iqr], spec_1d_smo[spec_1d_smo<med+1*iqr])
    pol2, _ = polynomial_smooth(x,y,y_err=None,n=3)
    spec_1d_smo = spec_1d_smo.copy() - pol2.eval(np.arange(len(spec_1d_smo)))
    return spec_1d_smo

def find_lines(spec_1d, mask, signi = 4): 
    spec_1d_smo = remove_bkg(spec_1d)
    iqr = np.percentile(spec_1d_smo,75) - np.percentile(spec_1d_smo,25)
    med = np.median(spec_1d_smo)
    spec_1d_smo[spec_1d_smo<med+signi*iqr] = 0
    peaks,_ = sig.find_peaks(spec_1d_smo)
    #prominences = sig.peak_prominences(spec_1d_smo, peaks)[0]
    peak_mask =[]
    peaks = peaks[(peaks>10)&(peaks<len(spec_1d)-9)]
    ind = np.where(np.diff(peaks)<4)[0]
    omit= np.minimum(peaks[ind],peaks[ind+1])
    peaks = peaks[[x not in omit for x in peaks]]
    #import pdb; pdb.set_trace()
    for x in peaks: peak_mask.append(x in np.where(mask)[0])
    peak_mask = np.array(peak_mask)
    return peaks, peak_mask 

def FWHM(spec_1d, peaks, half_width = 5,pix_size = 13, N = 1000): 
    fwhm_list =[]
    spec_1d_smo = remove_bkg(spec_1d)
    for peak in peaks: 
        if peak<half_width:
            s_red = spec_1d_smo[0: peak +half_width].copy()
        else: 
            s_red = spec_1d_smo[peak -half_width: peak +half_width].copy()
        grid = np.arange(N*len(s_red))/N
        inter = np.interp(grid, np.arange(len(s_red)), s_red)
        #cum = np.cumsum(inter)/np.sum(inter)
        inter_max = inter/np.max(inter)
        fwhm = np.where(inter_max>0.5)[0][-1] - np.where(inter_max>0.5)[0][0]+1
        if peak<half_width:
            fwhm = np.where(inter_max>0.5)[0][-1] - np.where(inter_max==np.max(inter_max))[0][0]+0.5
        fwhm = pix_size*fwhm/N
        fwhm_list.append(fwhm)
    return np.array(fwhm_list)


class Polynomial:
    def __init__(self, *coefficients):
        """ input: coefficients are in the form a_n, ...a_1, a_0 
        """
        self.coefficients = list(coefficients) # tuple is turned into a list
     
    def __repr__(self):
        """
        method to return the canonical string representation 
        of a polynomial.
        """
        return "Polynomial" + str(tuple(self.coefficients))
    def eval(self, x):    
        res = 0
        for index, coeff in enumerate(self.coefficients[::-1]):
            res += coeff * x** index
        return res 
    def poly_func(self, t,*coeffs):
        self.coefficients=coeffs
        M=self.eval(t)
        return M

def polynomial_smooth(x,y,y_err=None,n=5):
    coeffs_init=np.ones_like(np.arange(0,n))
    poly=Polynomial(*coeffs_init)
    popt, _ = curve_fit(poly.poly_func, x, y, sigma=y_err,p0=coeffs_init)
    return poly, popt     


def get_res_fwhm(sp, masks , half_width = 5,pix_size = 13):
    peaks, masks_img  = find_lines(sp, masks)
    fwhm = FWHM(sp, peaks, half_width = half_width,pix_size = pix_size)
    fwhm_valid = fwhm[(~masks_img)]
    peaks_valid = peaks[(~masks_img)]
    return fwhm_valid, peaks_valid


