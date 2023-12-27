# -*- coding: utf-8 -*-
"""
Exposure time calculator for DeepSpec
@author: idoi
"""

from astropy import constants
import numpy as np
import sys
#UTILS_PATH=os.environ['UTILS_PATH']
UTILS_PATH='C:\\Users\\idoi\\Dropbox\\Utils'
sys.path.insert(1,UTILS_PATH)
import matplotlib.pyplot as plt
import math
from scipy.optimize  import newton_krylov



#defualt parameters
R=600 # Resolution
T_exp=1500 # Max exposure time
wl=5000 # Wavelength for SNR calculation
d_lam=wl/R # Resolution element width in Ang d_lam = R/lambda
r_cm=61/2 #cm telescope radius in cm 
mu=13.0 # pixel size in microns
F_num=3.0   # F number
f_mm=2*F_num*r_cm*10 # focal length in mm
Sky_brightness_surface_den=20.5 #mag/arcsec**2 sky background surface brightness
bkg2src=1 # ratio of number of background to source pixels used for estimating the background in a single resolution element
n_tel=20 # number of telescopes 
fiber=25*235/150 # fiber image size on the focal plane, in microns
DC=0.0008 # Dark current  e-/s/pixel @ -80C for GreatEyes ELSEi
read_noise=2.8 # Read noise e-/pixel @ 50kHz pixels/sec for GreatEyes ELSEi
cosmetic_effects=0.005 #not implemented


Eff_spec= ((0.985)**4)*((0.99)**6)*((0.995)**6)*0.9 #(0.985)**number of surfaces in col* (0.99)**number of surfaces in cam* (0.995)**number of bonding layers* CCD QE *prism throughput
Eff_spec_tel= 0.8*Eff_spec # including reflection losses at the telescope (assumed to be 20%)
Eff= 0.83*Eff_spec_tel # including fiber throughput (measured at 83%)



#Theoretical resolution and throughput for prisms as a function of wavelength
Respath1 =       'C:\\Users\\idoi\\Dropbox\\DeepSpec\\single_prism\\350nm_440nm_F_SILICA_data.txt'
Respath2 =       'C:\\Users\\idoi\\Dropbox\\DeepSpec\\single_prism\\438nm_545nm_PBL35Y_data.txt'
Respath3 =       'C:\\Users\\idoi\\Dropbox\\DeepSpec\\single_prism\\545nm_680nm_L-LAH85V_data.txt'
Respath4 =       'C:\\Users\\idoi\\Dropbox\\DeepSpec\\single_prism\\680nm_900nm_S-LAH99_data.txt'

# not implemented yet: throughput for QE, reflection losses, fiber throughput, and lens losses as a function of wavelength
# change prisms to real prisms used in DS design


Res_dic = {}
wl_dic = {}
Trans_dic = {}

for i in range(4):
    Respath = eval('Respath'+'{0}'.format(i+1))
    Wlength_m, Res, Trans = np.genfromtxt(Respath)[0:3]
    Wl_Ang = Wlength_m*1e10
    Res_dic[i] =Res
    Trans_dic[i] = Trans
    wl_dic[i] = Wl_Ang

#stitch together Resolution and trans values for each prism by wavelength 
Res = np.hstack((Res_dic[0],Res_dic[1],Res_dic[2],Res_dic[3]))
Wave = np.hstack((wl_dic[0],wl_dic[1],wl_dic[2],wl_dic[3]))
Trans = np.hstack((Trans_dic[0],Trans_dic[1],Trans_dic[2],Trans_dic[3]))

#SNR_vec = SNR_wl_array(Wave,20,Res, 900, eff = Trans*Eff)




def Sky_bkg_density_per_pixel(mu,Sky_brightness_surface_den,f_mm=1000):
    '''
    Sky background surface brightness in mag/arcsec**2
    In: mu - pixel size in microns
    In: Sky_brightness_surface_den - sky background surface brightness in mag/arcsec**2
    In: f_mm - focal length in mm
    Out: Sky_brightness_mag - sky background surface brightness in mag/pixel
    '''
    Plate_Scale=206265*mu/1000/f_mm #arcsec/pixela
    Sky_brightness_flux_den=10**(-0.4*(Sky_brightness_surface_den))
    Sky_brightness_flux=Sky_brightness_flux_den*Plate_Scale**2
    Sky_brightness_mag=-2.5*np.log10(Sky_brightness_flux)
    return Sky_brightness_mag

def N_phot_in_element(M_ab,T_exp=500,d_lam=d_lam,r_cm=r_cm,wl_AA=wl):
    '''
    Number of photons in a single resolution element
    In: M_ab - source magnitude
    In: T_exp - exposure time
    In: d_lam - resolution element width in cm
    In: r_cm - telescope radius in cm
    In: wl_AA - wavelength in Angstrom
    Out: N_phot - number of photons in a single resolution element    
    '''
    h=constants.h.cgs.value
    c=constants.c.cgs.value
    wl_cm=wl_AA*1e-8
    nu_cgs=c/wl_cm
    photon_energy=h*nu_cgs
    f_lam=mag2flux(M_ab,wl_AA)
    luminosity=f_lam*d_lam*np.pi*r_cm**2
    photon_flux=luminosity/photon_energy
    N_phot=photon_flux*T_exp
    return N_phot

def SNR(M_src,T_exp=500,Sky_brightness= Sky_brightness_surface_den,seeing=2,DC=DC,Read=read_noise,G=1,Sig_ADU=0.289,eff=Eff,n_tel=n_tel,mu=mu,fiber=fiber,d_lam=d_lam,r_cm=r_cm,bkg2src=bkg2src,binning=[1,1],wl_AA=wl):
    '''
    Signal to noise ratio for a single resolution element, using the ccd equation (Merline & Howell 1995)
    In: M_src - source magnitude
    In: T_exp - exposure time, default is 500s
    In: Sky_brightness - sky background surface brightness in mag/arcsec**2, default is 20.5 mag/arcsec**2
    In: seeing - seeing in arcsec  #not implemented yet as slit losses
    In: DC - dark current in e-/s  default is 0.0008 e-/s for GreatEyes ELSEi at -80C
    In: Read - read noise in e-    default is 2.8 e- for GreatEyes ELSEi at 50kHz pixels/sec
    In: G - gain in e-/ADU         default is 1 e-/ADU for GreatEyes ELSEi at 50kHz pixels/sec
    In: Sig_ADU - rounding error   default is standard 0.289 
    In: eff - total efficiency     Default is  a conservative 0.51, assming 80% fiber and telescope throughput, 83% CCD QE, and accounting for prism, lens throughput and reflection losses
    In: n_tel - number of telescopes #default at 20 telescopes
    In: mu - pixel size in microns #default at 13 micron pixel
    In: fiber - fiber diameter in microns #default at 25 micron fiber
    In: d_lam - resolution element width in Ang d_lam = R/lambda #default is calculated for R=600 at 5000A 
    In: r_cm - telescope radius in cm # default at 30.5cm
    In: bkg2src - ratio of number of background to source pixels used for estimating the background in a single resolution element #default at 1. For longslit instruments this is typically -> infinity
    In: binning - binning in x and y directions #default at [1,1]
    In: wl_AA - wavelength in Angstrom #default at 5000A
    Out: SNR - signal to noise ratio for a single resolution element
    Out: noise_comp - variance components in the ccd equation for source, background, dark current and readout noise respectively. 
    '''
    w=math.ceil(fiber/mu)
    if w<3:
        w=3
    fill_factor=1 #ratio of used to unused pixels which are read out

    N_pix_tot=n_tel*w**2*(1/fill_factor)

    bin_factor=binning[0]*binning[1]
    N_pix_read = N_pix_tot/bin_factor
    N2_readout_Gain=(Read**2+G**2*Sig_ADU**2)*N_pix_read
    N2_Dark=N_pix_tot*DC*T_exp

    pixel_sky_mag=Sky_bkg_density_per_pixel(mu,Sky_brightness,f_mm=f_mm)
    N_bkg=eff*n_tel*N_phot_in_element(pixel_sky_mag,T_exp=T_exp,d_lam=d_lam,r_cm=r_cm,wl_AA=wl_AA)
    N_src=eff*n_tel*N_phot_in_element(M_src,T_exp=T_exp,d_lam=d_lam,r_cm=r_cm,wl_AA=wl_AA)

    noise=np.sqrt(N_src+(1+1/bkg2src)*(N_bkg+N2_Dark+N2_readout_Gain))
    SNR=N_src/noise

    noise_comp=[N_src,(1+1/bkg2src)*N_bkg, (1+1/bkg2src)*N2_Dark, (1+1/bkg2src)*N2_readout_Gain]
    return SNR, noise_comp

def SNR_wl_array(M_src, Wlength_AA = Wave,Resolution=Res,T_exp = 500,f_SNR = SNR,Sky_brightness= Sky_brightness_surface_den,seeing=2,DC=DC,Read=read_noise,G=2,Sig_ADU=0.289,eff=float(Eff),n_tel=n_tel,mu = mu,fiber=fiber,r_cm=r_cm,bkg2src=bkg2src,binning=[1,1]):
    '''
    Signal to noise ratio for a resolution element as a function of wavelength, using the ccd equation (Merline & Howell 1995)
    In: M_src - source magnitude
    In: Wlength_AA - wavelength array in Ang. Default is an array from 3500A to 9000A 
    In: Resolution - resolution array as a function of wavelength. Defualt is to use the predicted resolution for DeepSpec prisms
    In: T_exp - exposure time, default is 500s
    In: Sky_brightness - sky background surface brightness in mag/arcsec**2, default is 20.5 mag/arcsec**2
    In: seeing - seeing in arcsec  #not implemented yet as slit losses
    In: DC - dark current in e-/s  default is 0.0008 e-/s for GreatEyes ELSEi at -80C
    In: Read - read noise in e-    default is 2.8 e- for GreatEyes ELSEi at 50kHz pixels/sec
    In: G - gain in e-/ADU         default is 1 e-/ADU for GreatEyes ELSEi at 50kHz pixels/sec
    In: Sig_ADU - rounding error   default is standard 0.289
    In: eff - total efficiency     Default is  a conservative 0.51, assming 80% fiber and telescope throughput, 83% CCD QE, and accounting for prism, lens throughput and reflection losses
    In: n_tel - number of telescopes #default at 20 telescopes
    In: mu - pixel size in microns #default at 13 micron pixel
    In: fiber - fiber diameter in microns #default at 25 micron fiber
    In: r_cm - telescope radius in cm # default at 30.5cm
    In: bkg2src - ratio of number of background to source pixels used for estimating the background in a single resolution element #default at 1. For longslit instruments this is typically -> infinity
    In: binning - binning in x and y directions #default at [1,1]
    Out: SNR - signal to noise ratio for a single resolution element
    Out: noise_comp - variance components in the ccd equation for source, background, dark current and readout noise respectively.
    '''
    SNR = []
    noise_comp = []

    D_lam_vec = Wlength_AA/Resolution
    for i in range(len(Wlength_AA)):
        if type(eff) == float:
            Ef = eff
        elif len(eff)>1:
            Ef = eff[i]
        SNR_wl,n_comp = f_SNR(M_src,T_exp=T_exp,Sky_brightness=Sky_brightness,seeing=seeing,DC=DC,Read=read_noise,G=G,Sig_ADU=Sig_ADU,eff=Ef,
                    n_tel=n_tel,mu=mu,fiber=fiber,d_lam=D_lam_vec[i],r_cm=r_cm,bkg2src=bkg2src,binning=binning,wl_AA=Wlength_AA[i])
        SNR.append(SNR_wl)
        noise_comp.append(n_comp)


    SNR = np.array(SNR)
    noise_comp = np.array(noise_comp)
    return SNR
def find_limiting_mag(T_vec,f_SNR,SNR_lim,n_tel=n_tel,mu = mu,fiber=25,Sky_brightness= Sky_brightness_surface_den,Xin=[15,19],binning=[2,2],wl_AA=wl, d_lam=d_lam, read_noise = read_noise ):
    '''
    Find the limiting magnitude for a given exposure time and SNR limit
    In: T_vec - array of exposure times
    In: f_SNR - SNR function. default is SNR
    In: SNR_lim - SNR limit for the limiting magnitude. default is 10 sigma 
    In: n_tel - number of telescopes #default at 20 telescopes. default at 20 telescopes
    In: mu - pixel size in microns #default at 13 micron pixel
    In: fiber - fiber diameter in microns #default at 25 micron fiber
    In: Sky_brightness - sky background surface brightness in mag/arcsec**2, default is 20.5 mag/arcsec**2
    In: Xin - initial guess for the range containing the limiting magnitude. default is [15,19]
    In: binning - binning in x and y directions #default at [1,1]
    In: wl_AA - wavelength in Angstrom #default at 5000A
    In: d_lam - resolution element width in Ang d_lam = R/lambda #default is calculated for R=600 at 5000A
    In: read_noise - read noise in e-    default is 2.8 e- for GreatEyes ELSEi at 50kHz pixels/sec
    Out: sol - limiting magnitude
    Out: noise_comp - variance components in the ccd equation for source, background, dark current and readout noise respectively.
    '''
    def func(M,T):
        res=f_SNR(M_src=M,T_exp=T,n_tel=n_tel,mu=mu,fiber=fiber,Sky_brightness=Sky_brightness,binning=binning,wl_AA=wl_AA, d_lam=d_lam,  Read = read_noise)[0]-SNR_lim
        return res
    
    if (isinstance(T_vec,float)|(isinstance(T_vec,int))):
        f=lambda m: func(m,T=T_vec)
        if T_vec>50:
            xin=Xin[1]

        elif T_vec<50:
            xin=Xin[0]
        try:
            try:
                sol=newton_krylov(f,xin=xin,f_tol=1e-2)
            except:
                sol=newton_krylov(f,xin=xin-3,f_tol=1e-2)
        except:
            sol=newton_krylov(f,xin=xin+3,f_tol=1e-2)
        _ , noise_comp =f_SNR(M_src=sol,T_exp=T_vec,n_tel=n_tel,mu=mu,fiber=fiber,Sky_brightness=Sky_brightness,binning=binning,wl_AA=wl_AA, d_lam=d_lam,  Read = read_noise)
    else:
        sol=[]
        noise_comp=np.empty_like(np.array([0,0,0,0]))
        for T in T_vec:
            f=lambda m: func(m,T=T)
            if T>50:
                xin=Xin[1]
            elif T<50:
                xin=Xin[0]
            try:
                try:
                    res=newton_krylov(f,xin=xin,f_tol=1e-2)
                except:
                    res=newton_krylov(f,xin=xin-3,f_tol=1e-2)
            except:
                res=newton_krylov(f,xin=xin+3,f_tol=1e-2)
            sol.append(res)
            last=np.array(f_SNR(M_src=sol[0],T_exp=T,n_tel=n_tel,mu=mu,fiber=fiber,Sky_brightness=Sky_brightness,binning=binning,wl_AA=wl_AA, d_lam=d_lam,  Read = read_noise)[1])
            noise_comp=np.vstack([noise_comp,last])
        sol=np.array(sol)

    noise_comp=np.delete(noise_comp,0,0)
    return sol, noise_comp

def zp_piv_wl_ab(piv_wl):
    '''
    Zero point magnitude for a given wavelength in AB system
    '''
    # in angstrom
    c=constants.c.value*1e10
    zp_band_ab=-2.5*np.log10((piv_wl**2)/c)-48.6
    return  zp_band_ab
def mag2flux(mag,piv_wl):
    '''
    Convert magnitude to flux (f_lambda) at a given wavelength
    '''
    zp=zp_piv_wl_ab(piv_wl)
    flux=10**(-0.4*(mag-zp)) 
    return flux
def flux2mag(flux,piv_wl):
    '''
    Convert flux (f_lambda) to magnitude at a given wavelength
    '''
    zp=zp_piv_wl_ab(piv_wl)
    mag=-2.5*np.log10(flux)+zp 
    return mag
def magerr2fluxerr(magerr,flux):
    '''
    Convert flux error (f_lambda) to magnitude error at a given wavelength
    '''
    fluxerr=np.abs(-2.303/2.5*magerr*flux)
    return fluxerr
def fluxerr2magerr(fluxerr,flux):
    '''
    Convert magnitude error to flux error (f_lambda) at a given wavelength
    '''
    magerr=np.abs(-2.5/2.303*fluxerr/flux)
    return magerr
if __name__ == '__main__':
    log_exp=np.log10(T_exp)
    T_exp_vec=10**np.linspace(0,log_exp,100)
    
    SNR_readout={}
    
    color={
        '16':'m',
        '17':'c',
        '18':'y',
        '19':'b',
        '20':'g',
        '21':'r',
        '22':'k',
    }
    
    limmag_readout_10_sig,_    =find_limiting_mag(T_exp_vec,SNR,SNR_lim=10,n_tel=n_tel,mu=mu,fiber=fiber,Sky_brightness=Sky_brightness_surface_den,Xin=[14,17],binning=[2,2])
    limmag_readout_10_sig_single,_    =find_limiting_mag(T_exp_vec,SNR,SNR_lim=10,n_tel=1,mu=mu,fiber=fiber,Sky_brightness=Sky_brightness_surface_den,Xin=[14,17],binning=[2,2])
    limmag_readout_10_sig_4,_    =find_limiting_mag(T_exp_vec,SNR,SNR_lim=10,n_tel=4,mu=mu,fiber=fiber,Sky_brightness=Sky_brightness_surface_den,Xin=[14,17],binning=[2,2])

    plt.figure(figsize=[7,7])
    plt.plot([T_exp,T_exp],[10,22],'k--',alpha=0.25)
    plt.plot([300,300],[10,22],'k--',alpha=0.25)
    
    plt.plot([1,T_exp+200],[19.5,19.5],'k--',alpha=0.25)
    
    plt.plot(T_exp_vec,limmag_readout_10_sig,'r')
    plt.plot(T_exp_vec,limmag_readout_10_sig_single,'r--')
    plt.plot(T_exp_vec,limmag_readout_10_sig_4,'r.-')
        
    plt.title('Limiting magnitude for R={0}'.format(R),fontsize=15)
    plt.ylabel('10-sigma limiting magnitude (AB)',fontsize=13)
    plt.xlabel('Exposure time (s)',fontsize=13)
    plt.xscale('log')
    plt.ylim((10,22))
    plt.xlim((1,T_exp+200))
    plt.gca().invert_yaxis()
    plt.gca().tick_params(labelsize=12)
    plt.text(11,19.35,'19.5 mag')
    plt.text(1300,15,'1500 s',rotation = 90)
    plt.text(260,15,'300 s',rotation = 90)
    
    plt.text(73,np.mean(limmag_readout_10_sig)+1 - 1.75,'single trace',rotation = -35, color = '#990000',fontsize = 14)
    plt.text(73,np.mean(limmag_readout_10_sig)+1.45 - 1.25,'4 traces',rotation = -35, color = '#990000',fontsize = 14)
    plt.text(73,np.mean(limmag_readout_10_sig)+1.7 - 5/8,'20 traces',rotation = -35, color = '#990000',fontsize = 14)
    
    plt.text(260,15,'300 s',rotation = 90)
    plt.legend()
    plt.xlim(10,1600)
    plt.ylim(22,14.5)
    plt.tight_layout()
    plt.show()

