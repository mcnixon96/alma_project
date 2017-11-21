import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling import functional_models
import random

class RightAscension:
    def __init__(self,angle):
        self.hours=int(angle/(np.pi/12))
        self.minutes=int((angle-self.hours*(np.pi/12))/(np.pi/720))
        self.seconds=(angle-self.hours*(np.pi/12)-self.minutes*(np.pi/720))/(np.pi/43200)

class Declination:
    def __init__(self,angle):
        self.degrees=int(angle/(np.pi/180))
        self.minutes=int((angle-self.degrees*(np.pi/180))/(np.pi/(180*60)))
        self.seconds=(angle-self.degrees*(np.pi/180)-self.minutes*(np.pi/(180*60)))/(np.pi/(180*60*60))

        
def disk_distance_plot(band_no):

    if band_no==3:
        T_sys=67 #Band 3 parameters
        freq_lower=8.4*10**10
        freq_upper=1.16*10**11
        aperture_eff=0.69
        ref_flux=2.6*10**-28 #band 3 flux for HD 100546 in Wm-2Hz-1 from SED
    elif band_no==7:
        T_sys=251
        freq_lower=2.75*10**11
        freq_upper=3.7*10**11
        aperture_eff=0.74
        ref_flux=8*10**-24 #band 7 flux for HD 100546 in Wm-2Hz-1 from SED

    obs_time=3600 #observing time in seconds

    sigma=T_sys/((freq_upper-freq_lower)*obs_time)**0.5 #rms noise

    ref_temp=aperture_eff*12**2*ref_flux*10**26/3514.0

    ref_snr=ref_temp/sigma

    distance=np.logspace(2,6,50)
    snr=np.zeros(len(distance))

    for d in range(len(distance)):
        snr[d]=ref_snr*10**4/float(distance[d])**2

    plt.plot(distance,snr)
    plt.show()
    plt.ylim(0,10)
    plt.xscale('log')
    plt.xlabel('distance / pc')
    plt.ylabel('signal to noise ratio / sigma')

def disk_contour(distance,centre_ra=0,centre_dec=0):
    
    distance=float(distance)
    ra=np.arange(centre_ra-(8.73*10**-8*100),centre_ra+(8.73*10**-8*100),8.73*10**-8)
    dec=np.arange(centre_dec-(8.73*10**-8*100),centre_dec+(8.73*10**-8*100),8.73*10**-8)

    xx,yy=np.meshgrid(ra,dec)

    peak_flux=800*10**4/float(distance**2) #peak flux in mJy as function of distance (using band 7 here)

    fwhm=50. #in au
    fwhm_angle=abs(np.arctan((fwhm*4.848*10**-6)/distance))
    disk_ra=ra[int(len(ra)/2)]
    disk_dec=dec[int(len(dec)/2)]
    
    G=functional_models.Gaussian2D(amplitude=peak_flux,x_mean=disk_ra,y_mean=disk_dec,x_stddev=fwhm_angle/2.355,y_stddev=fwhm_angle/2.355)

    #Creating map of Extragalactic Background Light

    flux_bin_limits=[0.015,0.027,0.047,0.084,0.15,0.267,0.474,0.843,1.124]
    log_bin_count=np.array([5.5,5.6,5.3,5,5.1,4.9,4.7,4.2,3.4])
    bin_count=10**log_bin_count
    flux_bins=np.zeros((len(flux_bin_limits)-1,10))
    for i in range(len(flux_bin_limits)-1):
        flux_bins[i,:]=np.logspace(np.log10(flux_bin_limits[i]),np.log10(flux_bin_limits[i+1]),10)

    area=(np.amax(ra)-np.amin(ra))*(np.amax(dec)-np.amin(dec))
    fractional_area=area/(3.046*10**-4)
    adjusted_bin_count=fractional_area*bin_count
    
    mean_counts=np.random.poisson(adjusted_bin_count)

    for i in range(len(mean_counts)):
        if mean_counts[i]!=0:
            for j in range(mean_counts[i]):
                source_ra=random.choice(ra)
                source_dec=random.choice(dec)
                source_flux=random.choice(flux_bins[i,:])
                source=functional_models.Gaussian2D(amplitude=peak_flux,x_mean=source_ra,y_mean=source_dec,x_stddev=(8.73*10**-8)/2.355,y_stddev=(8.73*10**-8)/2.355)
                G=G+source
    
    levels=np.linspace(0,peak_flux,20)

    ra_units=[]
    ra_array=np.zeros(len(ra))
    dec_units=[]
    dec_array=np.zeros(len(dec))
    for i in range(len(ra)):
        ra_units.append(RightAscension(ra[i]))
        dec_units.append(Declination(dec[i]))
        ra_array[i]=ra_units[i].seconds
        dec_array[i]=dec_units[i].seconds


    z=G(xx,yy)
    fig=plt.figure(figsize=(8,6))
    ax1=fig.add_subplot(111)
    im=plt.contourf(ra,dec,z,levels)
    ax2=ax1.twinx()
    ax2.plot(dec,ra_array,alpha=0.0)
    ax3=ax1.twiny()
    ax3.plot(dec_array,ra,alpha=0.0)

    ax1.set_xlim(np.amin(ra),np.amax(ra))
    ax1.set_ylim(np.amin(dec),np.amax(dec))
    ax1.set_xlabel('Right ascension',labelpad=30)
    ax1.set_ylabel('Declination',labelpad=30)
    ax1.tick_params('both',labelbottom=0,labeltop=0,labelleft=0,labelright=0)
    ax2.tick_params('y',labelleft='on',labelright='off')
    ax3.tick_params('x',labeltop='off',labelbottom='on')
    fig.colorbar(im,ax=[ax1,ax2,ax3],pad=0.1)
    plt.show()

    return peak_flux
