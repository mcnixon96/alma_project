import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FormatStrFormatter
from astropy.modeling import functional_models
from scipy import signal, integrate
import random
from astropy.convolution import Gaussian2DKernel, convolve
import matplotlib.colors as clr
from astropy.io import fits

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

def disk_distance_plot(band_nos):

    for band_no in band_nos:
        if band_no==3:
            T_sys=60.0 #Band 3 parameters
            freq_lower=8.4*10**10
            freq_upper=1.16*10**11
            aperture_eff=0.69
            ref_flux=2.6*10**-28 #band 3 flux for HD 100546 in Wm-2Hz-1 from SED
        elif band_no==7:
            T_sys=219.0
            freq_lower=2.75*10**11
            freq_upper=3.7*10**11
            aperture_eff=0.74
            ref_flux=8*10**-27 #band 7 flux for HD 100546 in Wm-2Hz-1 from SED

        obs_time=3600 #observing time in seconds

        sigma=T_sys/((freq_upper-freq_lower)*obs_time)**0.5 #rms noise

        ref_temp=aperture_eff*12**2*ref_flux*10**26/3514.0

        ref_snr=ref_temp/sigma
        print ref_snr*10**4, band_no

        distance=np.logspace(0,6,50)
        snr=np.zeros(len(distance))

        for d in range(len(distance)):
            snr[d]=ref_snr*10**4/float(distance[d])**2

        plt.plot(distance,snr,linewidth=2)
        plt.xlim(1,np.amax(distance))
        plt.ylim(0,np.amax(snr))
        plt.xscale('log')
        plt.yscale('log')
        plt.xlabel('distance / pc')
        plt.ylabel('signal to noise ratio / sigma')
    plt.title('SNR of HD 100546 at different distances')
    plt.plot(distance,5*np.ones(len(distance)),'r')
    plt.legend(['Band 3','Band 7'])
    plt.savefig('HD_SNR.png')
    plt.show()

def disk_time_plot(distance):

    """distance in parsecs"""

    T_sys=219.0 #working in band 7 for now
    freq_lower=2.75*10**11
    freq_upper=3.7*10**11
    aperture_eff=0.74
    ref_flux=8*10**-27 #band 7 flux for HD 100546 in Wm-2Hz-1 from SED

    flux=ref_flux*10**4/(distance**2)
    print flux

    obs_time=np.linspace(1,20,20)
    snr=np.zeros_like(obs_time)
    
    for i in range(len(obs_time)):
        sigma=T_sys/((freq_upper-freq_lower)*3600*obs_time[i])**0.5 #rms noise
        temp=aperture_eff*12**2*flux*10**26/3514.0
        snr[i]=temp/sigma

    plt.plot(obs_time,snr,linewidth=2)
    plt.xlim(1,20)
    plt.xlabel('observing time / hours')
    plt.ylabel('signal-to-noise ratio')
    plt.title('SNR of disc at '+str(int(distance))+'pc with varying observing time')
    plt.plot(obs_time,5*np.ones(len(obs_time)),'r')
    plt.savefig('HD_time.png')
    plt.show()

def time_dist_plot():

    T_sys=219.0 #working in band 7 for now
    freq_lower=2.75*10**11
    freq_upper=3.7*10**11
    aperture_eff=0.74
    ref_flux=8*10**-1 #band 7 flux for HD 100546 in Wm-2Hz-1 from SED

    obs_time=np.linspace(1,50,50)
    distance=np.zeros_like(obs_time)
    
    for i in range(len(obs_time)):
        sigma=T_sys/((freq_upper-freq_lower)*3600*obs_time[i])**0.5 #rms noise
        temp=5*sigma
        flux=3514.0*temp/(aperture_eff*144) #in Jy
        distance[i]=(ref_flux*10**4/flux)**0.5

    plt.plot(obs_time,distance/1000,linewidth=2)
    plt.xlim(1,50)
    plt.xlabel('observing time / hours')
    plt.ylabel('distance / kpc')
    plt.title('Maximum distance for 5 sigma detection for given observing time')
    plt.savefig('time_dist.png')
    plt.show()
    
def measure_sigma(band_no):

    if band_no==3:
        T_sys=60.0 #Band 3 parameters
        freq_lower=8.4*10**10
        freq_upper=1.16*10**11
        aperture_eff=0.69
        ref_flux=2.6*10**(-2) #band 3 flux for HD 100546 in Jy from SED
    elif band_no==7:
        T_sys=219.0
        freq_lower=2.75*10**11
        freq_upper=3.7*10**11
        aperture_eff=0.74
        ref_flux=8*10**-1 #band 7 flux for HD 100546 in Jy from SED

    sigma_1hr=T_sys/((freq_upper-freq_lower)*3600)**0.5 #rms noise in 1hr
    sigma_10hr=T_sys/((freq_upper-freq_lower)*36000)**0.5 #rms noise in 10hr

    sigma_1hr_Jy=3514*sigma_1hr/(aperture_eff*144)    #sigma in Jy
    sigma_10hr_Jy=3514*sigma_10hr/(aperture_eff*144)
    print sigma_1hr_Jy, sigma_10hr_Jy

    distance=np.logspace(2,6,50)
    flux=np.zeros(len(distance))

    for d in range(len(distance)):
        flux[d]=ref_flux*10**4/float(distance[d])**2

    plt.plot(distance,flux,linewidth=2)
    plt.plot(distance,sigma_1hr_Jy*np.ones(len(distance)))
    plt.plot(distance,sigma_10hr_Jy*np.ones(len(distance)))

    plt.ylim(0,2*np.amax(np.array([sigma_1hr_Jy,sigma_10hr_Jy])))
    plt.legend(['Source flux','Sensitivity (1hr)','Sensitivity(10hr)'],loc=0)
    plt.xlabel('distance / pc')
    plt.ylabel('Flux / Jy')
    plt.xscale('log')
    plt.title('Flux versus distance for HD 100546 in Band '+str(band_no))
    plt.show()


def disk_contour(distance,centre_ra=3.0256,centre_dec=-1.22513,convolve='true'):
    
    #Note: for HD 100546 use ra=3.0256, dec=-1.22513

    T_sys=219.0
    freq_lower=2.75*10**11 #in Hz
    freq_upper=3.7*10**11
    aperture_eff=0.74
    dishsize=12 #dish diameter in metres
    baseline=1600 #max baseline in metres (15-16000)

    resolution=0.299615052/(baseline*((0.5*(freq_upper-freq_lower))/(10**9))) #in radians
    field_size=3.2433*10**7/(0.5*(freq_upper-freq_lower))

    distance=float(distance)

    #this one works

    step_no=np.ceil(field_size/resolution)
    ra=np.linspace(centre_ra-(field_size/2),centre_ra+(field_size/2),step_no)
    dec=np.linspace(centre_dec-(field_size/2),centre_dec+(field_size/2),step_no)

    xx,yy=np.meshgrid(ra,dec)

    peak_flux=800*10/float(distance**2) #peak flux in Jy as function of distance (using band 7 here)

    fwhm=100. #in au
    fwhm_angle=abs(np.arctan((fwhm*4.848*10**-6)/distance))
    disk_ra=ra[int(len(ra)/2)]
    disk_dec=dec[int(len(dec)/2)]
    
    G=functional_models.Gaussian2D(amplitude=peak_flux,x_mean=disk_ra,y_mean=disk_dec,x_stddev=fwhm_angle/2.355,y_stddev=fwhm_angle/2.355)
    #to put the disk in the middle set x_mean=disk_ra,y_mean=disk_dec
    #to randomly position the disk set x_mean=random.choice(ra), y_mean=random.choice(dec)

    z=G(xx,yy)
    
    #Creating map of Extragalactic Background Light

    flux_bin_limits=10**-3*np.array([0.015,0.027,0.047,0.084,0.15,0.267,0.474,0.843,1.124])
    log_bin_count=np.array([5.5,5.6,5.3,5.1,4.9,4.7,4.2,3.4])
    bin_count=10**log_bin_count
    flux_bins=np.zeros((len(bin_count),10))
    for i in range(len(flux_bin_limits)-1):
        flux_bins[i,:]=np.logspace(np.log10(flux_bin_limits[i]),np.log10(flux_bin_limits[i+1]),10)

    area=(np.amax(ra)-np.amin(ra))*(np.amax(dec)-np.amin(dec))
    fractional_area=area/(3.046*10**-4)
    adjusted_bin_count=fractional_area*bin_count
    
    mean_counts=np.random.poisson(adjusted_bin_count)
    source_map=np.zeros_like(z)
    
    for i in range(len(mean_counts)):
        if mean_counts[i]!=0:
            for j in range(int(mean_counts[i])):
                source_ra=random.choice(ra)
                source_dec=random.choice(dec)
                source_flux=random.choice(flux_bins[i,:])
                s_x=random.randint(0,step_no-1)
                s_y=random.randint(0,step_no-1)
                source_map[s_x,s_y]=source_map[s_x,s_y]+source_flux
                
    z=z+source_map
    
    obs_time=3600 #observing time in seconds
    
    sigma=T_sys/((freq_upper-freq_lower)*obs_time)**0.5 #rms noise
    sigma_Jy=3514*sigma/(aperture_eff*dishsize**2)    #sigma in Jy

    ra_units=[]
    ra_array=np.zeros(len(ra))
    dec_units=[]
    dec_array=np.zeros(len(dec))
    for i in range(len(ra)):
        ra_units.append(RightAscension(ra[i]))
        dec_units.append(Declination(dec[i]))
        ra_array[i]=ra_units[i].seconds
        dec_array[i]=dec_units[i].seconds

    #adding noise

    noise=abs(np.random.normal(loc=sigma_Jy,scale=sigma_Jy,size=np.shape(z)))
    z=z+noise


    band7_freq=3.0*10**11
    beam_fwhm=1.02*((3*10**8/band7_freq)/dishsize)
    
    if convolve is 'true':
        #convolving with primary beam
        B_Gaussian=functional_models.Gaussian2D(x_mean=disk_ra,y_mean=disk_dec,x_stddev=beam_fwhm/2.355,y_stddev=beam_fwhm/2.355)

        Beam=B_Gaussian(xx,yy)

        image=signal.convolve2d(z,Beam,mode='same',boundary='wrap')
        image=np.log10(image)
    else:
        image=np.log10(z)

    #levels=np.linspace(np.amin(z),np.amax(z),40)
    levels=np.linspace(np.amin(image),np.amax(image),40)
    #levels=np.append(np.array([0]),np.logspace(-4,np.log10(np.amax(image)),9))
    cmap=clr.ListedColormap(['#400040','#4a1454','#542868','#5e3c7c','#685090','#7264a4','#7c78b8','#868cc2','#90a0cc','#9ab4d6'])

    
    fig=plt.figure(figsize=(10,7.5))
    ax1=fig.add_subplot(111)
    im1=plt.contourf(ra,dec,image,levels)
    ax2=ax1.twinx()
    ax2.plot(dec,ra_array,alpha=0.0)
    ax3=ax1.twiny()
    ax3.plot(dec_array,ra,alpha=0.0)

    ax1.set_xlim(np.amin(ra),np.amax(ra))
    ax1.set_ylim(np.amin(dec),np.amax(dec))
    ax1.set_xlabel('Right ascension / seconds',labelpad=25)
    ax1.set_ylabel('Declination / seconds',labelpad=40)
    ax1.tick_params('both',labelbottom=0,labeltop=0,labelleft=0,labelright=0)
    ax2.tick_params('y',labelleft='on',labelright='off')
    ax3.tick_params('x',labeltop='off',labelbottom='on')

    fig.text(0.07,0.93,str(dec_units[0].degrees)+' deg '+str(dec_units[0].minutes)+' min')
    fig.text(0.7,0.05,str(ra_units[0].hours)+' hrs '+str(ra_units[0].minutes)+' min')
    fig.text(0.9,0.5,'log10(flux/Jy beam-1)',rotation=90)
    ax2.xaxis.set_major_formatter(FormatStrFormatter('%.2f'))
    ax3.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))

    fig.colorbar(im1,ax=[ax1,ax2,ax3],pad=0.1)

    if convolve is 'true':
        beam_circ=plt.Circle((np.amin(ra)+beam_fwhm/2+2*resolution,np.amin(dec)+beam_fwhm/2+2*resolution),radius=beam_fwhm/2,color='white',fill=False)
        ax1.add_patch(beam_circ)
    
    plt.title('HD 100546 in Band 7 at '+str(distance)+' pc')

    #fig.savefig('convolved_disk_'+str(int(distance))+'pc.png')
    plt.show()
    
    return [z,ra,dec]

def fits_image(distance):
    
    data=disk_contour(distance,convolve='false')
    skymap=data[0]

    fits_data=skymap[:,:,np.newaxis,np.newaxis]
    
    ra=data[1]*(180.0/np.pi)
    dec=data[2]*(180.0/np.pi)
    hdu=fits.PrimaryHDU(fits_data)

    prihdr=fits.Header()

    prihdr['BSCALE'] = 1.0
    prihdr['BZERO'] = 0.0
    prihdr['BTYPE'] = 'Intensity'
    prihdr['OBJECT'] = ''
    prihdr['BUNIT'] = 'Jy/pixel'
    prihdr['LONPOLE'] = ra[int(np.ceil(len(ra)/2))]
    prihdr['LATPOLE'] = dec[int(np.ceil(len(dec)/2))]

    prihdr['RESTFRQ'] = 3.45*10**11
    
    prihdr['CTYPE1'] = 'RA---SIN'
    prihdr['CRVAL1'] = ra[int(np.ceil(len(ra)/2))]
    prihdr['CDELT1'] = ra[1]-ra[0]
    prihdr['CROTA1'] = 0.0
    prihdr['CRPIX1'] = np.ceil(len(ra)/2)
    prihdr['CUNIT1'] = 'deg '

    prihdr['CTYPE2'] = 'DEC---SIN'
    prihdr['CRVAL2'] = dec[int(np.ceil(len(dec)/2))]
    prihdr['CDELT2'] = dec[1]-dec[0]
    prihdr['CROTA2'] = 0.0
    prihdr['CRPIX2'] = np.ceil(len(dec)/2)
    prihdr['CUNIT2'] = 'deg '

    prihdr['CTYPE3'] = 'FREQ '
    prihdr['CRVAL3'] = 3.45*10**11
    prihdr['CDELT3'] = 7.5*10**9
    prihdr['CROTA3'] = 0.0
    prihdr['CRPIX3'] = 1.0
    prihdr['CUNIT3'] = 'Hz '

    prihdr['CTYPE4'] = 'STOKES '
    prihdr['CRVAL4'] = 1.0
    prihdr['CDELT4'] = 1.0
    prihdr['CROTA4'] = 0.0
    prihdr['CRPIX4'] = 1.0
    prihdr['CUNIT4'] = ' '


    hdu=fits.PrimaryHDU(fits_data,header=prihdr)

    hdulist=fits.HDUList([hdu])

    hdulist.writeto('sky_map_'+str(int(distance))+'pc.fits')

def change_fits(distance):

    data=disk_contour(distance,convolve='false')
    ra=data[1]*(180.0/np.pi)
    dec=data[2]*(180.0/np.pi)

    skymap=data[0]

    fits_data=skymap[np.newaxis,np.newaxis,:,:]
    new_hdu=fits.PrimaryHDU(fits_data)

    hdulist = fits.open('modelsky.fits')
    in_header=hdulist[0].header


    in_header['LONPOLE'] = ra[int(np.ceil(len(ra)/2))]
    in_header['LATPOLE'] = dec[int(np.ceil(len(dec)/2))]

    in_header['RESTFRQ'] = 3.45*10**11

    in_header['CRVAL1'] = ra[int(np.ceil(len(ra)/2))]
    in_header['CDELT1'] = ra[1]-ra[0]

    in_header['CRVAL2'] = dec[int(np.ceil(len(dec)/2))]
    in_header['CDELT2'] = dec[1]-dec[0]

    hdulist[0].data=fits_data

    naxisvals=fits_data.shape

    print naxisvals

    in_header['NAXIS1']=naxisvals[0]
    in_header['NAXIS2']=naxisvals[1]
    print in_header['NAXIS1']

    hdulist.close()
    hdulist.writeto('sky_map_'+str(int(distance))+'pc.fits')
