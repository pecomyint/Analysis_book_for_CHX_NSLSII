# -*- coding: utf-8 -*-
"""
Last modified date 10/2019

@author: Peco


Can be used to do:
extract intensity vs time
extract amplification factors
linear theory analysis: first extract R and then fit them with long wave form and with orchart term
visualize the R vs q plots and fitting them
chi square optimization for the linear theory
"""

from scipy import optimize #for curve fitting
from tqdm import tqdm_notebook as tqdm #for progress bar

'''
The following code is use latex computer modern font for all plots and enable latex..
doesn't work on the CHX cluster but should work on your local computer
'''
# from matplotlib import rc
# import matplotlib.pylab as plt
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.family'] = "serif"
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'figure.figsize': (7.2,4.45)})
# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=True)





mpl.rcParams.update({'xtick.labelsize': 15})

mpl.rcParams.update({'ytick.labelsize': 15})

mpl.rcParams.update({'font.size': 14})

mpl.rcParams.update({'figure.autolayout': True})

mpl.rcParams.update({'axes.titlesize': 15})

mpl.rcParams.update({'axes.labelsize': 15})

mpl.rcParams.update({'lines.linewidth': 2})

mpl.rcParams.update({'lines.markersize': 6})

mpl.rcParams.update({'legend.fontsize': 13})

mpl.rcParams.update({'legend.fontsize': 13})


#define your sample's materials critical angle
alpha_c = 0.185*(np.pi/180) #silicon .

#ccd positions during data acquisition
ccdx=0
ccdz=0

#incident beam info
x0 = center[::-1][0] #the beam centers are reculated. the python notebook has to be run before you assign this value  #md['beam_center_x'] ##
y0 = center[::-1][1]# md['beam_center_y'] ##
ccdx0=0
ccdz0=0

#the following are detector distance, Energy and incident angles.
#it depends on how the metadata is written in your HDF5 file. disabled for now
# d=md['det_distance'] ##
# E=md['photon_energy']*1.602*1E-19 ##

d=md['eiger4m_single_det_distance']*1000 ##
Energy=md['eiger4m_single_photon_energy']*1.602*1E-19 ##

alpha_i=np.float(md['incident angle'])*(np.pi/180) ##

#time used to take a single frame
datatakingtime = np.float(md['acquire period'])  #md['cam_acquire_period'] ##


############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx  q calculation start   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################

def plot_a_detector_image(n, color_min=0.01, color_max=0, crop = (None,None,None,None), show_q = True, color_map=plt.cm.jet, save = 0, show_subtitle = True): #cmap=cmap_albula
    '''
    in crop .. (vetical_min,vertical_max,horizontal_min,horizontal_max)
    Plot a detector image ... summing available. 
    '''
    #put n as an array
    image = FD.rdframe(n[0])
    original_dimension = np.shape(image)[0]
    for i in tqdm(n[1:]):
        image += FD.rdframe(i)
    max_intensity_scattering = np.max(image[0:2100,0:900])
    image  = image[crop[0]:crop[1],crop[2]:crop[3]]
#     image[image == 0.] = np.nan #convert 0 to nan beause we will be doing log scale 
#     image = image + np.ones(np.shape(image)) #* 1E-3
    if color_max == 0:
        color_max = max_intensity_scattering #//2
    
    
    #determine the dimensions of the cropped image
    dimension = np.shape(image)
    if crop[0] == None:
        start_q_z = 0
    else:
        start_q_z = crop[0]
    if crop[1] == None:
        end_q_z = dimension[0]
    else:
        end_q_z = crop[1]
    if crop[2] == None:
        start_q_parallel = 0
    else:
        start_q_parallel = crop[2]
    if crop[3] == None:
        end_q_parallel = dimension[1]
    else:
        end_q_parallel = crop[3]
    
    q_parallel_min = np.round(convert_pixel_to_q(start_q_parallel,dimension[1])[1],decimals=3)
    q_parallel_max = np.round(convert_pixel_to_q(end_q_parallel,dimension[1])[1],decimals=3)
    q_z_min = np.round(convert_pixel_to_q(dimension[0],start_q_z)[2],decimals=3)
    q_z_max = np.round(convert_pixel_to_q(dimension[0],end_q_z)[2],decimals=3)
# extent=[horizontal_min,horizontal_max,vertical_min,vertical_max].
    
    
    #determine the a proper size for the image
    x, y = 8*(dimension[0]/original_dimension),8*(dimension[1]/original_dimension)
    if x < 8:
        multiply_factor = 8/x
    else:
        multiply_factor = 1
    x, y = x*multiply_factor, y*multiply_factor
#     if save == 1:
#         plt.figure(figsize=np.round((3.5,3.5*(y/(x+y))),decimals=0))
#     else:
    x,y = 4,4
    plt.figure(figsize=np.round((x,y),decimals=0))
    ax = plt.gca()
    if show_q == True:
        im = ax.imshow(image, origin = 'lower',\
                       cmap=color_map, norm=LogNorm(vmin=color_min, vmax=color_max),\
                       extent = [q_parallel_min,q_parallel_max,q_z_min,q_z_max],\
                      interpolation='none') #insert after image, --- norm=LogNorm(vmin=color_min, vmax=color_max),       vmin = color_min , vmax = color_max,
        ax.set_ylabel(r'$q_z$ [$nm^{-1}$]')
        ax.set_xlabel(r'$q_{//}$ [$nm^{-1}$]')
    else:
        im = ax.imshow(image, origin = 'lower',\
                   cmap=color_map, norm=LogNorm(vmin=color_min, vmax=color_max),\
                  interpolation='none') #insert after image, --- norm=LogNorm(vmin=color_min, vmax=color_max),       vmin = color_min , vmax = color_max,
        ax.set_ylabel(r'Pixels')
        ax.set_xlabel(r'Pixels')
    
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    
    if show_subtitle == True:
        None
#         ax.set_title('Detector images summed over frames No. ' + np.str(n[0]) + ' to No. ' + np.str(n[-1]) + ' ')
    
    if save == 1:
        # current date and time
        now = datetime.now()
        # timestamp = datetime.timestamp(now)
        now = np.str_(now)
        now = now.replace('-','')
        now = now.replace(':','')
        now = now.replace(' ','')
        index = now.index('.')
        now = now[:index]

        try:
            uid
        except NameError:
            print("well, uid WASN'T defined! saving in primary folder")
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name +'Detector_image_' + np.str(n[0])  + '_to_' + np.str(n[-1]) + '_'+ now + '.pdf'
        else:
            print("uid was defined. Saving in respective folder")
            directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
            if directory_exists == False:
                os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' +'Detector_image_' + np.str(n[0])  + '_to_' + np.str(n[-1]) + '_'+ now + '.pdf'
        
        plt.savefig(fp,bbox_inches='tight', format='eps', dpi=50) # 
    plt.show()

def convert_pixel_to_q(x,y):#,ccdx,ccdz,x0,y0,ccdx0,ccdz0,d,Energy,alpha_i):
    """Notice that x comes first here. sorry for inconsistency. normally Y is taken first because of 
    the way data is sliced. but for computing I would like to keep things in x and y order"""
    Psi = ((ccdx0 - ccdx)*(1E-3) + (x - x0)*(75E-6)) / (d*1E-3)
    alpha_f = ((ccdz - ccdz0)*(1E-3) + (y - y0)*(75E-6)) / (d*1E-3)
    
    h = 6.62607 * 1E-34
    c = 3E8
    lambda_ = (h*c)/Energy
    q_x = ((2*np.pi)/lambda_)*(np.cos(alpha_f)*np.cos(Psi)-np.cos(alpha_i))
    q_y = ((2*np.pi)/lambda_)*(np.cos(alpha_f)*np.sin(Psi))
    q_z = ((2*np.pi)/lambda_)*(np.sin(alpha_f)+np.sin(alpha_i))#+ omitted  selecting qz = 0 at incident beam also sin i + sin f = sin if in small angle
    
    
    #fix the alpha_f since now you have to count pixels from the plane of sample. not from incident beam
    alpha_f = (((ccdz - ccdz0)*(1E-3) + (y - y0)*(75E-6)) / (d*1E-3)) - alpha_i
    
    #to calculate qprime
    if alpha_f < alpha_c:
        alpha_f_prime = 0
    else:
        alpha_f_prime = np.sqrt(alpha_f**2-alpha_c**2)
    alpha_i_prime = np.sqrt(alpha_i**2-alpha_c**2)
    
    #calcuate qprime
    qx_prime = ((2*np.pi)/lambda_)*(np.cos(alpha_f_prime)*np.cos(Psi)-np.cos(alpha_i))
    qy_prime = ((2*np.pi)/lambda_)*(np.cos(alpha_f_prime)*np.sin(Psi))
    qz_prime = ((2*np.pi)/lambda_)*(np.sin(alpha_f_prime) + np.sin(alpha_i_prime)) #np.sin(alpha_i)+ omitted  selecting qz = 0 at incident beam
    
#     print(alpha_f_prime,alpha_f )
    
    #convert m^-1 to nm^-1
    q_x,q_y,q_z,qx_prime,qy_prime,qz_prime = q_x*1E-9,q_y*1E-9,q_z*1E-9, qx_prime*1E-9,qy_prime*1E-9,qz_prime*1E-9
    #calculate q magnitude
    q = (q_x**2+q_y**2+q_z**2)**(0.5)
    
    return(q_x,q_y,q_z,q,alpha_f,qx_prime,qy_prime,qz_prime) 

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx  q calculation end   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################


############################################################################################################
############################################################################################################
############################Simple Intensity Retriaval codes below #########################################
############################################################################################################

def get_intensity_for_pixel_box(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe):
    """
    simple code to get I vs t for a known pixel ROI
    """
    Is=np.zeros((endframe-startframe,2))
    pos = 0
    for i in tqdm(range(startframe,endframe)):
        img = FD.rdframe(i)#imgsa[i]#DD.FD.rdframe(i)
        Is[pos,0] = pos*datatakingtime #labelling time since it will be useful for fitting with log later .. use either i or pos here
        Is[pos,1] = (np.sum(img[yslicestart:ysliceend,xslicestart:xsliceend]))/((ysliceend-yslicestart)*(xsliceend-xslicestart))
        pos += 1
    return(Is)

def plot_intensity_for_pixel_box(Is):
    """
    simple code to plot I vs t for a known pixel ROI
    """
    plt.figure()
    plt.scatter(Is[:,0], Is[:,1], color = 'black', marker = '.', label='Data')
    plt.ylabel('Intensity [A.U.]')
    plt.xlabel('Time [s]')
    plt.legend(loc='best')
    plt.show()

def plot_intensity_adv_for_pixel_box(Is,yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe,ylimmin,ylimmax):
    """
    same as plot_intensity_for_pixel_box but can label things and save the output plot 
    plots the intensity that is given, plots the fitted function, save the image to the folder
    """
    plt.figure()
    plt.scatter(Is[:,0], Is[:,1], color = 'black', marker = '.', label='Data')
    plt.ylabel('Intensity [A.U.]')
    plt.xlabel('Time [s]')
    plt.legend(loc='best')
    plt.ylim(ylimmin, ylimmax)
    plt.suptitle(' \n yslicestart = {0}'.format(yslicestart) + 
                 ' \n ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' \n xsliceend = {0}'.format(xsliceend) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) , fontsize=8)
    
    plt.show()
    plt.savefig(' yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) + 
                 ' startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe))

def exponential_function(x, a, b, c):
    """
    used by functions for fitting R, b is R here. 
    """
    return a * np.exp(2*b*x) + c

def Residual_expfun(x, t, y):
    """
    used for soft L1 fits.. this computs residuals 
    """
    return x[2] + x[0] * np.exp(2* x[1] * t) - y

def plot_intensity_andfit_for_pixel_box(Is,yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe,fitsign):
    """
    plots the intensity that is given, plots the fitted function, save the image to the folder
    
    make sure to change the p0 guess when appropriate
    """
    try:
        params, params_covariance = optimize.curve_fit(exponential_function, Is[:,0], Is[:,1], p0=[1.3, fitsign*0.04,0.1]) 
    except RuntimeError:
        params, params_covariance = optimize.curve_fit(exponential_function, Is[:,0], Is[:,1], p0=[1.3, fitsign*(-1)*0.04,0.1]) 
    plt.figure()
    plt.scatter(Is[:,0], Is[:,1], color = 'black', marker = '.', label='Data')
    plt.plot(Is[:,0], exponential_function(Is[:,0], params[0], params[1],params[2]),label='Fitted function')
    plt.ylabel('Intensity [A.U.]')
    plt.xlabel('Time [s]')
    plt.legend(loc='best')
    plt.suptitle('parameters for a*exp(2*b*x)+c are {0}'.format(params) + ' \n yslicestart = {0}'.format(yslicestart) + 
                 ' \n ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' \n xsliceend = {0}'.format(xsliceend) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) , fontsize=8)
    
    plt.show()
    plt.savefig(' yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) + 
                 ' startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe))

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxx I Vs t retriving codes  START    xxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################
def read_and_plot_intensity_pixel(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe):
    Is = get_intensity_for_pixel_box(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe)
    plot_intensity_for_pixel_box(Is)
    
def read_and_plot_fit_intensity_pixel(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe):
    Is = get_intensity_for_pixel_box(yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe)
    plot_intensity_andfit_for_pixel_box(Is,yslicestart,ysliceend,xslicestart,xsliceend,startframe,endframe,fitsign)
    
def read_and_plot_fit_intensity_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe,fitsign):
    loop = np.arange(xslicestart,xsliceend,delx)
    Rs = np.zeros((len(loop),2))
    pos = 0
    for i in tqdm(loop):
        Is = get_intensity_for_pixel_box(yslicestart,ysliceend,i,i+delx,startframe,endframe)
        plot_intensity_andfit_for_pixel_box(Is,yslicestart,ysliceend,i,i+delx,startframe,endframe,fitsign)

def read_and_plot_intensity_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe,ylimmin,ylimmax):
    loop = np.arange(xslicestart,xsliceend,delx)
    Rs = np.zeros((len(loop),2))
    pos = 0
    for i in tqdm(loop):
        Is = get_intensity_for_pixel_box(yslicestart,ysliceend,i,i+delx,startframe,endframe)
        plot_intensity_adv_for_pixel_box(Is,yslicestart,ysliceend,i,i+delx,startframe,endframe,ylimmin,ylimmax)
    
############################################################################################################
########   xxxxxxxxxxxxxxxxxx I Vs t retriving codes  END    xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################


############################################################################################################
########   xxxxxxxxxxxxxxxxxx Amplification Factor retriving codes  START    xxxxxxxxxxxxxxxxxxx ############
############################################################################################################

def read_and_get_R_along_q_parallel(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe, p0_values = [1.3,0.04,0.1], interpolate=False, Structure_factor_correction = 1, robust_fit = False):
    """
     interpolating the modular gaps may yield better results
     
     this is a code where the R values are extracted by first plotting I vs t. here I vs t is for the same time frame for all wavenumbers
     
     in this:
     You will be able to get not only R, but also other fitting parameters. 
     it will save the final results as numpy files
    
    """
    #define an empty list for params, params-std and chi. those will correspond to R. the outputs will have the same data structure as the read_and_get_R_along_q_parallel_chi_optimized
    import scipy as sp
    from scipy.optimize import least_squares
    
    
    loop = np.arange(xslicestart,xsliceend,delx)
    Is=np.zeros((len(loop),endframe-startframe,2))
    pos = 0
    for i in tqdm(range(startframe,endframe)):
        img = FD.rdframe(i)#imgsa[i]#DD.FD.rdframe(i)
        
        #interpolate the gaps
        if interpolate:
            import pandas
            df = pandas.DataFrame(img)
            mask_df = pandas.DataFrame(mask)
            mask_df = mask_df.replace(0,np.nan)
            df = df*mask_df
            df = df.interpolate(method='linear', limit_direction='forward', limit = 15, axis=1) #along x
            df = df.interpolate(method='linear', limit_direction='backward', axis=0) #along y
            df = df.replace(np.nan,0)
            img = df.values
        
        
        for j in range(len(loop)):
            Is[j,pos,0] = pos*datatakingtime #labelling time since it will be useful for fitting with log later .. use either i or pos here
            Is[j,pos,1] = Structure_factor_correction*(np.sum(img[yslicestart:ysliceend,loop[j]:loop[j]+delx]))#/((ysliceend-yslicestart)*((loop[j]+delx)-loop[j]))
        pos += 1
    Rs = np.zeros((len(loop),4))
    qs = np.zeros((len(loop),1))
    params_new = []
    params_std_new = []
    chi_new = []
    to_delete_ = []
    pos = 0
    params_init = [1.3, 0.04,0.1]
    for i in tqdm(range(len(loop))):
        qs[pos,0] = convert_pixel_to_q(loop[i]+(delx)/2,(ysliceend-yslicestart)/2)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q
        try:
            if robust_fit == True:
                res_soft_l1 = least_squares(Residual_expfun, np.array(p0_values), loss='soft_l1', f_scale=0.8, args=(Is[i,:,0], Is[i,:,1])) #for loss choose either 'soft_l1' or 'cauchy'
                params, std_deviations = [*res_soft_l1.x], [0,0] # can't compute error by fitting this way
            else:
                params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:,0], Is[i,:,1], sigma = np.sqrt(Is[i,:,1]), p0=p0_values)
                std_deviations = np.sqrt(np.diag(params_covariance))
            Rs[pos,1] = params[1]
            Rs[pos,0] = convert_pixel_to_q(loop[i]+(delx)/2,(ysliceend-yslicestart)/2)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q
            Rs[pos,2] = std_deviations[1]
            
            params_new.append(params)
            params_std_new.append(std_deviations)
            expected_frequences= exponential_function(Is[i,:,0], params[0], params[1], params[2])
            normalized_chi = (sp.stats.chisquare(Is[i,:,1], f_exp=expected_frequences)[0])/(len(Is[i,:,0])-3)
            chi_new.append(normalized_chi)
#             p0_values = params # use the previous result
            
            pos +=  1
        except RuntimeError:
            try:
                if robust_fit == True:
                    res_soft_l1 = least_squares(Residual_expfun, np.array([p0_values[0], (-1)*p0_values[1],p0_values[2]]), loss='soft_l1', f_scale=0.8, args=(Is[i,:,0], Is[i,:,1])) #for loss choose either 'soft_l1' or 'cauchy'
                    params, std_deviations = [*res_soft_l1.x], [0,0] # can't compute error by fitting this way
                else:
                    params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:,0], Is[i,:,1], sigma = np.sqrt(Is[i,:,1]), p0=[p0_values[0], (-1)*p0_values[1],p0_values[2]])
                    std_deviations = np.sqrt(np.diag(params_covariance))
                Rs[pos,1] = params[1]
                Rs[pos,0] = convert_pixel_to_q(loop[i]+(delx)/2,(ysliceend-yslicestart)/2)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q
                Rs[pos,2] = std_deviations[1]
                
                params_new.append(params)
                params_std_new.append(std_deviations)
                expected_frequences= exponential_function(Is[i,:,0], params[0], params[1], params[2])
                normalized_chi = (sp.stats.chisquare(Is[i,:,1], f_exp=expected_frequences)[0])/(len(Is[i,:,0])-3)
                chi_new.append(normalized_chi)
#                 p0_values = params # use the previous result
                
                pos +=  1
            except RuntimeError:
#                 None #let's not delete the Rs so that we keep the same dimension for Rs, Is, qs
                to_delete_.append(pos)
                pos +=  1
    Is = np.delete(Is, to_delete_, 0)# delete at last to avoid interfering with the loop
    Rs = np.delete(Rs, to_delete_, 0)
    qs = np.delete(qs, to_delete_, 0)
    #write the set parameters into your Rs
    Rs[0,3] = yslicestart
    Rs[1,3] = ysliceend
    Rs[2,3] = xslicestart
    Rs[3,3] = xsliceend
    Rs[4,3] = delx
    Rs[5,3] = startframe
    Rs[6,3] = endframe
    Rs[7,3] = interpolate*1
    
    directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
    ion_angle = md['ion angle'].replace(' ', '_')
    name = md['sample'] +'_'+ ion_angle +'_' + md['gas'] + '_'+'Energy'+  md['beam voltage']  +'_%s_'%(uid)
    if directory_exists == False:
        os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
    
    
    from datetime import datetime
    # current date and time
    now = datetime.now()
    # timestamp = datetime.timestamp(now)
    now = np.str_(now)
    now = now.replace('-','')
    now = now.replace(':','')
    now = now.replace(' ','')
    index = now.index('.')
    now = now[:index]
    
    fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + name
    np.save(fp + '_Rs_' + now, Rs)
    np.save(fp + '_Is_' + now, Is)
    np.save(fp + '_qs_' + now, qs)
    np.save(fp + '_params_' + now, np.asarray(params_new))
    np.save(fp + '_params_std_' + now, np.asarray(params_std_new))
    np.save(fp + '_chi_' + now, np.asarray(chi_new))
    return(Rs,Is,qs,np.asarray(params_new),np.asarray(params_std_new),np.asarray(chi_new))

"""a new extractions of Rs below... with Chi square optimization"""

def read_and_get_R_along_q_parallel_chi_optimized(yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe, interpolate=False, chi_end_frame = 80, increment = 20, Structure_factor_correction =1):
    """
     interpolating the modular gaps may yield better results
     
     this is a code where the R values are extracted by first plotting I vs t. here I vs t is for the DIFFERENT time frame for all wavenumbers. the difference in time is determined by first plotting the chi squre values for each fit and then optimizing the right time scale for each wave number. this may or may not work well. It may perhaps introduce some stastical noise since the chi optimization will forcefully select normalized chi squares around 1.
     
     in this:
     You will be able to get not only R, but also other fitting parameters. 
     it will save the final results as numpy files
     I vs t values will not have the same length, depending on what wavenumber you are looking at
   
    
    chi_end_frame is the starting point for doing chi sqare test.
    
    For example, startframe = 100, endframe = 1000, chi_end_frame = 100 , where datatakingtime = 1 #read from metadata
    you will see chi square vs frame plot starting from 100 (correspond to actual frame #200) until 900 correspond to actual frame #1000
    
    you will see I vs time starting from 0 (correspond to actual frame #100) until 1000 (correspond to actual frame #1000)
    """
    loop = np.arange(xslicestart,xsliceend,delx)
    Is=np.zeros((len(loop),endframe-startframe,2))
    pos = 0
    qs = np.zeros((len(loop),1))
    for i in tqdm(range(startframe,endframe)):
        img = FD.rdframe(i)
        #interpolate the gaps
        if interpolate:
            import pandas
            df = pandas.DataFrame(img)
            mask_df = pandas.DataFrame(mask)
            mask_df = mask_df.replace(0,np.nan)
            df = df*mask_df
            df = df.interpolate(method='linear', limit_direction='forward', limit = 15, axis=1) #along x
            df = df.interpolate(method='linear', limit_direction='backward', axis=0) #along y
            df = df.replace(np.nan,0)
            img = df.values
        
        
        for j in range(len(loop)):
            Is[j,pos,0] = pos*datatakingtime #labelling time since it will be useful for fitting with log later .. use either i or pos here
            Is[j,pos,1] = Structure_factor_correction*(np.sum(img[yslicestart:ysliceend,loop[j]:loop[j]+delx]))#averaging not necessary .. summing is better for chi square test
            qs[j,0] = convert_pixel_to_q(loop[j]+(delx)/2,(ysliceend-yslicestart)/2)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q
        pos += 1
    
    Is, qs_new, chi_now_new, params_new, params_std_new = chi_sqare_optimize(Is, qs, end_frame = chi_end_frame, frame_increment_for_chi = increment)
    #The following is to put the Rs in the fomat that most functions take
    Rs = np.zeros((len(qs_new),4))
    pos = 0
    for i in tqdm(range(len(qs_new))):
        Rs[pos,1] = params_new[i][1]
        Rs[pos,0] = qs_new[i]
        Rs[pos,2] = params_std_new[i][1]
        pos +=  1
    #write the set parameters into your Rs
    Rs[0,3] = yslicestart
    Rs[1,3] = ysliceend
    Rs[2,3] = xslicestart
    Rs[3,3] = xsliceend
    Rs[4,3] = delx
    Rs[5,3] = startframe
    Rs[6,3] = endframe
    Rs[7,3] = interpolate*1
    Rs[8,3] = chi_end_frame 
    Rs[9,3] = increment
 
    from datetime import datetime
    # current date and time
    now = datetime.now()
    # timestamp = datetime.timestamp(now)
    now = np.str_(now)
    now = now.replace('-','')
    now = now.replace(':','')
    now = now.replace(' ','')
    index = now.index('.')
    now = now[:index]
    
    directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
    ion_angle = md['ion angle'].replace(' ', '_')
    name = md['sample'] +'_'+ ion_angle +'_' + md['gas'] + '_'+'Energy'+  md['beam voltage']  +'_%s_'%(uid)
    if directory_exists == False:
        os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
    fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + name
    np.save(fp + '_Chi_optimized_Rs_' + now, Rs)
    np.save(fp + '_Chi_optimized_Is_' + now, Is)
    np.save(fp + '_Chi_optimized_qs_' + now, qs_new)
    np.save(fp+ '_Chi_optimized_chi_' + now , chi_now_new)
    np.save(fp+ '_Chi_optimized_params_' + now , params_new)
    np.save(fp+'_Chi_optimized_params_std_new_'  + now,params_std_new)
    
    return(Rs,Is,qs_new,params_new, params_std_new, chi_now_new) 

def chi_sqare_optimize(Is, qs, end_frame=200, frame_increment_for_chi = 20):
    import scipy as sp
    #the following loop goes through each ROI
    I_new = []
    q_new = []
    chi_new = []
    params_new = []
    params_std_new = []
    Num_ROIs = np.int(np.shape(Is)[0])
    
    param_bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf])
    for i in tqdm(range(Num_ROIs)): #loops through all the ROIs
        chi_end_frame = end_frame
        append_ready = True #we will append the new values for the above .._new = [] only if fits are successful        
        try:
            params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:chi_end_frame,0], Is[i,:chi_end_frame,1], sigma = np.sqrt(Is[i,:chi_end_frame,1]), p0=[1.3, (1)*0.04,0.1], bounds = param_bounds)#
            std_deviations = np.sqrt(np.diag(params_covariance))
        except RuntimeError:
            try:
                params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:chi_end_frame,0], Is[i,:chi_end_frame,1], sigma = np.sqrt(Is[i,:chi_end_frame,1]), p0=[1.3, (-1)*0.04,0.1], bounds = param_bounds) #
                std_deviations = np.sqrt(np.diag(params_covariance))
            except RuntimeError:

                append_ready = False
            except ValueError:
                append_ready = False
        except ValueError:
            append_ready = False
    

        if append_ready:
            expected_frequences= exponential_function(Is[i,:chi_end_frame,0], params[0], params[1], params[2])
            normalized_chi = (sp.stats.chisquare(Is[i,:chi_end_frame,1], f_exp=expected_frequences)[0])/(len(Is[i,:chi_end_frame,0])-3)
            # the following lists for to be only used in the while loops and will be emptied again for a new ROIs
            chi = []
            frame_for_chi = []
            params_for_chi = []
            params_std_for_chi = []
            
            chi.append(normalized_chi)
            frame_for_chi.append(chi_end_frame)
            params_for_chi.append(params)
            params_std_for_chi.append(std_deviations)
            
            append_ready_2 = True #this is to stop looking if the fit fails in the while loop .. enable in the line below.. i don't know if it is worth it... maybe after some frames you may be able to fit again
            while chi_end_frame <= np.shape(Is)[1] and chi[-1] < 2.5:# and append_ready_2:
                chi_end_frame = chi_end_frame + frame_increment_for_chi 
                try:
                    params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:chi_end_frame,0], Is[i,:chi_end_frame,1], sigma = np.sqrt(Is[i,:chi_end_frame,1]),p0=[1.3, 0.04,0.1], bounds = param_bounds) # 
                    std_deviations = np.sqrt(np.diag(params_covariance))
                except RuntimeError:
                    try:
                        params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:chi_end_frame,0], Is[i,:chi_end_frame,1], sigma = np.sqrt(Is[i,:chi_end_frame,1]),p0=[1.3, (-1)*0.04,0.1], bounds = param_bounds) # 
                        std_deviations = np.sqrt(np.diag(params_covariance))
                    except RuntimeError:

                        append_ready_2 = False
                    except ValueError:
                        append_ready_2 = False
                except ValueError:
                    append_ready_2 = False
                
                if append_ready_2:
                    expected_frequences= exponential_function(Is[i,:chi_end_frame,0], params[0], params[1], params[2])
                    normalized_chi = (sp.stats.chisquare(Is[i,:chi_end_frame,1], f_exp=expected_frequences)[0])/(len(Is[i,:chi_end_frame,0])-3)

                    chi.append(normalized_chi)
                    frame_for_chi.append(chi_end_frame)
                    params_for_chi.append(params)
                    params_std_for_chi.append(std_deviations)
        if append_ready:
            last_chi = chi[-1]
            choice_of_chi = last_chi
            if len(chi) > 5:
                from scipy.stats import mode
                search_chi = np.asarray(chi)
                search_chi = np.round(search_chi*20,decimals = 0)
                modeis= mode(search_chi.tolist())[0]
                boolianlist = search_chi == modeis
                the_index = np.flatnonzero(boolianlist * search_chi)[-1] #take the last index of the chosen chi squares
             
                #make sure the chosen chi squre doesn't have a huge error bar
                chosen_indexes = np.flatnonzero(boolianlist * search_chi)
                #now params for the particular index
                last_index = -1
                the_params = params_for_chi[np.int(the_index)]
                the_params_std = params_std_for_chi[np.int(the_index)]
                while (np.abs(the_params[1])/4 < np.abs(the_params_std[1])) and (last_index > (-1)*len(chosen_indexes)): #for high errors ... those with high errors here choosing 1/3 
                       
                    last_index = last_index-1 
                    the_index = chosen_indexes[last_index] # change the index and go to a smaller index
                    the_params = params_for_chi[np.int(the_index)]
                
                if last_index != -1:
                    print('See graph below. The last index of mode chi has high error in Standard Deviation. A new  selected frame... = ' + np.str_(frame_for_chi[the_index]) )
                       
                 
                #now after doing all of these effort, if the chi square is smaller than 0.6 and bigger than 1.4 use find nearest to 1            
                if np.float(chi[np.int(the_index)]) < 0.8 or np.float(chi[np.int(the_index)]) > 1.3:
                    chosen_chi, chosen_indexes = find_nearest(chi,[0.95,0.975,1,1.025,1.05,1.075,1.1,1.125,1.25])
                    #make sure the chosen chi squre doesn't have a huge error bar
                    #now params for the particular index
                    last_index = -1
                    the_index = chosen_indexes[-1]
                    the_params = params_for_chi[np.int(the_index)]
                    the_params_std = params_std_for_chi[np.int(the_index)]
                    while (np.abs(the_params[1])/4 < np.abs(the_params_std[1])) and (last_index > (-1)*len(chosen_indexes)): #for high errors ... those with high errors here choosing 1/3 
                        last_index = last_index-1 
                        the_index = chosen_indexes[last_index] # change the index and go to a smaller index
                        the_params = params_for_chi[np.int(the_index)]
                    print('See graph below. The Mode chi search gave chi square which is smaller than 0.6 and bigger than 1.4. Therefore, searched chi square around 1 manually. A new  selected frame... = ' + np.str_(frame_for_chi[the_index]))
#                     if last_index != -1:
#                         print('See graph below. The last index of manual chi square pick did not work. A new index selected... index = ' + np.str_(last_index) )

                  
                plt.figure()
                plt.scatter(frame_for_chi[:-1],chi[:-1],label = '{0}'.format(np.round(np.float(qs[i]), decimals = 4)) + r' $nm^{-1}$')
                frames_modified = np.asarray(frame_for_chi)
                chi_modified = np.asarray(chi)
                plt.scatter(frames_modified.ravel()[np.flatnonzero(boolianlist * search_chi)],chi_modified.ravel()[np.flatnonzero(boolianlist * search_chi)],c = 'r', label = 'Desired range')
                plt.ylabel(r'$\chi^2_\nu$')
                plt.xlabel('frame #')
                plt.legend(loc='best')
                plt.show()
                
                #after all of this effort, if b is positive in a is - in ae^(-2bt)+c fit, disregard this fit
                #this is unphysical scenario. the intensity drops even lower than the very intial value
                if (params_for_chi[np.int(the_index)][1] > 0 ) and (params_for_chi[np.int(the_index)][0] < 0):
                    print('Neglected the fit since R > 0 and a < 0')
                else:
                    I_new.append(Is[i,:frame_for_chi[np.int(the_index)],:])
                    q_new.append(qs[i])
                    chi_new.append(np.float(chi[np.int(the_index)]))
                    params_new.append(params_for_chi[np.int(the_index)])
                    params_std_new.append(params_std_for_chi[np.int(the_index)])

                plt.figure()
                plt.scatter(Is[i,:frame_for_chi[np.int(the_index)],0], Is[i,:frame_for_chi[np.int(the_index)],1], marker = '.', label='q = {0}'.format(np.round(np.float(qs[i]), decimals = 4)) + r' $nm^{-1}. \chi^2_v$' + ' ={0}'.format(np.round(np.float(chi[np.int(the_index)]), decimals = 3)) + '\nR = {0}'.format(np.round(np.float(params_for_chi[np.int(the_index)][1]), decimals = 3) *1E3) + ' ' + r'$ 10^{-3} s^{-1}$' + ' frame # = {0}'.format(frame_for_chi[np.int(the_index)]) + '\na = {0}'.format(np.round(np.float(params_for_chi[np.int(the_index)][0]), decimals = 3) *1E3) + ' c = {0}'.format(np.round(np.float(params_for_chi[np.int(the_index)][2]), decimals = 3) *1E3))
                plt.plot(Is[i,:frame_for_chi[np.int(the_index)],0], exponential_function(Is[i,:frame_for_chi[np.int(the_index)],0], params_for_chi[np.int(the_index)][0], params_for_chi[np.int(the_index)][1],params_for_chi[np.int(the_index)][2]))
                plt.ylabel('Intensity [A.U.]')
                plt.xlabel('Time [s]')
                plt.legend(loc='best')
                plt.show()
    return(np.asarray(I_new), np.asarray(q_new), np.asarray(chi_new), np.asarray(params_new),np.asarray(params_std_new))


############################################################################################################
########   xxxxxxxxxxxxxxxxxx Amplification Factor retriving codes  END    xxxxxxxxxxxxxxxxxxx ############
############################################################################################################


############################################################################################################
########   xxxxxxxxxxxxxxxxxx Amplification Factor plotting codes  START    xxxxxxxxxxxxxxxxxxx ############
############################################################################################################

def stitch(a,b):
    """
    used to stich the R values generated by read_and_get_R_along_q_parallel or the optimized one
    """
    c = np.zeros((len(a)+len(b),4))
    c[0:len(a),:] = a
    c[len(a):,:]=b
    return(c)

def plotRs(Rs, manual_name=''):#,yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe,chitest=" no chi test "):
    '''
    This plots the R values vs q. 
    '''
    #write the set parameters into your Rs
    yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe,interpolate, chi_end_frame, increment = 0,0,0,0,0,0,0,0,0,0
    yslicestart = Rs[0,3]
    ysliceend = Rs[1,3]
    xslicestart = Rs[2,3] 
    xsliceend = Rs[3,3] 
    delx = Rs[4,3] 
    startframe = Rs[5,3] 
    endframe = Rs[6,3] 
    interpolate = Rs[7,3] 
    chi_end_frame = Rs[8,3]  # won't be included in all Rs
    increment = Rs[9,3]  # won't be included in all Rs
 
    
    R_copy = Rs
    to_delete = []
    for i in range(np.int(len(R_copy[:,1]))):
        if (np.abs(R_copy[:,1][i]))/4 < np.abs(R_copy[:,2][i]): #delete those with high errors here choosing 1/3 
            to_delete.append(i)

    R_copy = np.delete(R_copy,to_delete,0)
    
    
    plt.figure()
    plt.errorbar(R_copy[:,0], R_copy[:,1]*1E3,yerr=R_copy[:,2]*1E3, color='black', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label= 'R at ' + r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+np.absolute(ysliceend-yslicestart)//2)[2],decimals=3)) + r' $nm^{-1}$' )
#     plt.ylim(-50, 40)
    plt.ylabel(r'R($q_{y}$) $10^{-3}$ [$s^{-1}$]')
    plt.xlabel(r'$q_{//}$ [$nm^{-1}$]')
    plt.legend(loc='best')
#     if chi_end_frame != 0:
#         more_label = ', chi_end_frame = {0}'.format(chi_end_frame) + ', increment = {0}'.format(increment)
#     plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
#                  ', ysliceend = {0}'.format(ysliceend) + 
#                  ' \n xslicestart = {0}'.format(xslicestart) + 
#                  ', xsliceend = {0}'.format(xsliceend) + 
#                  ' \n delx = {0}'.format(delx) + 
#                  ' \n startframe = {0}'.format(startframe) + 
#                  ', endframe = {0}'.format(endframe) +more_label, fontsize=8)
    
    from datetime import datetime
    # current date and time
    now = datetime.now()
    # timestamp = datetime.timestamp(now)
    now = np.str_(now)
    now = now.replace('-','')
    now = now.replace(':','')
    now = now.replace(' ','')
    index = now.index('.')
    now = now[:index]
    
    try:
        uid
    except NameError:
        print("well, uid WASN'T defined! saving in primary folder")
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name +'_Rs_extracted_at_qz_'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3)) + '_' + now + '.pdf'
    else:
        print("uid was defined. Saving in respective folder")
        directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
        if directory_exists == False:
            os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
        ion_angle = md['ion angle'].replace(' ', '_')
        name = md['sample'] +'_'+ ion_angle +'_' + md['gas'] + '_'+'Energy'+  md['beam voltage']  +'_%s_'%(uid)
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + name +'_Rs_extracted_at_qz_'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3)) + '_' + now + '.pdf'
    
    plt.savefig(fp,bbox_inches='tight')
    plt.show()
    
    
def plotRs_manuallabel(Rs,text,ylim):
    '''
    this is a simple code just to change the suptitle of the plot
    '''
    plt.figure()
    plt.scatter(Rs[:,0], Rs[:,1]*1E3, color = 'black', marker = '.', label='R')
    plt.ylabel(r'R($q_{y}$) ($10^{-3} s^{-1}$)')
    plt.xlabel(r'q ($nm^{-1}$)')
    #plt.ylim(bottom=-10.5)
    plt.legend(loc='best')
    plt.suptitle(text, fontsize=8)
    
    plt.show()
    print('the figure is not saved automatically. Please save manually')

def plot_intensity_for_all_pixel_box(Rs, Is, qs,separation=1,manual_name = '', select = False, select_q = [-0.22,0.22]):
    """will overlap the plots for all the ROIs' intensity vs time plots. 
    
    use either separation or select_q = [,insert wave numbers]..this is activted when select = True
    
    manual_name is for saving the pdf purposes.
    
    you can use this as long as Rs Is qs are generated from the same R extracting codes"""
    
    #write the set parameters into your Rs
    yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe,interpolate, chi_end_frame, increment = 0,0,0,0,0,0,0,0,0,0
    yslicestart = Rs[0,3]
    ysliceend = Rs[1,3]
    xslicestart = Rs[2,3] 
    xsliceend = Rs[3,3] 
    delx = Rs[4,3] 
    startframe = Rs[5,3] 
    endframe = Rs[6,3] 
    interpolate = Rs[7,3] 
    chi_end_frame = Rs[8,3]  # won't be included in all Rs
    increment = Rs[9,3]  # won't be included in all Rs
    
    
    
    q_interest, q_index =find_nearest(qs, select_q)
    plt.figure()
    try:
    
        if select:
            for i in range(len(q_interest)):
                plt.scatter(Is[q_index[i],:,0], Is[q_index[i],:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[q_index[i]]),decimals=3)) + r' $nm^{-1}$' )
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
#             plt.show()

        else:
            Num_ROIs = np.shape(Is)[0]
            for i in range(np.int(np.round(Num_ROIs/separation))):
                plt.scatter(Is[i*separation,:,0], Is[i*separation,:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[i*separation]),decimals=3)) + r' $nm^{-1}$' )
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
#             plt.show()
    except IndexError:
        if select:
            for i in range(len(q_interest)):
                plt.scatter(Is[q_index[i]][:,0], Is[q_index[i]][:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[q_index[i]]),decimals=3)) + r' $nm^{-1}$' )
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
#             plt.show()
#             print(a)

        else:
            Num_ROIs = np.shape(Is)[0]
            for i in range(np.int(np.round(Num_ROIs/separation))):
                plt.scatter(Is[i*separation][:,0], Is[i*separation][:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[i*separation]),decimals=3)) + r' $nm^{-1}$' )
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
    
    from datetime import datetime
    # current date and time
    now = datetime.now()
    # timestamp = datetime.timestamp(now)
    now = np.str_(now)
    now = now.replace('-','')
    now = now.replace(':','')
    now = now.replace(' ','')
    index = now.index('.')
    now = now[:index]
    
    try:
        uid
    except NameError:
        print("well, uid WASN'T defined! saving in primary folder")
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name +'_Is_extracted'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3))+ '_' + now  + '.pdf'
    else:
        print("uid was defined. Saving in respective folder")
        directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
        if directory_exists == False:
            os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
        ion_angle = md['ion angle'].replace(' ', '_')
        name = md['sample'] +'_'+ ion_angle +'_' + md['gas'] + '_'+'Energy'+  md['beam voltage']  +'_%s_'%(uid)
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + name +'_Is_extracted'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3)) + '_' + now + '.pdf'
    
    
    plt.savefig(fp,bbox_inches='tight')
    plt.show()
    
            
def fit_intensity_for_all_pixel_box(Rs, Is, qs, params, manual_name = '', separation=1 , select = False, select_q = [-0.22,0.22]):
    """will overlap the plots for all the ROIs' intensity vs time plots and fit them for you. 
    
    use either separation or select_q = [,insert wave numbers]..this is activted when select = True
    
    manual_name is for saving the pdf purposes.
    
    you can use this as long as Rs Is qs params are generated from the same R extracting codes"""
    
    #write the set parameters into your Rs
    yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe,interpolate, chi_end_frame, increment = 0,0,0,0,0,0,0,0,0,0
    yslicestart = Rs[0,3]
    ysliceend = Rs[1,3]
    xslicestart = Rs[2,3] 
    xsliceend = Rs[3,3] 
    delx = Rs[4,3] 
    startframe = Rs[5,3] 
    endframe = Rs[6,3] 
    interpolate = Rs[7,3] 
    chi_end_frame = Rs[8,3]  # won't be included in all Rs
    increment = Rs[9,3]  # won't be included in all Rs
    
    q_interest, q_index =find_nearest(qs, select_q)
    plt.figure()
    try:
    
        if select:
            for i in range(len(q_interest)):
                plt.scatter(Is[q_index[i],:,0], Is[q_index[i],:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[q_index[i]]),decimals=3)) + r' $nm^{-1}$' , alpha = 0.4)
                plt.plot(Is[q_index[i],:,0], exponential_function(Is[q_index[i],:,0], params[q_index[i]][0], params[q_index[i]][1],params[q_index[i]][2]))
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
#             plt.show()

        else:
            Num_ROIs = np.shape(Is)[0]
            for i in range(np.int(np.round(Num_ROIs/separation))):
                plt.scatter(Is[i*separation,:,0], Is[i*separation,:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[i*separation]),decimals=3)) + r' $nm^{-1}$' , alpha = 0.7)
                plt.plot(Is[i*separation,:,0], exponential_function(Is[i*separation,:,0], params[i*separation][0], params[i*separation][1],params[i*separation][2]))
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
#             plt.show()
    except IndexError:
        if select:
            for i in range(len(q_interest)):
                plt.scatter(Is[q_index[i]][:,0], Is[q_index[i]][:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[q_index[i]]),decimals=3)) + r' $nm^{-1}$' , alpha = 0.7)
                plt.plot(Is[q_index[i]][:,0], exponential_function(Is[q_index[i]][:,0], params[q_index[i]][0], params[q_index[i]][1],params[q_index[i]][2]))
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
#             plt.show()

        else:
            Num_ROIs = np.shape(Is)[0]
            for i in range(np.int(np.round(Num_ROIs/separation))):
                plt.scatter(Is[i*separation][:,0], Is[i*separation][:,1], marker = '.', label='' + r'$q_{//}$' + ' = {0}'.format(np.round(np.float(qs[i*separation]),decimals=3)) + r' $nm^{-1}$' , alpha = 0.7)
                plt.plot(Is[i*separation][:,0], exponential_function(Is[i*separation][:,0], params[i*separation][0], params[i*separation][1],params[i*separation][2]))
            plt.ylabel('Intensity [A.U.]')
            plt.xlabel('Time [s]')
            plt.legend(loc='best')
    
    
    
    from datetime import datetime
    # current date and time
    now = datetime.now()
    # timestamp = datetime.timestamp(now)
    now = np.str_(now)
    now = now.replace('-','')
    now = now.replace(':','')
    now = now.replace(' ','')
    index = now.index('.')
    now = now[:index]
    
    try:
        uid
    except NameError:
        print("well, uid WASN'T defined! saving in primary folder")
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name +'_Is_extracted_fitted'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3))+ '_' + now  + '.pdf'
    else:
        print("uid was defined. Saving in respective folder")
        directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
        if directory_exists == False:
            os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
        ion_angle = md['ion angle'].replace(' ', '_')
        name = md['sample'] +'_'+ ion_angle +'_' + md['gas'] + '_'+'Energy'+  md['beam voltage']  +'_%s_'%(uid)
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + name +'_Is_extracted_fitted'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3)) + '_' + now + '.pdf'
    
    
    plt.savefig(fp)#,bbox_inches='tight')
    plt.show()

def read_and_plot_fit_intensity_along_q_parallel_on_computed_data(Is,yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe, interpolation):
    '''
    Just input Is and it will try to fit them and plot them out for you. 
    '''
    
    for i in tqdm(range(np.shape(Is)[0])):
        plt.figure()
        plt.scatter(Is[i,:,0], Is[i,:,1], color = 'black', marker = '.', label='Data')
        try:
            params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:,0], Is[i,:,1], p0=[1.3, 1*0.04,0.1])
        except RuntimeError:
            try:
                params, params_covariance = optimize.curve_fit(exponential_function, Is[i,:,0], Is[i,:,1], p0=[1.3, 1*(-1)*0.04,0.1])
            except RuntimeError:
                print('there was a fitting error. All variables set to zero')
                params[0], params[1],params[2] = 0, 0, 0
        plt.plot(Is[i,:,0], exponential_function(Is[i,:,0], params[0], params[1],params[2]),label='Fitted function')
        plt.ylabel('Intensity [A.U.]')
        plt.xlabel('Time')
        plt.legend(loc='best')
        if interpolation == False:
            text_toprint = 'no interpolation'
        else:
            text_toprint = 'gaps are interpolated'
        plt.suptitle('parameters for a*exp(2*b*x)+c are {0}'.format(params) + ' \n yslicestart = {0}'.format(yslicestart) + 
                     ', ysliceend = {0}'.format(ysliceend) + 
                     ' \n xslicestart = {0}'.format(xslicestart) + 
                     ', xsliceend = {0}'.format(xsliceend) + 
                     ' \n startframe = {0}'.format(startframe) + 
                     ', endframe = {0}'.format(endframe) +
                     text_toprint, fontsize=8)
        plt.show()
        
#         plt.savefig(' yslicestart = {0}'.format(yslicestart) + 
#                      ' ysliceend = {0}'.format(ysliceend) + 
#                      ' xslicestart = {0}'.format(xslicestart) + 
#                      ' xsliceend = {0}'.format(xsliceend) + 
#                      ' startframe = {0}'.format(startframe) + 
#                      ' endframe = {0}'.format(endframe)+
#                      text_toprint)


############################################################################################################
########   xxxxxxxxxxxxxxxxxx Amplification Factor plotting codes  START    xxxxxxxxxxxxxxxxxxx ############
############################################################################################################   

def find_nearest(array, value):
    '''
    
    Insert an array and find the nearest value to the one you specified. it will return the found value and 
    its index
    
    '''
    q_select = []
    q_select_indexes = []
    for i in range(len(value)):
        array = np.asarray(array)
        idx = (np.abs(array - value[i])).argmin()
        q_select.append(array[idx])
        q_select_indexes.append(idx)
    return (np.asarray(q_select),np.asarray(q_select_indexes))


            
############################################################################################################
########   xxxxxxxxxxxxxxxxxxxx plot Sq vs q START   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################

from matplotlib.ticker import FormatStrFormatter
from scipy import stats

def savitzky_golay_1d(y, window_size, order, deriv=0, rate=1):
    r"""Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')
from scipy import interpolate

def read_and_get_sq_q_at_different_t(select_frames = [], startframe = 0 , endframe = 0, frame_seperation = 0, bin_num = 1, xslice = [0,2160], yslice = 1450, interpolate = True, SG = False, box_size  = 59, polynomial_order = 1, logscale = True, norm_corrector = 1, ylim = [], feed_results=False, pixelinfo = np.zeros(2),sqs = [], i_s = [], q_z = 0, peakinfo = np.zeros((2,2)), loglogpeaktime = False, peakinset= False, linreg_onpeak = False, leftbottomwidthheight = [0.26, 0.7, 0.20, 0.15], save = 0):
    """
    
    I haven't implemented just plotting selected frames.
    
    this is to plot the intensity along q// at different times.
    you have the abolity to apply SG here ... set SG = True and then box_size, polynomial_order will have to be specified.
    
   log scale = True will plot the intensity in y axis in log scale
   
   y slice is along horizontal axis of the detector
   you may want to specify xslice so that you are only looking around the peak. 
   
   i haven't developed such that the metadat are not displayed 
    """
    if feed_results==False:
        sqs=[]
        i_s = []
        loop = np.arange(startframe,endframe,frame_seperation)
        if isinstance(yslice, list):
            y_avg = (yslice[1]+yslice[0])//2
        else:
            y_avg = yslice
        fig, ax = plt.subplots()
        ax.grid(color='grey', linestyle='-', which='both', axis = 'y', linewidth=1, alpha=0.2)
        for i in tqdm(loop):
            img=FD.rdframe(i)
            for j in range(bin_num):
                img+=FD.rdframe(i+j-bin_num//2)

            pixelinfo = np.linspace(np.round(convert_pixel_to_q(xslice[0],y_avg)[1],decimals=2), np.round(convert_pixel_to_q(xslice[1],y_avg)[1],decimals=2),len(img[y_avg,xslice[0]:xslice[1]])) #- 0.017
            #get the mask
            masker = np.average(mask[yslice[0]:yslice[1],xslice[0]:xslice[1]],axis=0).tolist() 
            masker = np.asarray(masker)
            masker[:] = 1
            masker[990:1030] = 0
            
# #             print(np.where(masker==0)[0])  ....then multiply 
            if SG:
                if interpolate:
                    df = pandas.DataFrame(img)
                    mask_df = pandas.DataFrame(mask)
                    mask_df = mask_df.replace(0,np.nan)
                    df = df*mask_df
                    df = df.interpolate(method='linear', limit_direction='forward', limit = 15, axis=1) #along x
                    df = df.interpolate(method='linear', limit_direction='backward', axis=0) #along y
                    df = df.replace(np.nan,0)
                    img = df.values
                
                Zf = sgolay2d( img, window_size=box_size, order=polynomial_order)
                if isinstance(yslice, list):
                    sq = np.average(Zf[yslice[0]:yslice[1],xslice[0]:xslice[1]],axis=0)
#                     f = interpolate.interp1d(pixelinfo, sq*norm_corrector*masker, kind = 'cubic')
# #                     print('here')
#                     sq = f(pixelinfo)
                    ax.scatter(pixelinfo, sq, marker = '.', label='t = ' + np.str_(np.int(i*datatakingtime)) + ' s.') 
                else:
                    sq = Zf[yslice,xslice[0]:xslice[1]]
                    ax.scatter(pixelinfo, sq*norm_corrector*masker, marker = '.', label='t = ' + np.str_(np.int(i*datatakingtime)) + ' s.') 
            else:
                if isinstance(yslice, list):
                    sq = np.average(img[yslice[0]:yslice[1],xslice[0]:xslice[1]],axis=0)
                    ax.scatter(pixelinfo, sq*norm_corrector*masker, marker = '.', label='t = ' + np.str_(np.int(i*datatakingtime)) + ' s.')
                else:
                    sq = img[yslice,xslice[0]:xslice[1]]
                    ax.scatter(pixelinfo, sq*norm_corrector*masker, marker = '.', label='t = ' + np.str_(np.int(i*datatakingtime)) + ' s.')
#             sqs.append(sq*norm_corrector*masker)
            sqs.append(sq*norm_corrector)
            i_s.append(i)
        if ylim == [5,3.3E3]:
            ax.set_title(r'$q_z′$ = ' + '0.21' + r' $nm^{-1}$')
        else:
            ax.set_title(r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[7],decimals=2)) + r' $nm^{-1}$')
        q_z = np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2))
        print('qz and y_avg are %s %s' %(q_z,y_avg))
    else: # if the data is fed in
        fig, ax = plt.subplots()
        ax.grid(color='grey', linestyle='-', which='both', axis = 'y', linewidth=1, alpha=0.2)
        for i in tqdm(range(len(sqs))):
            ax.scatter(pixelinfo, sqs[i], marker = '.', label='t = {0:.0f}'.format(i_s[i]*datatakingtime) + ' s.')
        if peakinset == True:
            left, bottom, width, height = leftbottomwidthheight
            ax2 = fig.add_axes([left, bottom, width, height])
            ax2.set_title('Cross Section')
            if loglogpeaktime:
                ax2.scatter(peakinfo[0,:],np.abs(peakinfo[1,:]), label='Peak positions')
                tt = peakinfo[0,:]
                dt = np.abs(peakinfo[1,:])
                ax2.set_xscale('log')
                ax2.set_yscale('log')
                ax2.set_ylim(dt[-1]-dt[-1]/30,dt[0]+dt[0]/30)
                ax2.set_yticks(np.round(np.linspace(dt[-1]-dt[-1]/30,dt[0]+dt[0]/30,num=4),decimals=3).tolist())
#                 ax2.get_yaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
                
                #do linear regression
                if linreg_onpeak:
                    try:
                        slope, intercept, r_value, p_value, std_err = stats.linregress(np.log(tt[1:]),np.log(dt[1:]))
                        # print(slope, intercept, r_value, p_value, std_err )
                        xvals, yvals = np.exp(np.asarray(np.log(tt))),np.exp(slope*np.asarray(np.log(tt))+intercept)
                        ax2.plot(xvals, yvals)#label='t^({0:.2f}) power law with \nintercept = {1:.2f} nm^(-1)'.format(slope,np.exp(intercept))
        #                 ax.set_xlabel('Time [s]')
        #                 ax.set_ylabel(r'q$_{//}$ [nm$^{-1}$]')
                        power_law_slope = '{0:.2f}'.format(slope)
                        power_law_slope =  '$^{'+ power_law_slope+ '}$'
                        ax2.annotate(r't' + power_law_slope, (xvals[len(xvals)//2],yvals[(len(xvals)//3)]), fontsize=16)
                    except:
                        None
                ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
#                 ax2.xaxis.set_major_formatter(FormatStrFormatter('%.d'))
                ax2.set_ylabel(r'$q_{//}$ [$nm^{-1}$]' )
                ax2.set_xlabel(r'Time [s]')
                if peakinfo[1,:][len(dt)//2] < 0:
                    ax2.set_title('Location of \nfirst harmonic peak\n on the left side.')
                else:
                    ax2.set_title('Location of \nfirst harmonic peak\n on the right side.')
        
#                 labels = [item.get_text() for item in ax2.get_xticklabels()]
#                 labels[1] = 'Testing'
#                 print(labels)

#                 ax2.set_xticklabels([])
#                 ax2.get_xaxis().set_ticks([])
#                 ax2.set_xticks([1000])
#                 ax2.xaxis.set_major_formatter(FormatStrFormatter('%.0E'))
#                 ax2.get_xaxis().set_major_formatter(matplotlib.ticker.ScalarFormatter())
            else:
                ax2.plot(peakinfo[0,:],peakinfo[1,:], '.k-')
                tt = peakinfo[0,:]
                dt = np.abs(peakinfo[1,:])
                ax2.yaxis.set_major_formatter(FormatStrFormatter('%.2f'))
            # plt.title(r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2)) + r' $nm^{-1}$')
                ax2.set_ylabel(r'$q_{//}$ [$nm^{-1}$]' )
                ax2.set_xlabel(r'Time [s]')
                if peakinfo[1,:][len(dt)//2] < 0:
                    ax2.set_title('Location of \nfirst harmonic peak\n on the left side.')
                else:
                    ax2.set_title('Location of \nfirst harmonic peak\n on the right side.')
                
#                 ax2.set_ylim(dt[-1]-dt[-1]/30,dt[0]+dt[0]/30)
#             ax2.set_ylim()
#             ax2.set_xlim()
#             ax2.set_ylim(rms-4*std-rms,rms+3*std-rms)
            
        ax.set_title(r'$q_z′$ = ' + q_z + r' $nm^{-1}$')
    
    ax.set_ylabel(r'Normalized Intensity [$nm^4$]')
    ax.set_xlabel(r'$q_{//}$ [$nm^{-1}$]')
    ax.legend(loc='lower left', ncol=2)
    
    ax.set_xlim(pixelinfo[30],pixelinfo[-25])
    
#     plt.setp(ax2.get_xticklabels(), rotation=30, horizontalalignment='right')
    if logscale:
        if ylim:
            ax.set_ylim(ylim[0], ylim[1])
        ax.set_yscale('log')
    if save == 1:
        plt.savefig('S_vs_q_timeslice_'+np.str(uid)+'.pdf')
       
    plt.show()
    metadata = select_frames + [startframe, endframe, frame_seperation, bin_num , xslice, yslice, SG , box_size, polynomial_order]
    return(pixelinfo,sqs, i_s, q_z,metadata)
    
def read_and_get_sq_peak_q_at_different_t(select_frames = [], startframe = 0 , endframe = 0, logframe_separate = False, frame_seperation = 0, bin_num = 1, xslice = [0,2160], yslice = 1450, spline = False, SG = False, box_size  = 59, polynomial_order = 1, logscale = False):
    """
    this is to plot the intensity along q// at different times.
    you have the abolity to apply SG here ... set SG = True and then box_size, polynomial_order will have to be specified.
    
   log scale = True will plot the intensity in y axis in log scale
   
   y slice is along horizontal axis of the detector
   you may want to specify xslice so that you are only looking around the peak. 
   
   if spline is enabled, spline fits will be done and 10 more data points are filled in between to find peaks with more precision
   
   eg if 
   start frame = 1, endframe = 10
   frame_sepeartion = 3
   you will have 1st, 4th, 7th, 10th frames to look from
   then, if bin_num = 3 ||| 0,1,2 will be averaged for 1st ||| 3,4,5 will be avearged for 4th, etc
   
    """
    from scipy import interpolate as sci_interpolate
    from scipy.signal import savgol_filter as sg1
    if logframe_separate:
        loop = np.round(np.geomspace(startframe,endframe,((endframe-startframe)//frame_seperation)),decimals=0).astype(int)
    else:
        loop = np.arange(startframe,endframe,frame_seperation)
    times = []
    qs = []
    max_I_sq = []
    
    lowbounds = []
    highbounds = []
    if isinstance(yslice, list):
        y_avg = (yslice[1]-yslice[0])//2
    else:
        y_avg = yslice
    
    for i in tqdm(loop):
        img=FD.rdframe(i)
        for j in range(bin_num):
            img+=FD.rdframe(i+j-bin_num//2)
        pixelinfo = np.linspace(np.round(convert_pixel_to_q(xslice[0],y_avg)[1],decimals=2), np.round(convert_pixel_to_q(xslice[1],y_avg)[1],decimals=2),len(img[y_avg,xslice[0]:xslice[1]]))
        
        if SG:
            img = sgolay2d( img, window_size=box_size, order=polynomial_order)   
        
        if spline: # we will smoothened the splien fit again with SG if SG is enabled
            img = img[yslice[0]:yslice[1],xslice[0]:xslice[1]]
            img_shape = np.shape(img)
            x = np.arange(0, img_shape[1], 1)
            y = np.arange(0, img_shape[0], 1)
            xx, yy = np.meshgrid(x, y)
            img_raw = img[yy,xx] #np.sin(xx**2+yy**2)
            spline_f = sci_interpolate.interp2d(x, y, img_raw, kind='cubic')
            #yslice needs to be a list here!
            xnew = np.arange(0, img_shape[1], 0.08)
            ynew = np.arange(0, img_shape[0], 0.08)
#                 print(centery)
            ROI = spline_f(xnew, ynew)
            pixelinfo = np.linspace(np.round(convert_pixel_to_q(xslice[0],y_avg)[1],decimals=2), np.round(convert_pixel_to_q(xslice[1],y_avg)[1],decimals=2),len(xnew))
            The_slice = np.sum(ROI,axis=0)
            if SG:
                The_slice = sg1(The_slice, window_length = box_size, polyorder = polynomial_order)
#                 print('splined and smoothened')
        else:    # if spline is not enabled ... img which will be smoothed if SG is enabled  will be used instead 
            if isinstance(yslice, list):
                The_slice = np.average(img[yslice[0]:yslice[1],xslice[0]:xslice[1]],axis=0)
            else:
                The_slice = img[yslice,xslice[0]:xslice[1]]

        
        
        index_max_sq = np.argmax(The_slice)
        max_I = The_slice[index_max_sq]
        q = pixelinfo[index_max_sq]
        #write the peak information
        times.append(i*datatakingtime)
        qs.append(q)
        max_I_sq.append(max_I)
        
        #wrtie the half max information
        the_half_max_index_high_bound = find_nearest(The_slice[index_max_sq:], [max_I//2])[1]
        the_half_max_index_low_bound = find_nearest(The_slice[:index_max_sq], [max_I//2])[1]
        lowbound = pixelinfo[the_half_max_index_low_bound[0]]
        highbound = pixelinfo[the_half_max_index_high_bound[0]+index_max_sq]
        lowbounds.append(lowbound)
        highbounds.append(highbound)
    fig, ax = plt.subplots()
    ax.scatter(times, qs, marker = '.', label=r'Position of the peak along $q_{//}$') 
    ax.fill_between(np.asarray(times), np.asarray(lowbounds), np.asarray(highbounds), alpha=0.2, edgecolor='#1B2ACC', facecolor='#089FFF',linewidth=1, linestyle='dashdot', antialiased=True) 
    plt.title(r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2)) + r' $nm^{-1}$')
    ax.set_ylabel(r'$q_{//}$ [$nm^{-1}$]') 
    ax.set_xlabel(r'Time [s]')
    ax.legend(loc='best')
#     ax.xlim(times[0],times[-1])
    
    from matplotlib.ticker import ScalarFormatter
    if logscale:
#         plt.gca().set_xlim(times[0],times[-1])
        ax.set_xscale('log')
        ax.xaxis.set_major_formatter(ScalarFormatter())
    plt.show()
    
    plt.figure()
    plt.scatter(times, max_I_sq, marker = '.', label=r'Intensity of the peak along $q_{//}$') 
    plt.title(r'$q_z$ = ' + np.str_(np.round(convert_pixel_to_q(0,y_avg)[2],decimals=2)) + r' $nm^{-1}$')
    plt.ylabel('Intensity [A.U.]')
    plt.xlabel(r'Time [s]')
    plt.legend(loc='best')
#     plt.gca().set_ylim(1, 40000.0)
#     plt.yscale('log')
    plt.show()
    return(times,qs,lowbounds,highbounds,max_I_sq)

        

    

            
############################################################################################################
########   xxxxxxxxxxxxxxxxxxxx plot Sq vs q END   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################
def sgolay2d ( z, window_size, order, derivative=None):
    """
    """
    # number of terms in the polynomial expression
    n_terms = ( order + 1 ) * ( order + 2)  / 2.0

    if  window_size % 2 == 0:
        raise ValueError('window_size must be odd')

    if window_size**2 < n_terms:
        raise ValueError('order is too high for the window size')

    half_size = window_size // 2

    # exponents of the polynomial. 
    # p(x,y) = a0 + a1*x + a2*y + a3*x^2 + a4*y^2 + a5*x*y + ... 
    # this line gives a list of two item tuple. Each tuple contains 
    # the exponents of the k-th term. First element of tuple is for x
    # second element for y.
    # Ex. exps = [(0,0), (1,0), (0,1), (2,0), (1,1), (0,2), ...]
    exps = [ (k-n, n) for k in range(order+1) for n in range(k+1) ]

    # coordinates of points
    ind = np.arange(-half_size, half_size+1, dtype=np.float64)
    dx = np.repeat( ind, window_size )
    dy = np.tile( ind, [window_size, 1]).reshape(window_size**2, )

    # build matrix of system of equation
    A = np.empty( (window_size**2, len(exps)) )
    for i, exp in enumerate( exps ):
        A[:,i] = (dx**exp[0]) * (dy**exp[1])

    # pad input array with appropriate values at the four borders
    new_shape = z.shape[0] + 2*half_size, z.shape[1] + 2*half_size
    Z = np.zeros( (new_shape) )
    # top band
    band = z[0, :]
    Z[:half_size, half_size:-half_size] =  band -  np.abs( np.flipud( z[1:half_size+1, :] ) - band )
    # bottom band
    band = z[-1, :]
    Z[-half_size:, half_size:-half_size] = band  + np.abs( np.flipud( z[-half_size-1:-1, :] )  -band )
    # left band
    band = np.tile( z[:,0].reshape(-1,1), [1,half_size])
    Z[half_size:-half_size, :half_size] = band - np.abs( np.fliplr( z[:, 1:half_size+1] ) - band )
    # right band
    band = np.tile( z[:,-1].reshape(-1,1), [1,half_size] )
    Z[half_size:-half_size, -half_size:] =  band + np.abs( np.fliplr( z[:, -half_size-1:-1] ) - band )
    # central band
    Z[half_size:-half_size, half_size:-half_size] = z

    # top left corner
    band = z[0,0]
    Z[:half_size,:half_size] = band - np.abs( np.flipud(np.fliplr(z[1:half_size+1,1:half_size+1]) ) - band )
    # bottom right corner
    band = z[-1,-1]
    Z[-half_size:,-half_size:] = band + np.abs( np.flipud(np.fliplr(z[-half_size-1:-1,-half_size-1:-1]) ) - band )

    # top right corner
    band = Z[half_size,-half_size:]
    Z[:half_size,-half_size:] = band - np.abs( np.flipud(Z[half_size+1:2*half_size+1,-half_size:]) - band )
    # bottom left corner
    band = Z[-half_size:,half_size].reshape(-1,1)
    Z[-half_size:,:half_size] = band - np.abs( np.fliplr(Z[-half_size:, half_size+1:2*half_size+1]) - band )

    # solve system and convolve
    if derivative == None:
        m = np.linalg.pinv(A)[0].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, m, mode='valid')
    elif derivative == 'col':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -c, mode='valid')
    elif derivative == 'row':
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid')
    elif derivative == 'both':
        c = np.linalg.pinv(A)[1].reshape((window_size, -1))
        r = np.linalg.pinv(A)[2].reshape((window_size, -1))
        return scipy.signal.fftconvolve(Z, -r, mode='valid'), scipy.signal.fftconvolve(Z, -c, mode='valid')