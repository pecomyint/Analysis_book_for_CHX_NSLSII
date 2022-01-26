"""
Last modified date 01/12/2020

@author: Peco

Can be used to do:
g2 analysis with or without multi tau method
extract the exponents from the generated g2 plots. 
plot the exponents

TTCF analysis. Plot Two-time plots, do cross sections
and fit the cross sections as a function of T
The ablity to do 4 different normalizations 
and the ability to bin frames together for higher intensity 
"""
from scipy import optimize
from tqdm import tqdm_notebook as tqdm #for progress bar
import scipy as sp

'''
The following code is use latex computer modern font for all plots and enable latex..
doesn't work on the CHX cluster but should work on your local computer
'''
# from matplotlib import rc
# import matplotlib.pylab as plt

# rc('font', **{'family': 'serif', 'serif': ['Computer Modern']})
# rc('text', usetex=False)
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

# mpl.rcParams.update({'legend.fontsize': 13})
mpl.rcParams.update({'legend.fontsize': 13})

mpl.rcParams.update({'legend.fontsize': 13})

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx  define functions start   xxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################
def KWW(deltaT,beta,tau,n):
    """for a given q"""
    return 1+beta*np.exp(-2*((deltaT)/tau)**n)

def KWW(deltaT,beta,tau,n,B):
    """for a given q"""
    return B+beta*np.exp(-2*((deltaT)/tau)**n)

def linth4TTCF(deltaT,beta,R,T,B):
    """M. Mokhtarzadeh and K. Ludwig, J. Synchrotron Radiat. 24, 1187 (2017)"""
    return B+beta*((np.exp(2*R*deltaT) - np.exp(R*(deltaT + T)))/(1 - np.exp(R*(deltaT + T))))

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx  define functions start   xxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################




############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx  plot g2 functions start   xxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################

def plotandfit_g2(timeinfo, g2, p0_input, startframe, endframe, centery, centerx, boxsize, bounds_input=([-1,-np.inf,-3.0,0.999],[1,np.inf,3.0,1]),save=0, plot = 0):
    #Sy = 0
    params, params_covariance = optimize.curve_fit(KWW, timeinfo, g2,bounds=bounds_input, p0=p0_input)#[0.15,50,1.2,1.0]) 
    if plot == 1:
        plt.figure()
        plt.scatter(timeinfo, g2, color = 'black', marker = '.', label= r'$g_2$ at ' + r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$, ' + r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$' ) 
        plt.plot(timeinfo, KWW(timeinfo, params[0],params[1],params[2],params[3]),label=r'$B + \beta (q) e^{-2 (\frac{\Delta t}{\tau (q)})^{n(q)}}$ fit') 
    #     observed_values= g2
    #     expected_values = KWW(timeinfo, params[0],params[1],params[2])
    #     test = sp.stats.chisquare(observed_values, f_exp=expected_values)
    #     p_value = test[1]
    #     chi_square = test[0]
        plt.ylabel(r'$g_{2}(q,\Delta t)$')
        plt.xlabel(r'$\Delta t$ [s]')
        plt.legend(loc='best')
        if save != 1:
            plt.suptitle('startframe = {0}'.format(startframe) + 
                         ', endframe = {0}'.format(endframe) + 
                         ', center Y = {0}'.format(centery) + 
                         ', center X = {0}'.format(centerx) + 
                         ', Delx = {0}'.format(boxsize) +
                         '\n ' +
                         r'$B(q)$ = {0}'.format(params[3]) +
                         r',$\beta (q)$ = {0}'.format(params[0]) 
                         + r', $\tau (q)$ = {0}'.format(params[1]) +
                         r', $n (q)$ = {0}'.format(params[2]), fontsize=8, y = 1.1)
        plt.rcParams.update({'font.size': 14})
        plt.rcParams['axes.facecolor'] = 'white'
        plt.grid(False)
        plt.box(True)
        plt.xscale('log')
        if save == 1:
            plt.savefig('g2.eps')
        plt.show()
    return(params, params_covariance)

#     plt.savefig(DD.DD['filename'][:-4]+'_start{0}'.format(startframe) + 
#                  '_end{0}'.format(endframe) + 
#                  '_Y{0}'.format(centery) +
#                  '_X{0}'.format(centerx) + 
#                  '_box{0}'.format(boxsize) +
#                  '_beta{0}'.format(params[0]) + 
#                 '_tau{0}'.format(params[1])  +
#                  '_n{0}'.format(params[2]) + '_KWW parameters_extracted')# '_q{0}'.format(convert_pixel_to_q(centerx,centery)[1]))

    
    
#     print(chi_square,p_value)
#plotandfit_g2(timeinfo[1:], g2[1:],700, 900, 360, 160,11)
# plt.figure()
# plt.plot(np.linspace(1,50,50), KWW(np.linspace(1,50,50), 2,2.5,3.5,5),label=r'$1 + \beta (q) e^{-2 (\frac{\Delta t}{\tau (q)})^{n(q)}}$ fit')

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxx  plot g2 functions start   xxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################



############################################################################################################
########   xxxxxxxxxxxxx extract g2 WITH or WITHOUT multi tau start   xxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################

def g2_with_SG(startframe, endframe, centery, centerx, boxsize):
    """
    Peco's code
    with plain SG smoothening

    """
    if boxsize//2 == 0:
        print('boxsize must be odd')
    else:
        halflength = np.int_((boxsize-1)/2)
    images = []
    g2 = []
    tau_val = []
    Itimeavg = np.zeros((boxsize,boxsize))
    Is_avg = np.zeros((boxsize,boxsize))
    for frame in tqdm(np.arange(startframe,endframe)):
        image = DD.FD.rdframe(frame)
        images.append(image)
#         images_smooth = sgolay2d(eliminate_module_beamstop(image)[1], window_size=boxsize, order=1) #interpolate and SG
        images_smooth = sgolay2d(image, window_size=boxsize, order=1) # just SG -- faster
        Is_avg += images_smooth[centery-halflength:centery+halflength+1,centerx-halflength:centerx+halflength+1]
        Itimeavg += image[centery-halflength:centery+halflength+1,centerx-halflength:centerx+halflength+1]
    N = len(images)
    Is_avgval = Is_avg/N
    Itimeavg = Itimeavg/N
    print('smoothened avg I ={0}'.format(np.average(Is_avgval)) + '. time averaged I ={0}'.format(np.average(Itimeavg)))
    for tau in tqdm(range(N)):
        Summation = np.zeros((boxsize,boxsize))
        for i in range(N-tau):
            Summation += images[i][centery-halflength:centery+halflength+1,centerx-halflength:centerx+halflength+1] * images[i+ tau][centery-halflength:centery+halflength+1,centerx-halflength:centerx+halflength+1]     
        g2s = (Summation/(N-tau))/(Is_avgval**2)
        g2.append(np.average(g2s))
        tau_val.append(tau) 
    #print(np.shape(Is_avgval),np.shape(Summation))
    return(np.asarray(tau_val)*datatakingtime, np.asarray(g2)) 

def read_and_get_exponents_along_q_parallel(startframe, endframe, centery, centerxstart, centerxend, boxsize, tolerence):
    """
    peco's code. It does not follow the multi-tau method
    make sure you run "read_and_plot_fit_intensity_pixel" first to make sure the fit is decent then run this
    
    also make sure to change p0 guess whenver appropriate so that the fitting converges
    
    tolerence below 0.01 is a good choice
    """
    loop = np.arange(centerxstart,centerxend,boxsize)
    Vals = np.zeros((len(loop),4))
    STD = np.zeros((len(loop),4))
    pos = 0
    for i in tqdm(loop):
        timeinfo, g2 = g2_with_SG(startframe, endframe, centery, i, boxsize)
        #######
        timeinfo, g2 = timeinfo[3:-30], g2[3:-30] # dunno why the last thirty points are bad.. my g2 calculation code is not right ?
#         try:
        plotandfit_g2(timeinfo, g2,startframe, endframe, centery, i, boxsize)
        params, params_covariance = optimize.curve_fit(KWW, timeinfo, g2, p0=[1.1,1.1,1.2])
        observed_values= g2
        expected_values = KWW(timeinfo, params[0],params[1],params[2])
        test = sp.stats.chisquare(observed_values, f_exp=expected_values)
        p_value = test[1]
        chi_square = test[0]
        print(chi_square,p_value)
        if chi_square < tolerence:
            Vals[pos,3] = params[2]
            Vals[pos,2] = params[1]
            Vals[pos,1] = params[0]
            Vals[pos,0] = convert_pixel_to_q(i,centery)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q

            STD[pos,3] = std_deviations[2]
            STD[pos,2] = std_deviations[1]
            STD[pos,1] = std_deviations[0]
            STD[pos,0] = convert_pixel_to_q(xcenter_val,ycenter_val)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q
            pos +=  1
        else:
            Vals = np.delete(Vals,pos,0)
            STD = np.delete(STD,pos,0)
#         except:
#             pass
#             Vals = np.delete(Vals,pos,0)
#             STD = np.delete(STD,pos,0)
        
    return(Vals, STD)

########################################################################################################################
########################################################################################################################
########################################################################################################################
########################################################################################################################

def read_and_get_exponents_along_q_parallel_Mark(startY,endY,delY,startX,endX,delX,startframe,endframe, Xlogscale = False, Xlog_num = 100, p0=[0.15,50,1.2,1.0],bounds = ([-1,-np.inf,-3.0,0.999],[1,np.inf,3.0,1]),normalize = 1, bin_num = 1, bin_skip = 5, nobuf=8,nolev=8,skip=1, interpolater_on = False, SG_order = 1, SG_box_size = 23, save = 0, plot_g2 = [], saveg2fits = True):
    """
    will call the functions in autocor.py. make sure to run this first. the g2 calculation in here is done by the code written
    Mark Sutton. It follows the multi-tau method
    
    skip => skips the frame of interest
    bin_num and bin_skip are for normalization 
    bin_num bin over specified number of frames for normalization and bin_skip determines if you want to do binning over the bin_num frames around the next frame or skip to another one. However you bin_num or bin_skip it. the normalzation will match with the t1 or t2 you are considering. This computational structer was conceived first before TTCF binning method. 
    
    skip and bin_skip are not in sync necessarily. 
    
    comparing to TTCF terms. 
    skip --> bin_frame
    bin_num --> SG_Bin_frame
    bin_skip ---> no equivalence in TTCF codes... it lets to skip frames while binning for SG to save computation time
    
    
    star, end del X Y represent the pixel numbers usually from 0 to 2160 in 4M eiger
    start frame and end frame should be within available frame range of the data
    Xlogscale will increase delX in log scale, the increment is determined by Xlog_num
    p0 are the initial parameters for the curefit to try while fitting for the exponents
    normalizewith_Sg = 0 usues the average frame of the entire start-end frame range
    normalizewith_Sg = 1 usues SG smootheing per frame with SG-order and SG-Box_size
    normalizewith_Sg = 2 does binning frames and averaging them and using only the localized average.
    normalizewith_Sg = 3 does binning frames and averaging them and using the smoothened frame of the average.
    
    plot_g2 must be a list.
    
    
    """
    g2fits = []
    numberofframes = endframe - startframe
    
    pixellist=np.where(mask.ravel()==1)[0] #list of valid pixels
    #partition in bins using x,y (i,j) indices
    yy,xx=np.unravel_index(pixellist,mask.shape)
    
    if Xlogscale:
        ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
        xlist=mkpartlist(np.round(np.geomspace(startX,endX,Xlog_num),decimals=0).astype(int)) # break points for bins in x
        print(xlist)
        
#         ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
#         xlist=mkpartlist(np.round(np.logspace(np.log10(startX), np.log10(endX), num=Xlog_num),decimals=0).astype(int)) # break points for bins in x
#         print(xlist)
        
    else:
        ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
        xlist=mkpartlist(np.arange(startX,endX,delX)) # break points for bins in x
    plist,bind,xbind,ybind,noperbin,inpixellist,binlist,binxlist,binylist=\
        partition2d(pixellist,xx,xlist,yy,ylist)
    
    
    
    #loop = np.arange(centerxstart,centerxend,boxsize)
    Vals = np.zeros((max(bind)+1,5))
    STD = np.zeros((max(bind)+1,5))
    pos = 0
    timeinfo1, g21 = autocor(FD,startframe,numberofframes,plist,bind,nobuf,nolev,skip, norm=normalize, bin_frames = bin_num , bin_frame_skip = bin_skip ,  interpolater = interpolater_on, norm_order = SG_order, pixel_boxsize = SG_box_size)
    timeinfo1 = timeinfo1*datatakingtime*skip
    p0_previous = p0
    
    #make a new xlist since the g21 length is different from xlist length
    xlist = np.linspace(startX,endX,max(bind)+1,dtype=int) #linspace is used since we know how many increments which is max(bind)+1
    if plot_g2:
        f1 = plt.figure()
        ax1 = f1.add_subplot(111)
    for i in tqdm(range(max(bind)+1)):
#         xcenter_val = np.int_(round((xlist[i*2]+xlist[i*2+1]-1)/2))+1
        xcenter_val = xlist[i]
        ycenter_val = np.int_(round((ylist[0*2]+ylist[0*2+1]-1)/2))+1
        #####
        timeinfo, g2 = timeinfo1[2:], g21[:,i][2:]
#         print(i)
        try:
            if save == 0:
                params, params_covariance = plotandfit_g2(timeinfo, g2, p0_previous,startframe, startframe+numberofframes, ycenter_val, xcenter_val, np.int_(xlist[i]),bounds_input = bounds, plot = 0)
            else:
                if save != i:
                    params, params_covariance = plotandfit_g2(timeinfo, g2, p0_previous,startframe, startframe+numberofframes, ycenter_val, xcenter_val, np.int_(xlist[i]),bounds_input = bounds, plot = 0)
                else:
                    params, params_covariance = plotandfit_g2(timeinfo, g2, p0_previous,startframe, startframe+numberofframes, ycenter_val, xcenter_val, np.int_(xlist[i]),bounds_input = bounds,save = 1, plot = 0)
            p0_previous = p0 #params
            std_deviations = np.sqrt(np.diag(params_covariance))
            observed_values= g2
            expected_values = KWW(timeinfo, params[0],params[1],params[2],params[3])
            test = sp.stats.chisquare(observed_values, f_exp=expected_values)
            p_value = test[1]
    #         print(timeinfo1.size())
            chi_square = test[0]/(len(observed_values)-4)
#             print(chi_square,p_value)
#             if chi_square >0.9 and chi_square < 1.6:


            g2fits.append([timeinfo, g2, KWW(timeinfo, params[0],params[1],params[2],params[3]),np.round(convert_pixel_to_q(xcenter_val,ycenter_val)[1],decimals=2),np.round(convert_pixel_to_q(xcenter_val,ycenter_val)[7],decimals=2),params[0],params[1],params[2],params[3]])
            #now plot the select wavenumbers
            if plot_g2:
                qparallelnow = np.round(convert_pixel_to_q(xcenter_val,ycenter_val)[1],decimals=2)
                print('qparallelnow  plog2', qparallelnow, plot_g2 )
                if qparallelnow in plot_g2:
                    plot_g2.pop(plot_g2.index(qparallelnow))
                    print('qparallelnow chosen  plog2', qparallelnow, plot_g2 )
                    ax1.scatter(timeinfo, g2, marker = '.')#+ r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$' ) 
                    ax1.plot(timeinfo, KWW(timeinfo, params[0],params[1],params[2],params[3]), label= r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(xcenter_val,ycenter_val)[1],decimals=2)) + r' $nm^{-1}$.')#,label=r'$B + \beta (q) e^{-2 (\frac{\Delta t}{\tau (q)})^{n(q)}}$ fit') 
                    
        #     observed_values= g2
        #     expected_values = KWW(timeinfo, params[0],params[1],params[2])
        #     test = sp.stats.chisquare(observed_values, f_exp=expected_values)
        #     p_value = test[1]
        #     chi_square = test[0]
            
#             if save != 1:
#                 plt.suptitle('startframe = {0}'.format(startframe) + 
#                              ', endframe = {0}'.format(endframe) + 
#                              ', center Y = {0}'.format(centery) + 
#                              ', center X = {0}'.format(centerx) + 
#                              ', Delx = {0}'.format(boxsize) +
#                              '\n ' +
#                              r'$B(q)$ = {0}'.format(params[3]) +
#                              r',$\beta (q)$ = {0}'.format(params[0]) 
#                              + r', $\tau (q)$ = {0}'.format(params[1]) +
#                              r', $n (q)$ = {0}'.format(params[2]), fontsize=8, y = 1.1)
            
            
        
        
            Vals[pos,4] = params[3]
            Vals[pos,3] = params[2]
            Vals[pos,2] = params[1]
            Vals[pos,1] = params[0]
            Vals[pos,0] = convert_pixel_to_q(xcenter_val,ycenter_val)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q

            STD[pos,4] = std_deviations[3]
            STD[pos,3] = std_deviations[2]
            STD[pos,2] = std_deviations[1]
            STD[pos,1] = std_deviations[0]
            STD[pos,0] = convert_pixel_to_q(xcenter_val,ycenter_val)[1]#(i+(delx)/2) #extract q... select the center pixel of the box to calculate q
            pos +=  1
#             else:
#                 Vals = np.delete(Vals,pos,0)
#                 STD = np.delete(STD,pos,0)
        except RuntimeError:
            Vals = np.delete(Vals,pos,0)
            STD = np.delete(STD,pos,0)
        except ValueError:
            Vals = np.delete(Vals,pos,0)
            STD = np.delete(STD,pos,0)
        
    exponents = np.concatenate((Vals,STD),axis=1)  
    metadata = np.zeros((pos,1))
#     qz = np.round(convert_pixel_to_q(centerx,centery)[2],decimals=2)
    qz = np.round(convert_pixel_to_q(xcenter_val,ycenter_val)[2],decimals=2)
    metadata[:15,0] = startframe, endframe, ycenter_val, startX+delX//2, endX+delX//2, qz, nobuf,nolev,skip, interpolater_on*1, SG_order, SG_box_size, bin_num, bin_skip, normalize 
    exponents = np.concatenate((exponents,metadata),axis=1)
    if plot_g2:
        ax1.set_ylabel(r'$g_{2}(q,\Delta t)$')
        ax1.set_xlabel(r'$\Delta t$ [s]')
        ax1.set_ylim(0.9,1.4)
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
#             plt.rcParams.update({'font.size': 14})
#             plt.rcParams['axes.facecolor'] = 'white'
#             plt.grid(False)
#             plt.box(True)
        ax1.set_xscale('log')
        plt.tight_layout()
        if saveg2fits == True and plot_g2:
            plt.savefig(uid+'g2.pdf')
        plt.show()
    return(exponents,g2fits)


##########################################################################################################
########   xxxxxxxxxxxxx extract g2 WITH or WITHOUT multi tau end   xxxxxxxxxxxxxxxxxxxxxxxx ############
#########################################################################################################



############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxx plot exponents start  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################

def plot_exponents(exponents, manual_name = '', save = 0):
    """ this is plotting exponents vs q// .. not to use for exponents vs big T """
    Vals, STD, metadata = exponents[:,0:5], exponents[:,5:10], exponents[:,10:]
    startframe, endframe, centery, centerxstart, centerxend, qz, nobuf,nolev,skip, interpolater_on, SG_order, SG_box_size, bin_num, bin_skip, normalize = metadata[:15,0] 
    
    
    Vals_copy = Vals
    to_delete = []
#     #delete the data where beta's STD is too high
#     for i in range(np.int(len(Vals_copy[:,1]))):
#         if (np.abs(Vals_copy[:,1][i]))/2 < np.abs(STD[:,1][i]):
#             to_delete.append(i)

#     #delete the data where tau's STD is too high
#     for i in range(np.int(len(Vals_copy[:,1]))):
#         if (np.abs(Vals_copy[:,2][i]))/3 < np.abs(STD[:,2][i]):
#             to_delete.append(i)


    Vals_fixed = np.delete(Vals_copy,to_delete,0)
    STD_fixed = np.delete(STD,to_delete,0)
    
#     to_delete = []
#     #delete the data where n's STD is too high
#     for i in range(np.int(len(STD[:,3]))):
#         if (np.abs(Vals[:,3][i]))/3 < np.abs(STD[:,3][i]):
#             to_delete.append(i)

#     Vals_fixed = np.delete(Vals_copy,to_delete,0)
#     STD_fixed = np.delete(STD,to_delete,0)
    
    # if you don't want any filter choose the following two lines.
#     Vals_fixed = Vals
#     STD_fixed = STD

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
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name + np.str_(np.round(convert_pixel_to_q(centerxstart, centery)[7],decimals=3))+'_'+ now 
        else:
            print("uid was defined. Saving in respective folder")
            directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
            if directory_exists == False:
                os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/'+ uid +  np.str_(np.round(convert_pixel_to_q(centerxstart, centery)[7],decimals=3))+'_'+  now 


    plt.figure()
    plt.errorbar(Vals_fixed[:,0], Vals_fixed[:,1],yerr=STD_fixed[:,1], color='black', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label=r'$\beta (q)$ at ' + r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerxstart, centery)[7],decimals=3)) + r' $nm^{-1}$')
    plt.ylabel(r'$\beta (q)$')
    plt.xlabel(r'$q_{//}$ [$nm^{-1}$]')
    plt.legend(loc='best')
    plt.ylim(0.05,1.25)
    plt.suptitle('center Y = {0}'.format(centery) +
                 ', centerxstart = {0}'.format(centerxstart) + 
                 ', centerxend = {0}'.format(centerxend) +
                 ', Normalization Method = {0}'.format(normalize) +
                 ' \n startframe = {0}'.format(startframe) + 
                 ', endframe = {0}'.format(endframe)+
                 ', skip = {0}'.format(skip) +
                 ', interpolater_on = {0}'.format(interpolater_on) +
                 '\n SG_order = {0}'.format(SG_order) +
                 ', SG_box_size = {0}'.format(SG_box_size)+
                 ', SG_bin_num = {0}'.format(bin_num)+
                 ', SG_bin_skip = {0}'.format(bin_skip), fontsize=8, y = 1.1)
#     plt.rcParams.update({'font.size': 12})   
    if save == 1:
        plt.savefig(fp+ '_beta.pdf',bbox_inches='tight')
    plt.show()
    #saving figures doesn't work with some versions of maplotlib .. use matplotlib==2.2.3
#     plt.savefig(DD.DD['filename'][:-4] +'_cenY{0}'.format(centery) +
#                  '_xstrt{0}'.format(centerxstart) + 
#                  '_xemd{0}'.format(centerxend) + 
#                  '_box size{0}'.format(boxsize) + 
#                  '_startframe{0}'.format(startframe) + 
#                  '_endframe{0}'.format(endframe) + '_beta_extracted')
    
    plt.figure()
    plt.errorbar(Vals_fixed[:,0], Vals_fixed[:,2],yerr=STD_fixed[:,2], color = 'red', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label=r'$\tau (q)$ at ' + r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerxstart, centery)[7],decimals=3)) + r' $nm^{-1}$')
    plt.ylabel(r'$\tau (q)$ [s]')
    plt.xlabel(r'$q_{//}$ [$nm^{-1}$]')
    plt.legend(loc='best')
    plt.ylim(0,700)
#     plt.yscale('log')
    print(Vals_fixed[:,2])
    plt.suptitle('center Y = {0}'.format(centery) +
                 ', centerxstart = {0}'.format(centerxstart) + 
                 ', centerxend = {0}'.format(centerxend) +
                 ', Normalization Method = {0}'.format(normalize) +
                 ' \n startframe = {0}'.format(startframe) + 
                 ', endframe = {0}'.format(endframe)+
                 ', skip = {0}'.format(skip) +
                 ', interpolater_on = {0}'.format(interpolater_on) +
                 '\n SG_order = {0}'.format(SG_order) +
                 ', SG_box_size = {0}'.format(SG_box_size)+
                 ', SG_bin_num = {0}'.format(bin_num)+
                 ', SG_bin_skip = {0}'.format(bin_skip), fontsize=8, y = 1.1)
    if save == 1:
        plt.savefig(fp+ '_tau.pdf',bbox_inches='tight')
    plt.show()
#     plt.savefig(DD.DD['filename'][:-4] +'_cenY{0}'.format(centery) +
#                  '_xstrt{0}'.format(centerxstart) + 
#                  '_xemd{0}'.format(centerxend) + 
#                  '_box size{0}'.format(boxsize) + 
#                  '_startframe{0}'.format(startframe) + 
#                  '_endframe{0}'.format(endframe)+ '_tau_extracted')
    
    plt.figure()
    plt.errorbar(Vals_fixed[:,0], Vals_fixed[:,3],yerr=STD_fixed[:,3], color = 'blue', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label=r'$n(q)$ at ' + r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerxstart, centery)[7],decimals=3)) + r' $nm^{-1}$')
    plt.ylabel(r'$n(q)$')
    plt.xlabel(r'$q_{//}$ [$nm^{-1}$]')
    plt.legend(loc='best')
    plt.ylim(0.2,2.1)
    plt.suptitle('center Y = {0}'.format(centery) +
                 ', centerxstart = {0}'.format(centerxstart) + 
                 ', centerxend = {0}'.format(centerxend) +
                 ', Normalization Method = {0}'.format(normalize) +
                 ' \n startframe = {0}'.format(startframe) + 
                 ', endframe = {0}'.format(endframe)+
                 ', skip = {0}'.format(skip) +
                 ', interpolater_on = {0}'.format(interpolater_on) +
                 '\n SG_order = {0}'.format(SG_order) +
                 ', SG_box_size = {0}'.format(SG_box_size)+
                 ', SG_bin_num = {0}'.format(bin_num)+
                 ', SG_bin_skip = {0}'.format(bin_skip), fontsize=8, y = 1.1)
#     plt.rcParams.update({'font.size': 12})
    if save == 1:
        plt.savefig(fp+ '_n.pdf',bbox_inches='tight')
    plt.show()
#     plt.savefig(DD.DD['filename'][:-4] +'_cenY{0}'.format(centery) +
#                  '_xstrt{0}'.format(centerxstart) + 
#                  '_xemd{0}'.format(centerxend) + 
#                  '_box size{0}'.format(boxsize) + 
#                  '_startframe{0}'.format(startframe) + 
#                  '_endframe{0}'.format(endframe)+ '_n_extracted')
    plt.figure()
    plt.errorbar(Vals_fixed[:,0], Vals_fixed[:,4],yerr=STD_fixed[:,4], color = 'blue', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label='Baseline  at ' + r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerxstart, centery)[7],decimals=3)) + r' $nm^{-1}$')
    plt.ylabel(r'Baseline')
    plt.xlabel(r'$q_{//}$ [$nm^{-1}$]')
    plt.ylim(0.9,1.1)
    plt.legend(loc='best')
    plt.suptitle('center Y = {0}'.format(centery) +
                 ', centerxstart = {0}'.format(centerxstart) + 
                 ', centerxend = {0}'.format(centerxend) +
                 ', Normalization Method = {0}'.format(normalize) +
                 ' \n startframe = {0}'.format(startframe) + 
                 ', endframe = {0}'.format(endframe)+
                 ', skip = {0}'.format(skip) +
                 ', interpolater_on = {0}'.format(interpolater_on) +
                 '\n SG_order = {0}'.format(SG_order) +
                 ', SG_box_size = {0}'.format(SG_box_size)+
                 ', SG_bin_num = {0}'.format(bin_num)+
                 ', SG_bin_skip = {0}'.format(bin_skip), fontsize=8, y = 1.1)
#     plt.rcParams.update({'font.size': 12})
    if save == 1:
        plt.savefig(fp+ '_baseline.pdf',bbox_inches='tight')
    plt.show()
#     plt.savefig(DD.DD['filename'][:-4] +'_cenY{0}'.format(centery) +
#                  '_xstrt{0}'.format(centerxstart) + 
#                  '_xemd{0}'.format(centerxend) + 
#                  '_box size{0}'.format(boxsize) + 
#                  '_startframe{0}'.format(startframe) + 
#                  '_endframe{0}'.format(endframe)+ '_n_extracted')

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxx plot exponents start  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ############
############################################################################################################








############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxx TTC calculations start  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ##########
############################################################################################################
import sys
from io import StringIO
from IPython import get_ipython


class IpyExit(SystemExit):
    """Exit Exception for IPython.

    Exception temporarily redirects stderr to buffer.
    """
    def __init__(self):
        # print("exiting")  # optionally print some message to stdout, too
        # ... or do other stuff before exit
        sys.stderr = StringIO()

    def __del__(self):
        sys.stderr.close()
        sys.stderr = sys.__stderr__  # restore from backup


def ipy_exit():
    raise IpyExit


if get_ipython():    # ...run with IPython
    exit = ipy_exit  # rebind to custom exit
else:
    exit = exit      # just make exit importable
    


import gc
import pandas
from datetime import datetime
from mpl_toolkits.axes_grid1 import make_axes_locatable # to make the color bar fixed to the graphs
from ipywidgets import interactive #for interactive plotting / making a movie for the scatterin events
# from tqdm import tqdm as tqdm #for progress bar
from scipy import signal #sg 1 d to smoothen the SG windows 



def round_up_to_odd(f):
    converted = np.ceil(f-1) // 2 * 2 + 1
    return converted.astype(int)
        
  

def collect_normalize(startframe,endframe, centerx_start, centerx_end, centery,peak_in_X, boxsize = [21,21], bin_frame = 8, normalize = 1, SG_window = 39, SG_order = 1,  SG_bin_frame = 16, interpolate = 0, dynamical_ROI_relative_positions = np.zeros(10), discrete_rel_pos = True, moving_ref_frame = False):
    '''
    warning if your ROIs are a lot or have a lot of frames, numpy won't be able to handle the huge list and say can't put (21,101) into 101 or something like that.
    boxsize = (21,21) means (x_size,y_size)

    when you bin (by setting bin_frame>1), halves of bin_frame or SG_bin_frame are used to bin around your centered frame. for example consider [1,2,3,4,5,6,7] frames and bin_frame = 3 and SG_bin_frame = 6. center frame (frame of interest) = 4 --> bin {3,4,5} and SG bin {1,2,3,4,5,6,7}

    normalize 0,1,2,3 stands for using entire time average, using individual SG applied frame, binning beteween adjacent frames (use SG_bin_frame for this), binning between adjacent frames and then applying SG on it.

    bin_frame bins frames at the beginning for higher statistic
    SG_bin_frames only bins frames are the original frame rate .. they are independent of each other 
    
    
    if you are seeing "ValueError: could not broadcast input array from shape (121,25) into shape (121)" then exclude the extreme ends in the ROIs for TTCF computation
    
    '''

    CenterXs = np.arange(centerx_start,centerx_end,boxsize[0])
    gc.collect()
    image = None
    imagesmooth = None
    img_raw = None
    img_smooth = None
    raw_f = None
    smooth_f = None


    ### stop the program if the conditions are not met
    Inor_all = []
    if boxsize[0]//2 == 0 or boxsize[1]//2 == 0:
        print('boxsize must be odd')
        exit()
    else:
        halflength_x = np.int_((boxsize[0]-1)/2)
        halflength_y = np.int_((boxsize[1]-1)/2)
        
    if bin_frame != 1:
        if SG_bin_frame < bin_frame*2:
            print('SG_bin_frames should be bigger than bin_frames. SG_bin_frame > bin_frame*2 recommended ')
            exit()

        if ((startframe - SG_bin_frame//2) <6) or ((endframe + SG_bin_frame//2) > np.int(md['eiger4m_single_cam_num_images'])):
            print('Choose a bigger startframe and/or end frame if you are binning')
            exit()

        if bin_frame > 1: #this is ...the only odd bin_frame allowed is 1
            if ((SG_bin_frame%2) != 0) or ((bin_frame%2) != 0):
                print('Choose an even number for SG_bin_frame and/or SG_bin_frame if you are binning')
                exit()
    if len(dynamical_ROI_relative_positions)!= 10:
        if len(dynamical_ROI_relative_positions) != (endframe-startframe)//bin_frame:
            print('Dynamical ROI relative positions do not have the same dimension as the number of frames specified')
            print(len(dynamical_ROI_relative_positions), (endframe-startframe)//bin_frame)
            exit()

    #conditions check done
    if normalize == 0:
        print('normalizeation = 1 selected, calcuating entire time average.')
        temp_img =FD.rdframe(startframe)
        print('calculating entire time average')
        for i in tqdm(np.arange(startframe,endframe)):
            img = FD.rdframe(startframe+i)
            # if interpolate is on
            if interpolate == 1:
                df = pandas.DataFrame(img)
                mask_df = pandas.DataFrame(mask)
                mask_df = mask_df.replace(0,np.nan)
                df = df*mask_df
                df = df.interpolate(method='linear', limit_direction='forward', limit = 15, axis=1) #along x
                df = df.interpolate(method='linear', limit_direction='backward', axis=0) #along y
                df = df.replace(np.nan,0)
                img = df.values
            temp_img += img
        norm_avg = temp_img/noframes



#     print('collect the images, normalize and extract ROIs')
    processing_frames = [] #reset processing_frames for the next scan
    contrasts = [] #to keep track the contrast as a function of T ... cross correlation between adjacent frames !
    s_contrasts = [] #to keep track the smoothened contrast as a function of T
    SG_windows = [] #see how the recommended SG_window size changes over T

    #make a mask for correlation
    pixellist=np.where(mask.ravel()==1)[0] #list of valid pixels
    #partition in bins using x,y (i,j) indices
    yy,xx=np.unravel_index(pixellist,mask.shape)

    #make the partitions 
    Xlogscale = False
    if Xlogscale:
        ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
        xlist=mkpartlist(np.round(np.geomspace(startX,endX,Xlog_num),decimals=0).astype(int)) # break points for bins in x
        print(xlist)

#         ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
#         xlist=mkpartlist(np.round(np.logspace(np.log10(startX), np.log10(endX), num=Xlog_num),decimals=0).astype(int)) # break points for bins in x
#         print(xlist)
    #CenterXs = np.arange(centerx_start,centerx_end,boxsize[0])   
    else:
        ylist=mkpartlist(np.arange(centery-boxsize[1]//2,centery+boxsize[1]//2 + 2,boxsize[1])) # break points for bins in y
        xlist=mkpartlist(np.arange(centerx_start-boxsize[0]//2,centerx_end-boxsize[0]//2 +2,boxsize[0])) # break points for bins in x
    plist,bind,xbind,ybind,noperbin,inpixellist,binlist,binxlist,binylist=\
        partition2d(pixellist,xx,xlist,yy,ylist)






    for frame in tqdm(np.arange(startframe,endframe,bin_frame)):
        gc.collect()

#         if the list is empty, read frames
        if len(processing_frames) == 0:
            list_to_look = np.arange(frame- SG_bin_frame//2,frame+ SG_bin_frame//2 +1)
#             print(list_to_look)
            for i in list_to_look:
                processing_frames.append(FD.rdframe(i))
#                 print(i)
        else:# if len(processing_frames) == SG_bin_frame+1: #first read more frames if processing_frames is not empty
            processing_frames = processing_frames[bin_frame:] # delete the frames where you won't be using 
            list_to_look = np.arange(frame + 1, frame+ 1 + bin_frame)
#             print(list_to_look)
            for i in list_to_look: #you move forward with the bin_frame rate ... that is how many more frames you will read. 
                processing_frames.append(FD.rdframe(i))
#                 print(i)

        #binning frames
        image = None #reset the image
        center_index = len(processing_frames)//2 
        select_range_to_bin = processing_frames[center_index- bin_frame//2:center_index+ bin_frame//2+1]
        image = np.sum(np.asarray(select_range_to_bin),axis=0) #now summed


        # if interpolate is on
        if interpolate == 1:
            df = pandas.DataFrame(image)
            mask_df = pandas.DataFrame(mask)
            mask_df = mask_df.replace(0,np.nan)
            df = df*mask_df
            df = df.interpolate(method='linear', limit_direction='forward', limit = 15, axis=1) #along x
            df = df.interpolate(method='linear', limit_direction='backward', axis=0) #along y
            df = df.replace(np.nan,0)
            image = df.values

        #setup mask for ccorr
        ii=0*image
        ii[np.unravel_index(plist,ii.shape)]=bind+1 #bind==0 is not a bin

        #determine a range of ROIs to do contrast averaging over
        if SG_window == 0:
            location_peak = np.int((peak_in_X/(CenterXs[-1]-CenterXs[0]))*len(CenterXs))

        #normalize the image frames
        #you can either specify SG_window size or let it be detemined automatically by cross correlation (set SG_window =0)
        imagesmooth = None #reset
        if SG_window == 0: #determine this size automatically with cross-correlation if SG_window is zero
#             if len(SG_windows) > 31:
#                 SG_window = SG_windows[-1]-12 ### initial condition 
#                 if SG_window < 11:
#                     SG_window = 11
#             else:
#                 SG_window = 11 ### initial condition 
            SG_window = 7 
            contrast = 0.30 ### arbitary start... contrast can't be bigger than 30%
            continueloop = True
            last_contrast = 0.15
            while (contrast > 0.05) and (contrast < 0.35) and (SG_window < 151) and continueloop :
                if normalize == 1:
                    imagesmooth = sgolay2d(image, window_size=SG_window, order=SG_order)        
                if normalize == 0:
                    imagesmooth = norm_avg 
                if normalize == 2: 
                    imagesmooth = np.average(np.asarray(processing_frames),axis=0)#average of SG_bin_frames (SG_bin_frame > bin_frame)
                    imagesmooth = imagesmooth*(bin_frame+1) #because the binned frame will have the summed intensity of 41 frames
                if normalize == 3:
                    imagesmooth = np.average(np.asarray(processing_frames),axis=0)#average of SG_bin_frames (SG_bin_frame > bin_frame)
                    imagesmooth = imagesmooth*(bin_frame+1) #because the binned frame will have the summed intensity of 41 frames
                    imagesmooth = sgolay2d(imagesmooth, window_size=SG_window, order=SG_order) # apply SG

                CC = crosscor(imagesmooth.shape,mask=ii,normalization="symavg")#"regular")#
                rr=CC(imagesmooth,imagesmooth)

                #average contrast over all ROIs
                contrast_temp = []
                for j in np.arange(location_peak-3,location_peak+3): #only select around the peak area
                    try:
                        a=rr[j]
                #         p=r_[5,5]+np.unravel_index(np.argmax(a[5:-5,5:-5], axis=None), (10,10))
                        p=np.unravel_index(np.argmax(a, axis=None), a.shape) #position of max
                    #             print(a.shape)
                        ax=(a[p[0]-1,p[1]]+a[p[0]+1,p[1]]-2*a[p[0],p[1]])/2.
                        dx=(a[p[0]-1,p[1]]-a[p[0]+1,p[1]])/4./ax
                        ay=(a[p[0],p[1]-1]+a[p[0],p[1]+1]-2*a[p[0],p[1]])/2.
                        dy=(a[p[0],p[1]-1]-a[p[0],p[1]+1])/4./ay
                        cy=a[p[0],p[1]]-ay*dy*dy
                        cx=a[p[0],p[1]]-ax*dx*dx
                #         res2[j,(imgno-beg)//bin_frame,:]=r_[(cx+cy)/2.,p[0]+dx,p[1]+dy,np.sqrt(-cx/2/ax),np.sqrt(-cy/2/ay)]

                        contrast_temp_append = ((cx+cy)/2.)-1
                        contrast_temp.append(contrast_temp_append)
                    except:
                        pass
                contrast = np.average(np.asarray(contrast_temp))
                SG_window += 2
                if (contrast - last_contrast) > 0.03:
                    continueloop = False
                last_contrast = contrast
            print(contrast_temp,contrast,SG_window)


            #now use SG to make the transition smooth
            if len(SG_windows) > 31:

                SG_windows_temp = signal.savgol_filter(SG_windows[-7:], window_length = 3, polyorder = 2)
                SG_windows_temp = round_up_to_odd(SG_windows_temp)
                SG_window = SG_windows_temp[-1] #take the smoothened data for smoother transition 

                #following is a copy from above. .. it's just calculating the contrast again
                if normalize == 1:
                    imagesmooth = sgolay2d(image, window_size=SG_window, order=SG_order)        
                if normalize == 0:
                    imagesmooth = norm_avg 
                if normalize == 2: 
                    imagesmooth = np.average(np.asarray(processing_frames),axis=0)#average of SG_bin_frames (SG_bin_frame > bin_frame)
                    imagesmooth = imagesmooth*(bin_frame+1) #because the binned frame will have the summed intensity of 41 frames
                if normalize == 3:
                    imagesmooth = np.average(np.asarray(processing_frames),axis=0)#average of SG_bin_frames (SG_bin_frame > bin_frame)
                    imagesmooth = imagesmooth*(bin_frame+1) #because the binned frame will have the summed intensity of 41 frames
                    imagesmooth = sgolay2d(imagesmooth, window_size=SG_window, order=SG_order) # apply SG

                CC = crosscor(imagesmooth.shape,mask=ii,normalization="symavg")#"regular")#
                rr=CC(imagesmooth,imagesmooth)

                #average contrast over all ROIs
                contrast_temp = []
                for j in np.arange(location_peak-3,location_peak+3): #only select around the peak area
                    try:
                        a=rr[j]
                #         p=r_[5,5]+np.unravel_index(np.argmax(a[5:-5,5:-5], axis=None), (10,10))
                        p=np.unravel_index(np.argmax(a, axis=None), a.shape) #position of max
                    #             print(a.shape)
                        ax=(a[p[0]-1,p[1]]+a[p[0]+1,p[1]]-2*a[p[0],p[1]])/2.
                        dx=(a[p[0]-1,p[1]]-a[p[0]+1,p[1]])/4./ax
                        ay=(a[p[0],p[1]-1]+a[p[0],p[1]+1]-2*a[p[0],p[1]])/2.
                        dy=(a[p[0],p[1]-1]-a[p[0],p[1]+1])/4./ay
                        cy=a[p[0],p[1]]-ay*dy*dy
                        cx=a[p[0],p[1]]-ax*dx*dx
                #         res2[j,(imgno-beg)//bin_frame,:]=r_[(cx+cy)/2.,p[0]+dx,p[1]+dy,np.sqrt(-cx/2/ax),np.sqrt(-cy/2/ay)]

                        contrast_temp_append = ((cx+cy)/2.)-1
                        contrast_temp.append(contrast_temp_append)
                    except:
                        pass
            contrast = np.average(np.asarray(contrast_temp))
            s_contrasts.append(contrast)
            SG_windows.append(SG_window)
            SG_window = 0 
    
            #Calculate a cross correlation between two adjacent frames to see the contrast of the image. 
            #binning frames to create another binned frame, which is next to "image"
            image2 = None #reset the image
            center_index = len(processing_frames)//2  + bin_frame
            select_range_to_bin = processing_frames[center_index- bin_frame//2:center_index+ bin_frame//2+1]
            image2 = np.sum(np.asarray(select_range_to_bin),axis=0) #now summed
            
            CC = crosscor(image.shape,mask=ii,normalization="symavg")#"regular")#
            rr=CC(image,image2)

            #average contrast over all ROIs
            contrast_temp = []
            for j in np.arange(location_peak-3,location_peak+3): #only select around the peak area
                try:
                    a=rr[j]
            #         p=r_[5,5]+np.unravel_index(np.argmax(a[5:-5,5:-5], axis=None), (10,10))
                    p=np.unravel_index(np.argmax(a, axis=None), a.shape) #position of max
                #             print(a.shape)
                    ax=(a[p[0]-1,p[1]]+a[p[0]+1,p[1]]-2*a[p[0],p[1]])/2.
                    dx=(a[p[0]-1,p[1]]-a[p[0]+1,p[1]])/4./ax
                    ay=(a[p[0],p[1]-1]+a[p[0],p[1]+1]-2*a[p[0],p[1]])/2.
                    dy=(a[p[0],p[1]-1]-a[p[0],p[1]+1])/4./ay
                    cy=a[p[0],p[1]]-ay*dy*dy
                    cx=a[p[0],p[1]]-ax*dx*dx
            #         res2[j,(imgno-beg)//bin_frame,:]=r_[(cx+cy)/2.,p[0]+dx,p[1]+dy,np.sqrt(-cx/2/ax),np.sqrt(-cy/2/ay)]

                    contrast_temp_append = ((cx+cy)/2.)-1
                    contrast_temp.append(contrast_temp_append)
                except:
                    pass
            contrast = np.average(np.asarray(contrast_temp))
            contrasts.append(contrast)
        else: #below SGwindow >0 so we are not going to do cross-correlation to find SG-window size automatically
            if normalize == 1:
                imagesmooth = sgolay2d(image, window_size=SG_window, order=SG_order)        
            if normalize == 0:
                imagesmooth = norm_avg 
            if normalize == 2: 
                imagesmooth = np.average(np.asarray(processing_frames),axis=0)#average of SG_bin_frames (SG_bin_frame > bin_frame)
                imagesmooth = imagesmooth*(bin_frame+1) #because the binned frame will have the summed intensity of 41 frames
            if normalize == 3:
                imagesmooth = np.average(np.asarray(processing_frames),axis=0)#average of SG_bin_frames (SG_bin_frame > bin_frame)
                imagesmooth = imagesmooth*(bin_frame+1) #because the binned frame will have the summed intensity of 41 frames
                imagesmooth = sgolay2d(imagesmooth, window_size=SG_window, order=SG_order) # apply SG

        
        
        
        #Work on the spline fit so that in dynamical ROI, we can move subpixel resolution    
        if discrete_rel_pos == False:
            if dynamical_ROI_relative_positions != np.zeros(10):
                from scipy import interpolate as sci_interpolate
                img_shape = np.shape(image)
                x = np.arange(0, img_shape[1], 1)
                y = np.arange(0, img_shape[0], 1)
                xx, yy = np.meshgrid(x, y)
                img_raw = image[yy,xx] #np.sin(xx**2+yy**2)
                img_smooth = imagesmooth[yy,xx]
                raw_f = sci_interpolate.interp2d(x, y, img_raw, kind='quintic')
                smooth_f = sci_interpolate.interp2d(x, y, img_smooth, kind='quintic')
                del image,imagesmooth        

        #Going through each ROI
        #we are going to put dynamical_ROI_relative_positions, where it is assumed that the relative positions are in Y direction
        #movement of speckles in X direction, even if they are in unison, are probably due to other complicated causes.
        #we are chaning Y positions where the angle of the sample moves slightly due to heater 
        Inor = []
        if moving_ref_frame == False:
            centery_temp = centery
        if len(dynamical_ROI_relative_positions) != 10: #if you add reasonable relative positions for dynamic ROI, put them here. 
            if discrete_rel_pos: #if you want the relative positions as integers
                centery = centery + np.int(np.round(dynamical_ROI_relative_positions[(frame-startframe)//bin_frame],decimals=0))
            else:
                centery = centery + dynamical_ROI_relative_positions[(frame-startframe)//bin_frame]
#             print(centery)
        for i in range(len(CenterXs)):
            if discrete_rel_pos: #if you want the relative positions as integers
                centerx = CenterXs[i]
                ROI = image[centery-halflength_y:centery+halflength_y+1,centerx-halflength_x:centerx+halflength_x+1]
                ROI_nor = imagesmooth[centery-halflength_y:centery+halflength_y+1,centerx-halflength_x:centerx+halflength_x+1]
                ROI_normalized = ROI/ROI_nor
                Inor.append(ROI_normalized)
            else: # the relative pos movements are sub pixel and spline fit is necssary
                centerx = CenterXs[i]
#                 ROI = image[centery-halflength_y:centery+halflength_y+1,centerx-halflength_x:centerx+halflength_x+1]
#                 ROI_nor = imagesmooth[centery-halflength_y:centery+halflength_y+1,centerx-halflength_x:centerx+halflength_x+1]
#                 ROI_normalized = ROI/ROI_nor
#                 Inor.append(ROI_normalized)
                
                xnew = np.arange(centerx-halflength_x, centerx+halflength_x+1, 1)
                ynew = np.arange(centery-halflength_y, centery+halflength_y+1, 1)
#                 print(centery)
                ROI = raw_f(xnew, ynew)
                ROI_nor = smooth_f(xnew, ynew)
                ROI_normalized = ROI/ROI_nor
                Inor.append(ROI_normalized)
                
                
        Inor_all.append(Inor)
        Inor = None
        image = None
        imagesmooth = None
        if moving_ref_frame == False:
            centery = centery_temp
        img_raw = None
        img_smooth = None
        raw_f = None
        smooth_f = None
        
        del image,imagesmooth,img_raw, img_smooth, raw_f, smooth_f 

    # current date and time
    now = datetime.now()
    # timestamp = datetime.timestamp(now)
    now = np.str_(now)
    now = now.replace('-','')
    now = now.replace(':','')
    now = now.replace(' ','')
    index = now.index('.')
    now = now[:index]


#     now concatanate metadata onto the Inor_all 
    if discrete_rel_pos:
#         try:
        metadata = np.zeros(np.shape(Inor_all)[0])
#         except:
#             metadata = np.zeros(len(Inor_all))
    else:
        metadata = np.zeros(len(Inor_all))

    #change boxsize tuple into a float... for example (21,21) will be 21.21
    boxsize = np.float(np.str(boxsize[0]) + '.' + np.str(boxsize[1]))

    metadata[0],metadata[1],metadata[2],metadata[3],metadata[4],metadata[5], metadata[6], metadata[7],  metadata[8],  metadata[9],metadata[10],metadata[11]   = startframe,endframe, centerx_start, centerx_end, centery,boxsize, bin_frame, normalize, SG_window, SG_order,  SG_bin_frame, interpolate



#     Inor_all = np.concatenate(np.asarray(Inor_all),metadata.T, axis =1)
    Inor_all.append(metadata.tolist())



    try:
        uid
    except NameError:
#         print("well, uid WASN'T defined! saving in primary folder")
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+'Inormalized_all_'+  '_' + now + '.npy'
        f2 = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+'CenterXs_'+  '_' + now + '.npy'
    else:
#         print("uid was defined. Saving in respective folder")
        directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + np.str_(np.round(convert_pixel_to_q(centerx_start,centery)[7],decimals=3)) + '/')
        if directory_exists == False:
            os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=3)) + '/')
        fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/'+ np.str_(np.round(convert_pixel_to_q(centerx_start,centery)[7],decimals=3)) + '/'+'Inormalized_all_'+ np.str_(np.round(convert_pixel_to_q(centerx_start,centery)[7],decimals=3)) + '_' + np.str(centerx_start) + '_' + np.str(centerx_end) + '_'  + now + '.npy'
        fp2 = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/'+ np.str_(np.round(convert_pixel_to_q(centerx_start,centery)[7],decimals=3)) + '/'+'CenterXs_'+ np.str_(np.round(convert_pixel_to_q(centerx_start,centery)[7],decimals=3)) + '_' + np.str(centerx_start) + '_' + np.str(centerx_end) + '_' + now + '.npy'

    np.save(fp,np.asarray(Inor_all))
    np.save(fp2,np.asarray(CenterXs))
    
    #combine contrasts and s_contrasts
    all_contrasts = [np.asarray(contrasts),np.asarray(s_contrasts)]

    return(np.asarray(CenterXs),np.asarray(Inor_all),all_contrasts,SG_windows)

def plot_contrasts_SG_windows(Windows, contrasts, starttime,endtime,binframe):
    times = np.arange(starttime,endtime,binframe)*datatakingtime
    list_to_plot = [Windows, contrasts[0],contrasts[1],contrasts[0]/contrasts[1]]
    list_ylabels = ['SG window size', 'Contrast', 'Contrast', 'Arbitrary Unit']
    list_titles = ['Optimized SG windows', 'Contrast of raw images \n computed from cross-correlation','Contrast of a smoothened image \n computed from cross-correlation', 'Contrast of raw image \n divided by contrast of smoothened image']
    fig, ax = plt.subplots(nrows=2,ncols=2, figsize = (10,6) )
    track_index = 0
    for row in ax:
        for col in row:
            col.plot(times, list_to_plot[track_index])
            col.set_title(list_titles[track_index],fontsize=12)
            col.set_xlabel(r'Time [s]',fontsize=12)
            col.set_ylabel(list_ylabels[track_index],fontsize=12)
            track_index +=1
    fig.tight_layout()
    plt.show()

    
def TTCF(CenterXs,Inor_all):
    TTCFs = None
    TTCFs = []
    Inor_all = Inor_all.tolist()
    metadata = Inor_all[-1]
    Inor_all = Inor_all[:-1]
    timelength = len(Inor_all) #better to use this then time total to contruct TTCF
    #convert Inor_all to array so that it is easier to select the same ROI over time
    Inor_all = np.asarray(Inor_all)
    
    #unpack the metadata
    startframe,endframe, centerx_start, centerx_end, centery,boxsize, bin_frame, normalize, SG_window, SG_order,  SG_bin_frame, interpolate = metadata[0],metadata[1],metadata[2],metadata[3],metadata[4],metadata[5], metadata[6], metadata[7],  metadata[8],  metadata[9], metadata[10], metadata[11] 
    
    TTCF = None
    TTCF2 = None
    Inormalized_flat = None 
    
#     print('now compute TT for each ROI')
    for i in tqdm(range(len(CenterXs))):
        centerx = CenterXs[i]
        Inormalized = Inor_all[:,i]
        timetotal = np.int(endframe) - np.int(startframe)
        TTCF = np.zeros((timelength,timelength)) #np.zeros((timetotal,timetotal)) #note timelength from Inor_all is used
        TTCF2 = np.zeros((timelength,timelength)) # np.zeros((timetotal,timetotal))
        Inormalized_flat = Inormalized.ravel()
        
        #multiplying intensity matrices and averaing over ensemble
#         print('computing ROI ', np.str(i))
        numberof_pixels = np.shape(Inormalized)
        for j in tqdm(range(np.shape(Inormalized)[1]*np.shape(Inormalized)[2])): #np.shape(Inormalized)[1] and [2] are dimensions of each ROI, [0] represents how many frames there are
            #select intensity for a single pixel over a range of time
            I_vals = Inormalized_flat[np.int_(np.linspace(j+1,(np.shape(Inormalized)[0]-1)*(np.shape(Inormalized)[1]*np.shape(Inormalized)[2])+j+1,np.shape(Inormalized)[0])-1)] # here there shoule be np.shape(Inormalized)[0] amount of data points while skipping by all_pixels present in each ROI
            #convert the intensity into an array
            TTCF_temp = np.asarray(I_vals)
            #multiply I vs t with itself and add them over (pixels)
            TTCF += np.outer(TTCF_temp, TTCF_temp.T)
            if j == 2:
                TTCF2 = np.outer(TTCF_temp, TTCF_temp)
        
        TTCF = TTCF/(np.shape(Inormalized)[1]*np.shape(Inormalized)[2]) #divide by the number of pixels
        
        #now replace the diagnoal with its adjacent because the correlation is abnormally high
        diagonal_line = TTCF.diagonal(offset = 1)
        diagonal_line = diagonal_line.tolist()
        diagonal_line = diagonal_line + [0] #add an element so that it has the same length as diagonal
        diagonal_line = np.asarray(diagonal_line)
        np.fill_diagonal(TTCF,diagonal_line)

        
        #now concatanate metadata onto the TTCF 
        metadata = np.zeros(np.shape(TTCF)[0])
        metadata[0],metadata[1],metadata[2],metadata[3],metadata[4],metadata[5],metadata[6], metadata[7],  metadata[8], metadata[9], metadata[10] = startframe,endframe, CenterXs[i], centery,boxsize, bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate
#         TTCF = np.concatenate((TTCF,metadata.T), axis =1)
        TTCF = TTCF.tolist()
        TTCF.append(metadata)
    
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
#             print("well, uid WASN'T defined! saving in primary folder")
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+'TTCF_'+  + np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3))+'_'+ now + '.npy'
        else:
#             print("uid was defined. Saving in respective folder")
            directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=3)) + '/')
            if directory_exists == False:
                os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=3)) + '/')
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/'+ np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=3)) + '/'+'TTCF_'+ np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=3)) + '_' + now + '.npy'

        np.save(fp, TTCF)
        TTCFs.append(TTCF)

    
    return(TTCF2, TTCFs)

def plotTTCF(TTCF,manual_name = '', trim = [], color_min = 1 , color_max = 1.15, color_map=plt.cm.jet, save = 0):
    #choose cmap_albula
    metadata = TTCF[-1]
    TTCF = TTCF[:-1]
    TTCF = np.asarray(TTCF) #convert back to array
    if trim:
        TTCF = TTCF[trim[0]:trim[1],trim[2]:trim[3]]
    
    #unfold the metadata
    startframe,endframe, centerx, centery,boxsize, bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate = metadata[0],metadata[1],metadata[2],metadata[3],metadata[4],metadata[5],metadata[6], metadata[7],  metadata[8], metadata[9], metadata[10]
    
    plt.figure(figsize=(10,10))
    plt.rcParams.update({'xtick.labelsize': 20})

    plt.rcParams.update({'ytick.labelsize': 20})

    plt.rcParams.update({'font.size': 12})

    mpl.rcParams.update({'figure.autolayout': True})
    plt.rcParams.update({'axes.titlesize': 20})

    plt.rcParams.update({'axes.labelsize': 20})

    plt.rcParams.update({'lines.linewidth': 2})

    plt.rcParams.update({'lines.markersize': 6})

    # mpl.rcParams.update({'legend.fontsize': 13})

    plt.rcParams.update({'legend.fontsize': 14})
    dimension = np.shape(TTCF)[0]
    ax = plt.gca()
    im = ax.imshow(TTCF, vmin = color_min , vmax = color_max, origin = 'lower',cmap=color_map, extent = [(startframe)*datatakingtime,((startframe)*datatakingtime)+((dimension)*datatakingtime*(bin_frame)),(startframe)*datatakingtime,((startframe)*datatakingtime)+((dimension)*datatakingtime*(bin_frame))], interpolation = 'none')#cmap=cmap_albula) 
    ax.set_ylabel(r'$t_2$ [s]')
    ax.set_xlabel(r'$t_1$ [s]')
    # create an axes on the right side of ax. The width of cax will be 5%
    # of ax and the padding between cax and ax will be fixed at 0.05 inch.
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", size="5%", pad=0.05)
    plt.colorbar(im, cax=cax)
    
    ax.set_title('Two-time correlation at ' + r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + r' $nm^{-1}$, '+ r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3)) + r' $nm^{-1}$')
    
    #split the boxsize back into tuple
    boxsize = np.str_(boxsize)
    boxsize = boxsize.split('.')
    boxsize = boxsize[0] + ', ' + boxsize[1]
    if save == 0:
        plt.suptitle('Center X = {0}'.format(centerx) +
                     ', Center Y = {0}'.format(centery) +
                     ', ROI Box Size = {0}'.format(boxsize) +
                     ' \n startframe = {0}'.format(startframe) + 
                     ', endframe = {0}'.format(endframe)+
                     ', Bin frame = {0}'.format(bin_frame)+
                     ', Interpolater = {0}'.format(interpolate) +
                     '\n SG_order = {0}'.format(SG_order) +
                     ', SG_box_size = {0}'.format(SG_window)+
                     ', SG Bin frame = {0}'.format(SG_bin_frame)+
                     ', Normalization Method = {0}'.format(normalize), fontsize=12)
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
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name +'TTCF_plot_'+ np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3))+'_'+ now + '.pdf'
        else:
            print("uid was defined. Saving in respective folder")
            directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
            if directory_exists == False:
                os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
            fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/' +'TTCF_plot_'+  np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3))+'_'+  now + '.pdf'

        plt.savefig(fp,bbox_inches='tight')
    plt.show()

def loadTTC(TT_plots, Min_val, Max_val, offset_num,save_switch): #the better way here is to generate images first and store and show them fast
    file_path = dirs + TT_plots_list[TT_plots_list.index(TT_plots)]
    print('looking at:'+file_path)
    
    TTCF = np.load(file_path)
    plotTTCF(TTCF, color_min = Min_val , color_max = Max_val, color_map=plt.cm.jet, save = save_switch )
    
    #get the metadata
    metadata = TTCF[-1]
    TTCF = TTCF[:-1]
    TTCF = np.asarray(TTCF) #convert back to array
    
    #unfold the metadata
    startframe,endframe, centerx, centery,boxsize, bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate = metadata[0],metadata[1],metadata[2],metadata[3],metadata[4],metadata[5],metadata[6], metadata[7],  metadata[8], metadata[9], metadata[10]
    
    TT = TTCF # I used TT below
    
    #determine the size of the TT calculations 
    sizeofTTC = np.shape(TT)[0]
    biggest_offset = (sizeofTTC-1)
    
    flippeddata = np.flip(TT,0) # flip the matrix for slicing ... but the way I cut through the diagonal has to be consistent with this.
    
    #create Vals numpy, which includes beta, tau, n and baselines and their standard errors and metadata and big T times. 10 columns
    '''
    now it reads from top right to bottom left .. but it ensures better fitting ... we will
    reverse the generated time values later to get the right order where T increases
    '''    
    
    #remember which fits fail and delete them at last 
    offset_num = -offset_num/(datatakingtime*bin_frame) #back to pixel scale
    #slice through the TT plots 
    try:

        data = flippeddata.diagonal(offset = np.int(offset_num), axis1=1, axis2=0)
        big_T = ((biggest_offset-offset_num))*datatakingtime*bin_frame + (startframe)*datatakingtime
        sizeofdata = np.shape(data)[0]
        startslice = np.int((sizeofdata//2)) #np.int((sizeofdata//2)-((allowance_width/datatakingtime)//2)) #allowance width is in seconds

        #processed data
        delta_t = np.arange(-len(data)//2,len(data)//2)[startslice:]*datatakingtime*bin_frame*2  #time should be in seconds before fitting times two because we are using alternative coordinate axis
        data = data[startslice:]
    except IndexError:
        exit()
    plt.figure()#figsize=(9,9))
    plt.scatter(delta_t,data, marker = '.', color = 'black',  label= r'T = {0}'.format(np.int(big_T)) + ' s') #multiplied by 0.1*1e3/1e3 which mean it is per 1000 seconds
    plt.ylabel(r'C($t_1$, $t_2$) = C(T, $\Delta t$)')
    plt.xlabel(r'$\Delta t$ = ($t_2$ - $t_1$) (s)')
#     plt.xscale('log')
#     plt.xlim(1,1100)
    plt.legend(loc='best')
    plt.title(r'Perpendicular cut through the TTCF, where T = $t_1$ + $t_2$')#, fontsize=8)
    plt.grid(False)
    plt.show()        
        
        
#     filelocation = dirs + TT_plots_list[TT_plots_list.index(TT_plots)]
#     TT = np.load(filelocation)
#     plotTTCF(TT, color_min = Min_val , color_max = Max_val, color_map=plt.cm.jet, save = save_switch )
#     bin_frame = TT[-1][5]
#     sizeofTTC = np.shape(TT)[0]
#     biggest_offset = ((sizeofTTC*np.sqrt(2))//2)
#     big_T = ((biggest_offset-offset_num)/np.sqrt(2))*datatakingtime*bin_frame*2
#     #slice through the TT plots
#     offset_num = offset_num/(datatakingtime*bin_frame) #back to pixel scale
#     flippeddata = np.flip(TT,0)
    
#     plt.figure()#figsize=(9,9))
#     data = flippeddata.diagonal(offset = np.int(offset_num), axis1=1, axis2=0)
#     plt.scatter(np.arange(-len(data)//2,len(data)//2)[1:]*datatakingtime*bin_frame,data[1:], marker = '.', color = 'black',  label= r'T = {0}'.format(np.int(big_T)) + ' s') #multiplied by 0.1*1e3/1e3 which mean it is per 1000 seconds
#     plt.ylabel(r'C($t_1$, $t_2$) = C(T, $\Delta t$)')
#     plt.xlabel(r'$\Delta t$ = ($t_2$ - $t_1$) (s)')
# #     plt.xscale('log')
# #     plt.xlim(1,1100)
#     plt.legend(loc='best')
#     plt.title(r'Perpendicular cut through the TTCF, where T = $t_1$ + $t_2$')#, fontsize=8)
#     plt.grid(False)
#     plt.show()
############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxx TTC calculations end  xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ##########
############################################################################################################

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx TTC analysis start   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ##########
############################################################################################################

def slicethrough_TTC_and_extract_parameters(file_path, trim = [], skip_frame = 1, allowance_width=200, fitinitials = [0.15,100,1.0,1.0], fitbounds = ([-1,-np.inf,-3.0,1.00],[1,np.inf,3.0,1.03]), plot_fits = [], savefits = True):
    '''
    bigT is the unprocessed T (t1+t2)/2 ...
    starting time is the start of bombardment 
    graph will have label T = age since the start of bombardment. 
    
    provide the file_paht for the computed TTCF..npy. it has metada in it
    allowance_width gets rid of the noise at lage delta t.. see the diagonal cut plots in TT
    
    !!allowance_width is in seconds. See TTCF plot for a good estimate of it 
    
    returns contrast Beta, tau, n, and baseline 
    
    if you use skip_frame, the binning or averaging along the diagonal will happen when slicing. It should give better statistics. 
    
    '''
    TTCF = np.load(file_path)
    
    #get the metadata
    metadata = TTCF[-1]
    TTCF = TTCF[:-1]
    TTCF = np.asarray(TTCF) #convert back to array
    
    if trim:
        TTCF = TTCF[trim[0]:trim[1],trim[2]:trim[3]]
    
    #unfold the metadata
    startframe,endframe, centerx, centery,boxsize, bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate = metadata[0],metadata[1],metadata[2],metadata[3],metadata[4],metadata[5],metadata[6], metadata[7],  metadata[8], metadata[9], metadata[10]
    
    TT = TTCF # I used TT below
    
    #determine the size of the TT calculations 
    sizeofTTC = np.shape(TT)[0]
    biggest_offset = (sizeofTTC-sizeofTTC//10)
    
    flippeddata = np.flip(TT,0) # flip the matrix for slicing ... but the way I cut through the diagonal has to be consistent with this.
    
    #create Vals numpy, which includes beta, tau, n and baselines and their standard errors and metadata and big T times. 10 columns
    offset_vals = np.arange(-biggest_offset,biggest_offset, skip_frame) 
    '''
    now it reads from top right to bottom left .. but it ensures better fitting ... we will
    reverse the generated time values later to get the right order where T increases
    
    example: plot_fits = np.linspace(1700,150,num=10).tolist()+[1E5])
    '''
    Vals = np.zeros((len(offset_vals),10))
    
    
    #remember which fits fail and delete them at last 
    pos = len(offset_vals)-1
    
    
    if plot_fits:
        f1 = plt.figure()
        ax1 = f1.add_subplot(111)
    
    for offset_num in tqdm(offset_vals):
    
        #slice through the TT plots 
        try:
            
            data = flippeddata.diagonal(offset = np.int(offset_num), axis1=1, axis2=0)
            if skip_frame != 1:#> 1000:
                sumtimes = 1
                for i in np.arange((skip_frame//2)*(-1),(skip_frame//2),1):
                    if i != 0:
                        tempdata = flippeddata.diagonal(offset = np.int(offset_num+i), axis1=1, axis2=0)
                        if (np.shape(tempdata)[0]-np.shape(data)[0])%2 != 0:
                            chop = 1
                        else:
                            chop = (np.shape(tempdata)[0]-np.shape(data)[0])//2
#                         print(np.shape(tempdata),np.shape(data), chop)
                        if chop%2 == 0:
                                if chop > 0:
#                                     print(np.shape(tempdata),np.shape(data), chop, 'greater than 0')
#                                     print("before", data)#np.shape(np.asarray(data)), np.shape(np.asarray(tempdata[chop:-chop])))
                                    data = np.asarray(data) + np.asarray(tempdata[chop:-chop])
#                                     print("after",data)
                                    sumtimes += 1
                                if chop < 0:
#                                     print(np.shape(tempdata),np.shape(data), chop, 'less than 0')
                                    chop = np.abs(chop)
                                    datatemp = data
                                    temp2 = np.asarray(datatemp/sumtimes)
                                    temp2[chop:-chop] = np.asarray(tempdata)
#                                     print('data',np.shape(data),data[:10])
#                                     print('temp2',np.shape(temp2),temp2[:10])
#                                     print('choppedtemp2',np.shape(temp2[chop:-chop]))
#                                     print('tempdata',np.shape(tempdata),tempdata[:10])
                                    data = np.asarray(data) + np.asarray(temp2)
#                                     print("after",np.shape(data),data[:10])
                                    sumtimes += 1
#                 print('before avg', data )
#                 print(sumtimes)
                data = data/sumtimes #do averaging
#                 print('after avg', data )

#             big_T = ((biggest_offset-offset_num))*datatakingtime*bin_frame + (startframe)*datatakingtime 
            
        
            if trim:
                big_T = (biggest_offset+1-(offset_num)+1)*datatakingtime*bin_frame + (startframe)*datatakingtime*bin_frame + trim[0]*datatakingtime*bin_frame

            else:
                big_T = (biggest_offset+1-(offset_num)+1)*datatakingtime*bin_frame + (startframe)*datatakingtime*bin_frame
#                 
        
        
        
            sizeofdata = np.shape(data)[0]
            startslice = np.int((sizeofdata//2))+0 #np.int((sizeofdata//2)-((allowance_width/datatakingtime)//2)) #allowance width is in seconds
            if allowance_width == 0:
                endslice = None #(biggest_offset+1-(offset_num)+1)-100 #   #None #allowance width is in seconds      
            else:
                endslice = np.int((sizeofdata//2)+((allowance_width/(datatakingtime*bin_frame))//2))  #allowance width is in seconds      

            #processed data
            try:
                delta_t = np.arange(-len(data)//2,len(data)//2)[startslice:endslice]*datatakingtime*bin_frame*2 #time should be in seconds before fitting . times two because we are using alternative coordiate axis
                data = data[startslice:endslice]
            except IndexError:
                delta_t = np.arange(-len(data)//2,len(data)//2)[startslice:]*datatakingtime*bin_frame*2 #time should be in seconds before fitting
                data = data[startslice:]

            params, params_covariance = optimize.curve_fit(KWW, delta_t, data,bounds=fitbounds, p0=fitinitials) 
            p0_previous = params # use the previous params. 
            
              
            #plot the fit
            if plot_fits:
                bigT_now = big_T//2
#                 print('qparallelnow  plog2', bigT_now, plot_fits,np.abs(plot_fits[0]-bigT_now) < 10)
                if  np.abs(plot_fits[0]-bigT_now) < 10:# and plot_fits[0] > bigT_now-10 :
                    plot_fits.pop(0)
                    print('qparallelnow chosen  plog2', bigT_now, plot_fits )
                    ax1.scatter(delta_t, data, marker = '.')#+ r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$' ) 
                    newdelta_t = np.linspace(delta_t[0],delta_t[-1],num = 1000)
                    ax1.plot(newdelta_t, KWW(newdelta_t, params[0],params[1],params[2],params[3]),label= 'T = {0} s'.format(big_T/2))
#             if plot_fits == 1:
#                 #sliceandfit TTCF
#                 plt.figure()
                
                
#             #     observed_values= g2
#             #     expected_values = KWW(timeinfo, params[0],params[1],params[2])
#             #     test = sp.stats.chisquare(observed_values, f_exp=expected_values)
#             #     p_value = test[1]
#             #     chi_square = test[0]
#                 plt.rcParams.update({'font.size': 14})
#                 plt.ylabel(r'C($t_1$, $t_2$) = C(T, $\Delta t$)')
#                 plt.xlabel(r'$\Delta t$ = ($t_2$ - $t_1$) (s)')
#                 plt.legend(loc='best')
#                 plt.suptitle('startframe = {0}'.format(startframe) + 
#                              ', endframe = {0}'.format(endframe) + 
#                              ', center Y = {0}'.format(centery) + 
#                              ', center X = {0}'.format(centerx) + 
#                              ', Box Size = {0}'.format(boxsize) +
#                              '\n ' +
#                              r'$B(q)$ = {0}'.format(np.round(params[3],decimals=2)) +
#                              r',$\beta (q)$ = {0}'.format(np.round(params[0],decimals=2)) 
#                              + r', $\tau (q)$ = {0}'.format(np.round(params[1],decimals=2)) +
#                              r', $n (q)$ = {0}'.format(np.round(params[2]),decimals=2), fontsize=12)

#                 plt.rcParams['axes.facecolor'] = 'white'
#                 plt.grid(False)
#                 plt.box(True)
#                 plt.xscale('log')
#                 plt.show()

            Vals[pos,0] = big_T/2  # divided by two because biggestoffset - smallestoffset is 2t1 or 2t2 #time values should be in second. this is big T
            Vals[pos,1] = params[0]
            Vals[pos,2] = params[1]
            Vals[pos,3] = params[2]
            Vals[pos,4] = params[3]
            std_deviations = np.sqrt(np.diag(params_covariance))
            Vals[pos,5] = std_deviations[0]
            Vals[pos,6] = std_deviations[1]
            Vals[pos,7] = std_deviations[2]
            Vals[pos,8] = std_deviations[3]
            pos -= 1 
        except RuntimeError:
            Vals = np.delete(Vals,pos,0)
#             print('deleted', np.str(offset_num))
            pos -= 1 
        except ValueError: 
            Vals = np.delete(Vals,pos,0)
#             print('deleted', np.str(offset_num))
            pos -= 1 
            
    #finally store the last column for storing metadata
    if np.shape(Vals)[0] > 12: 
        Vals[:11,9] = startframe,endframe, centerx, centery,boxsize, bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate
    if plot_fits:
        ax1.set_ylabel(r'C($t_1$, $t_2$) = C(T, $\Delta t$)')
        ax1.set_xlabel(r'$\Delta t$ = ($t_2$ - $t_1$) (s)')
        ax1.legend(bbox_to_anchor=(1.05, 1.0), loc='upper left')
    #             plt.rcParams.update({'font.size': 14})
    #             plt.rcParams['axes.facecolor'] = 'white'
    #             plt.grid(False)
    #             plt.box(True)
#         ax1.set_xscale('log')
#         ax1.set_xlim(delta_t[0]+delta_t[1]-delta_t[0],delta_t[-1]+delta_t[-1]-delta_t[-2])
        plt.title(r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$, ' + r'$q′_z$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$')
        plt.tight_layout()
        if savefits == True and plot_fits:
            plt.savefig(uid+'TTCFslices_and_fits.pdf')
        plt.show()
    return(Vals)

def slicethrough_TTC_early_time(file_path,indexofslice = 1, trim = [], skip_frame = 1, startT=400, bigTparameters = [0,200,20], offset =1, allowance_width=200,p0_initials =[0.15,0.01,10,1.0], fitbounds = ([-1,-np.inf,10,1.00],[1,np.inf,10+1E-6,1.03]), figsizes = [7.2,4.45],plot_fits = 0):
    '''
    startT in second, allow_width in second
    trim in seconds ... converted by datataking time
    p0_initials arwe beta, R, T, B 
    
    if you are trmming the early static state by using trim = []. startT should be set to 0 
    if you are not trmming, set starT to relavent start T
    
    allowance width spans across both side of the diagonal ridge 
    
    '''
    big_Tstart, big_Tend, big_Tstep = bigTparameters[0],bigTparameters[1],bigTparameters[2]
    big_T_temp = -300 ###
    #plot the fit
    bigTs = []
    yvalsforbombardment = []
    
    TTCF = np.load(file_path)
    
    #get the metadata
    metadata = TTCF[-1]
    TTCF = TTCF[:-1]
    TTCF = np.asarray(TTCF) #convert back to array
    
    
    
    #unfold the metadata
    startframe,endframe, centerx, centery,boxsize, bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate = metadata[0],metadata[1],metadata[2],metadata[3],metadata[4],metadata[5],metadata[6], metadata[7],  metadata[8], metadata[9], metadata[10]
    
    if trim:
        trim = (np.asarray(trim)/(datatakingtime*bin_frame)).astype('int')
        TTCF = TTCF[trim[0]:trim[1],trim[2]:trim[3]]
    
    TT = TTCF # I used TT below
    
    #determine the size of the TT calculations 
    sizeofTTC = np.shape(TT)[0]
    biggest_offset = (sizeofTTC-1)
    
    flippeddata = np.flip(TT,0) # flip the matrix for slicing ... but the way I cut through the diagonal has to be consistent with this.
    p0_previous = p0_initials 
    
    #create Vals numpy, which includes beta, tau, n and baselines and their standard errors and metadata and big T times. 10 columns
    offset_vals = np.arange(-biggest_offset,biggest_offset+1, skip_frame) [::-1] # most positive offset is at T =0 
    '''
    now it reads from bottom left to top right
    '''
    Vals = np.zeros((len(offset_vals),10))
    
    
    #remember which fits fail and delete them at last 
    pos = len(offset_vals)-1
    print(offset_vals)
    if plot_fits == 1:
        plt.figure(figsize=(figsizes[0],figsizes[0]))
    count = 0
    for offset_num in tqdm(offset_vals):
        
#         print('here')
    
        #slice through the TT plots 
        try:
            
            data = flippeddata.diagonal(offset = np.int(offset_num), axis1=1, axis2=0)
            if skip_frame != 1:#> 1000:
                sumtimes = 1
                for i in np.arange((skip_frame//2)*(-1),(skip_frame//2),1):
                    if i != 0:
                        tempdata = flippeddata.diagonal(offset = np.int(offset_num+i), axis1=1, axis2=0)
                        if (np.shape(tempdata)[0]-np.shape(data)[0])%2 != 0:
                            chop = 1
                        else:
                            chop = (np.shape(tempdata)[0]-np.shape(data)[0])//2
#                         print(np.shape(tempdata),np.shape(data), chop)
                        if chop%2 == 0:
                                if chop > 0:
#                                     print(np.shape(tempdata),np.shape(data), chop, 'greater than 0')
#                                     print("before", data)#np.shape(np.asarray(data)), np.shape(np.asarray(tempdata[chop:-chop])))
                                    data = np.asarray(data) + np.asarray(tempdata[chop:-chop])
#                                     print("after",data)
                                    sumtimes += 1
                                if chop < 0:
#                                     print(np.shape(tempdata),np.shape(data), chop, 'less than 0')
                                    chop = np.abs(chop)
                                    datatemp = data
                                    temp2 = np.asarray(datatemp/sumtimes)
                                    temp2[chop:-chop] = np.asarray(tempdata)
#                                     print('data',np.shape(data),data[:10])
#                                     print('temp2',np.shape(temp2),temp2[:10])
#                                     print('choppedtemp2',np.shape(temp2[chop:-chop]))
#                                     print('tempdata',np.shape(tempdata),tempdata[:10])
                                    data = np.asarray(data) + np.asarray(temp2)
#                                     print("after",np.shape(data),data[:10])
                                    sumtimes += 1
#                 print('before avg', data )
#                 print(sumtimes)
                data = data/sumtimes #do averaging
#                 print('after avg', data )
            if trim:
                big_T = (biggest_offset+1-(offset_num)+1)*datatakingtime*bin_frame + (startframe)*datatakingtime*bin_frame + trim[0]*datatakingtime*bin_frame
        
#                 big_T = (((biggest_offset-offset_num))*datatakingtime*bin_frame + (startframe)*datatakingtime + (trim[0]*datatakingtime))/np.sqrt(2) *2
            else:
                big_T = (biggest_offset+1-(offset_num)+1)*datatakingtime*bin_frame + (startframe)*datatakingtime*bin_frame
#                 big_T = (((biggest_offset-offset_num))*datatakingtime*bin_frame + (startframe)*datatakingtime*bin_frame )/np.sqrt(2) *2
            
            #change the big_T inside the fitting parameters
            p0_initials[2] = big_T - startT
            big_T = p0_initials[2]
            fitbounds[0][2] = big_T - startT
            fitbounds[1][2] = big_T+1E-6 - startT
            
           
            sizeofdata = np.shape(data)[0]
            startslice = np.int((sizeofdata//2))+1 #np.int((sizeofdata//2)-((allowance_width/datatakingtime)//2)) #allowance width is in seconds
            if allowance_width == 0:
                endslice = None #allowance width is in seconds      
            else:
                endslice = np.int((sizeofdata//2)+((allowance_width/(datatakingtime*bin_frame))//2))  #allowance width is in seconds      

            #processed data
            try:
                delta_t = np.arange(-len(data)//2,len(data)//2)[startslice:endslice]*datatakingtime*bin_frame*2 #time should be in seconds before fitting  ... times 2 because it is diagonal ..it is not geometrical ... also big T is determined by offset, not through pythagorean theorem 
                data = data[startslice:endslice]
            except IndexError:
                delta_t = np.arange(-len(data)//2,len(data)//2)[startslice:]*datatakingtime*bin_frame*2 #time should be in seconds before fitting
                data = data[startslice+2:]
                
                
                
#             linth4TTCF(deltaT,beta,R,T,B)
#             params, params_covariance = optimize.curve_fit(linth4TTCF, delta_t, data,bounds=fitbounds, p0=p0_previous) 
#             p0_previous = params # use the previous params. 




            
            
#             print(big_T)
#             print(np.shape(data))
            
            if np.shape(data)[0] > 5: 
            
                if plot_fits == 1 and big_T > big_Tstart and big_T < big_Tend :
                    if big_T > big_T_temp + big_Tstep :
                        print(big_T, big_T_temp)

                        if big_T >= 1:
                            plt.scatter(delta_t, data+offset*count, marker = '.', label= 'T = {0}'.format(big_T) + ' s ') 
                            delta_t_at_Tzero, theindex = find_nearest(delta_t,[big_T])
    #                         print(delta_t_at_Tzero[0], theindex[0])

                            bigTs.append(big_T)
                            yvalsforbombardment.append(offset*count+1)

    #                         xvalforSTART_bombardment.append(delta_t_at_Tzero[0])
    #                         tempdata= data+1*count
    #                         yvalforSTART_bombardment.append(tempdata[theindex[0]])
    #                     print(delta_t)


    #                         params, params_covariance = optimize.curve_fit(linth4TTCF, delta_t[:theindex[0]], data[:theindex[0]],bounds=fitbounds, p0=p0_previous) 



    #                         plt.plot(delta_t[:theindex[0]], linth4TTCF(delta_t[:theindex[0]], params[0],params[1],params[2],params[3])[::-1]+0.51*(count)) 
#                             if big_T >= 5 and theindex[0]+1 < big_T:
                            plt.plot(delta_t[:theindex[0]], linth4TTCF(delta_t[:theindex[0]], p0_previous[0],p0_previous[1],p0_previous[2],p0_previous[3])+offset*(count))
                            #### take out a specific slice across the TTCF to do the linear theory fit 
                            if count == indexofslice:
                                toprint1, toprint2 = np.asarray(delta_t),np.asarray(data)
                                print('here is the index {0}'.format(theindex[0]), ' big T is {0}'.format(big_T))
                                thebigT,Timeendindex = big_T, theindex
                
                            big_T_temp = big_T
                            count += 1
                            lastdelta_t = delta_t[-1]
                
                            #fitting might not be successful
                            params, params_covariance = optimize.curve_fit(linth4TTCF, delta_t[:theindex[0]], data[:theindex[0]],bounds=fitbounds, p0=p0_previous) 
                            print(p0_previous)
                    #     observed_values= g2
                    #     expected_values = KWW(timeinfo, params[0],params[1],params[2])
                    #     test = sp.stats.chisquare(observed_values, f_exp=expected_values)
                    #     p_value = test[1]
                    #     chi_square = test[0]
                            
                            
                            Vals[pos,0] = big_T/2 #time values should be in second. this is big T
                            Vals[pos,1] = params[0]
                            Vals[pos,2] = params[1]
                            Vals[pos,3] = params[2]
                            Vals[pos,4] = params[3]
                            std_deviations = np.sqrt(np.diag(params_covariance))
                            Vals[pos,5] = std_deviations[0]
                            Vals[pos,6] = std_deviations[1]
                            Vals[pos,7] = std_deviations[2]
                            Vals[pos,8] = std_deviations[3]
                            pos -= 1 
                            
            else:
                Vals = np.delete(Vals,pos,0)
    #             print('deleted', np.str(offset_num))
                pos -= 1 
            
        except RuntimeError:
            Vals = np.delete(Vals,pos,0)
#             print('deleted', np.str(offset_num))
            pos -= 1 
        except ValueError: 
            Vals = np.delete(Vals,pos,0)
#             print('deleted', np.str(offset_num))
            pos -= 1 
            
    #finally store the last column for storing metadata
    if np.shape(Vals)[0] > 12: 
        Vals[:11,9] = startframe,endframe, centerx, centery,boxsize, bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate
    #add one more point for good display
#     bigTs.append((bigTs[-1]+(bigTs[-1]-bigTs[-2])))
#     yvalsforbombardment.append((yvalsforbombardment[-1]+(yvalsforbombardment[-1]-yvalsforbombardment[-2])))


#     plt.fill_between(bigTs,yvalsforbombardment,y2 = yvalsforbombardment[0], alpha=0.1,color='b', linewidth = 0)
#     plt.fill_betweenx(np.asarray(yvalsforbombardment),x1=np.asarray(bigTs)[-1]+2.960324E-1,x2 = lastdelta_t, linewidth = 0, alpha=0.1,color='b')
    
    plt.fill_betweenx(yvalsforbombardment,bigTs,x2 = lastdelta_t, linewidth = 0, alpha=0.1,color='b')
    
    
    
#     print(bigTs,yvalsforbombardment)
    plt.rcParams.update({'font.size': 14})
    plt.ylabel(r'C($t_1$, $t_2$) = C(T, $\Delta t$)')
    plt.xlabel(r'$\Delta t$ = ($t_2$ - $t_1$) (s)')
    plt.legend(loc='best')
#     plt.suptitle('startframe = {0}'.format(startframe) + 
#                  ', endframe = {0}'.format(endframe) + 
#                  ', center Y = {0}'.format(centery) + 
#                  ', center X = {0}'.format(centerx) + 
#                  ', Box Size = {0}'.format(boxsize) +
#                  '\n ' +
#                  r'$B(q)$ = {0}'.format(np.round(params[3],decimals=2)) +
#                  r',$\beta (q)$ = {0}'.format(np.round(params[0],decimals=2)) 
#                  + r', $R (q)$ = {0}'.format(np.round(params[1],decimals=2)) +
#                  r', $T$ = {0}'.format(np.round(params[2]),decimals=2), fontsize=12)
    plt.title(r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$, ' + r'$q′_z$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$',fontsize = 14 )

    plt.rcParams['axes.facecolor'] = 'white'
#                 plt.ylim(1.05,1.25)
    plt.grid(False)
    plt.box(True)
#                 plt.xscale('log')
    plt.savefig('Kr_earlytime_LinTh.pdf',bbox_inches='tight')
    plt.show()
    return(toprint1, toprint2,thebigT,Timeendindex)
#     return(Vals)



def plot_exponents_vs_T(exponents, manual_name = '', beta = 1, tau = 1, n = 1, baseline = 1, xlims = [], save = 0, startTime = 0):
    """ this is plotting exponents vs q// .. not to use for exponents vs big T """
    startmarkdone = [True,True,True,True]
    if type(exponents) is list:
        if beta == 1:
            f1 = plt.figure()
            ax1 = f1.add_subplot(111)
        if tau == 1:
            f2 = plt.figure()
            ax2 = f2.add_subplot(111)

        if n == 1:
            f3 = plt.figure()
            ax3 = f3.add_subplot(111)
        if baseline == 1:
            f4 = plt.figure()
            ax4 = f4.add_subplot(111)
        for i in tqdm(range(len(exponents))):
    #         f1 = plt.figure()
    #         ax1 = f1.add_subplot(111)
    #         f2 = plt.figure()
    #         ax2 = f2.add_subplot(111)
    #         ax1.plot(np.arange(1,10))
    #         ax2.plot(np.arange(10,20))
    #         ax1.plot(np.arange(1,10)+1)
            # plt.show()
            
            exponents_extracted = exponents[i] 
#             print(exponents_extracted)
            Vals, STD, metadata = exponents_extracted[:,0:5], exponents_extracted[:,5:9], exponents_extracted[:,9]
            startframe,endframe, centerx, centery,boxsize,bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate = metadata[:11]


            Vals_copy = Vals
            to_delete = []



#             #delete the data where beta's STD is too high
#             for i in range(np.int(len(Vals_copy[:,1]))):
#                 if (np.abs(Vals_copy[:,2][i]))/3 < np.abs(STD[:,1][i]):
#                     to_delete.append(i)
             #delete the data where beta's STD is too high
#             for i in range(np.int(len(Vals_copy[:,3]))):
#                 if (np.abs(Vals_copy[:,3][i])) > 1.9: #/3 < np.abs(STD[:,1][i]):
#                     to_delete.append(i)

            Vals_fixed = np.delete(Vals_copy,to_delete,0) #Vals
            STD_fixed = np.delete(STD,to_delete,0) #STD


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
                    fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name + np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3))+'_'+ now 
                else:
                    print("uid was defined. Saving in respective folder")
                    directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
                    if directory_exists == False:
                        os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
                    fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/'+ uid + np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3))+'_'+  now 


            if beta == 1:
#                 plt.figure()
#                 print(np.shape(Vals_fixed[:,0]), np.shape(Vals_fixed[:,1]))
                ax1.errorbar(Vals_fixed[:,0], Vals_fixed[:,1],yerr=STD_fixed[:,0], fmt = '.', elinewidth=0.5, capsize=2, label= r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$ ')
                ax1.set_ylabel(r'$\beta (q)$')
                ax1.set_xlabel(r'T ($t_1$ + $t_2$) [s] ')
                ax1.legend(loc='best')
            #     plt.ylim(0.05,0.25)
                if xlims:
                    ax1.set_xlim(xlims[0],xlims[1])
#                 plt.suptitle('Center X = {0}'.format(centerx) +
#                              ', Center Y = {0}'.format(centery) +
#                              ', ROI Box Size = {0}'.format(boxsize) +
#                              ' \n startframe = {0}'.format(startframe) + 
#                              ', endframe = {0}'.format(endframe)+
#                              ', Bin Frame = {0}'.format(bin_frame)+
#                              ', Interpolater = {0}'.format(interpolate) +
#                              '\n SG_order = {0}'.format(SG_order) +
#                              ', SG box size = {0}'.format(SG_window)+
#                              ', SG Bin frame = {0}'.format(SG_bin_frame)+
#                              ', Normalization Method = {0}'.format(normalize), fontsize=9)
                
                if startTime != 0 and startmarkdone[0]:
                    ax1.axvline(x=startTime, color='k', linestyle='--', label = 'Start of ion bombardment')
                    startmarkdone[0] = False
                if save == 1 and i == len(exponents) -1:
                    f1.savefig(fp+ '_beta.pdf',bbox_inches='tight')
#                 plt.show()

            if tau == 1:
#                 plt.figure()
                #change all the values before the bombardment start to zero
                list_of_indices = np.where(Vals_fixed[:,0]<startTime)
                Vals_fixed[:,2][list_of_indices] = 0
                STD_fixed[:,1][list_of_indices] = -1
                if startTime != 0 and startmarkdone[1]:
                    ax2.axvline(x=startTime, color='k', linestyle='--', label = 'Start of ion bombardment')
                    startmarkdone[1] = False
                #filter
                the_indexes = np.where(STD_fixed[:,1]<Vals_fixed[:,2])
                ax2.errorbar(Vals_fixed[:,0][the_indexes], Vals_fixed[:,2][the_indexes],yerr=STD_fixed[:,1][the_indexes],fmt ='--', elinewidth=1, capsize=2, label= r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$ ')
                ax2.set_ylabel(r'$\tau (q)$ [s]')
                ax2.set_xlabel(r'T ($t_1$ + $t_2$) [s] ')
                ax2.set_ylim(0,500)
                ax2.legend(loc='best')
                if xlims:
                    ax2.set_xlim(xlims[0],xlims[1])
        #         plt.suptitle('Center X = {0}'.format(centerx) +
        #                      ', Center Y = {0}'.format(centery) +
        #                      ', ROI Box Size = {0}'.format(boxsize) +
        #                      ' \n startframe = {0}'.format(startframe) + 
        #                      ', endframe = {0}'.format(endframe)+
        #                      ', Bin Frame = {0}'.format(bin_frame)+
        #                      ', Interpolater = {0}'.format(interpolate) +
        #                      '\n SG_order = {0}'.format(SG_order) +
        #                      ', SG box size = {0}'.format(SG_window)+
        #                      ', SG Bin frame = {0}'.format(SG_bin_frame)+
        #                      ', Normalization Method = {0}'.format(normalize), fontsize=9)
                
                
                if save == 1 and i == len(exponents) -1:
                    f2.savefig(fp+ '_tau.pdf',bbox_inches='tight')
#                 plt.show()

            if n == 1:
#                 plt.figure()
#                 ax3.errorbar(Vals_fixed[5:-15,0], Vals_fixed[5:-15,3],yerr=STD_fixed[5:-15,2], fmt='--', elinewidth=0.5, capsize=2, label= r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$ ')
                to_plot = Vals_fixed[5:-15,3]
                #filter to delete data that has uncertanties bigger than 0.3
                if np.abs(convert_pixel_to_q(centerx,centery))[1] < 0.10 or np.abs(convert_pixel_to_q(centerx,centery))[1] > 0.17: #  >.25high wavenumbers get special treatment and are allowed bigger uncertanity
                    the_indexes = np.where(STD_fixed[5:-15,2]<0.4)
                else:
                    the_indexes = np.where(STD_fixed[5:-15,2]<0.3)
                xvalforn, yvalforn = Vals_fixed[5:-15,0][the_indexes], Vals_fixed[5:-15,3][the_indexes]
                #filter to delete data points that are unreasonable
                the_indexes = np.where(yvalforn<1.9)
                xvalforn, yvalforn = xvalforn[the_indexes], yvalforn[the_indexes]
                
                #filter to delete data points before the bombardment started
                the_indexes = np.where(xvalforn>200)
                xvalforn, yvalforn = xvalforn[the_indexes], yvalforn[the_indexes] 
                if startTime != 0 and startmarkdone[2]:
                    ax3.axvline(x=startTime, color='k', linestyle='--', label = 'Start of ion bombardment')
                    startmarkdone[2] = False
                ax3.plot(xvalforn, yvalforn, '.-',label= r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$ ')
                ax3.set_ylabel(r'$n(q)$')
                ax3.set_xlabel(r'T ($t_1$ + $t_2$) [s] ')
                ax3.legend(loc='best')
                ax3.set_ylim(0.5,1.9)
                if xlims:
                    ax3.set_xlim(xlims[0],xlims[1])
        #         plt.suptitle('Center X = {0}'.format(centerx) +
        #                      ', Center Y = {0}'.format(centery) +
        #                      ', ROI Box Size = {0}'.format(boxsize) +
        #                      ' \n startframe = {0}'.format(startframe) + 
        #                      ', endframe = {0}'.format(endframe)+
        #                      ', Bin Frame = {0}'.format(bin_frame)+
        #                      ', Interpolater = {0}'.format(interpolate) +
        #                      '\n SG_order = {0}'.format(SG_order) +
        #                      ', SG box size = {0}'.format(SG_window)+
        #                      ', SG Bin frame = {0}'.format(SG_bin_frame)+
        #                      ', Normalization Method = {0}'.format(normalize), fontsize=9)
                
                
                if save == 1 and i == len(exponents) -1:
                    f3.savefig(fp+ '_n.pdf',bbox_inches='tight')
#                 plt.show()

            if baseline == 1:
#                 plt.figure()
                ax4.errorbar(Vals_fixed[:,0], Vals_fixed[:,4],yerr=STD_fixed[:,3], fmt='.', elinewidth=0.5, capsize=2, label= r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$ ')
                ax4.set_ylabel(r'Baseline')
                ax4.set_xlabel(r'T ($t_1$ + $t_2$) [s] ')
                ax4.legend(loc='best')
#                 ax4.suptitle('Center X = {0}'.format(centerx) +
#                              ', Center Y = {0}'.format(centery) +
#                              ', ROI Box Size = {0}'.format(boxsize) +
#                              ' \n startframe = {0}'.format(startframe) + 
#                              ', endframe = {0}'.format(endframe)+
#                              ', Bin Frame = {0}'.format(bin_frame)+
#                              ', Interpolater = {0}'.format(interpolate) +
#                              '\n SG_order = {0}'.format(SG_order) +
#                              ', SG box size = {0}'.format(SG_window)+
#                              ', SG Bin frame = {0}'.format(SG_bin_frame)+
#                              ', Normalization Method = {0}'.format(normalize), fontsize=9)
                ax4.set_ylim(0.9,1.1)
                if xlims:
                    ax4.set_xlim(xlims[0],xlims[1])
                
                if startTime != 0 and startmarkdone[3]:
                    ax4.axvline(x=startTime, color='k', linestyle='--', label = 'Start of ion bombardment')
                    startmarkdone[3] = False
                if save == 1 and i == len(exponents) -1:
                    f4.savefig(fp+ '_baseline.pdf',bbox_inches='tight')
        
        plt.legend(loc='best')
        plt.show()  
        
    if type(exponents) is np.ndarray:
        Vals, STD, metadata = exponents[:,0:5], exponents[:,5:9], exponents[:,9]
        startframe,endframe, centerx, centery,boxsize,bin_frame, normalize, SG_window , SG_order ,  SG_bin_frame, interpolate = metadata[:11]


        Vals_copy = Vals
        to_delete = []



    #     #delete the data where beta's STD is too high
    #     for i in range(np.int(len(Vals_copy[:,1]))):
    #         if (np.abs(Vals_copy[:,2][i]))/3 < np.abs(STD[:,1][i]):
    #             to_delete.append(i)

        Vals_fixed = np.delete(Vals_copy,to_delete,0) #Vals
        STD_fixed = np.delete(STD,to_delete,0) #STD


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
                fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/'+ manual_name + np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3))+'_'+ now 
            else:
                print("uid was defined. Saving in respective folder")
                directory_exists = os.path.isdir('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
                if directory_exists == False:
                    os.makedirs('/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/')
                fp = '/Users/pmyint/Downloads/TIFF_ALD_99b99be0/Analysis_Peco/' + '%s_'%(uid) +'/'+ uid + np.str_(np.round(convert_pixel_to_q(centerx, centery)[1],decimals=3)) + '_' + np.str_(np.round(convert_pixel_to_q(centerx, centery)[7],decimals=3))+'_'+  now 


        if beta == 1:
            plt.figure()
            plt.errorbar(Vals_fixed[:,0], Vals_fixed[:,1],yerr=STD_fixed[:,0], color='black', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label=r'$\beta (q)$ at ' + r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$, '+ r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$')
            plt.ylabel(r'$\beta (q)$')
            plt.xlabel(r'T ($t_1$ + $t_2$) [s] ')
            plt.legend(loc='best')
        #     plt.ylim(0.05,0.25)
            plt.suptitle('Center X = {0}'.format(centerx) +
                         ', Center Y = {0}'.format(centery) +
                         ', ROI Box Size = {0}'.format(boxsize) +
                         ' \n startframe = {0}'.format(startframe) + 
                         ', endframe = {0}'.format(endframe)+
                         ', Bin Frame = {0}'.format(bin_frame)+
                         ', Interpolater = {0}'.format(interpolate) +
                         '\n SG_order = {0}'.format(SG_order) +
                         ', SG box size = {0}'.format(SG_window)+
                         ', SG Bin frame = {0}'.format(SG_bin_frame)+
                         ', Normalization Method = {0}'.format(normalize), fontsize=9)
            if save == 1:
                plt.savefig(fp+ '_beta.pdf',bbox_inches='tight')
            plt.show()

        if tau == 1:
            plt.figure()
            plt.errorbar(Vals_fixed[:,0], Vals_fixed[:,2],yerr=STD_fixed[:,1], color = 'red', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label=r'$\tau (q)$ at ' + r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$, '+ r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$')
            plt.ylabel(r'$\tau (q)$ [s]')
            plt.xlabel(r'T ($t_1$ + $t_2$) [s] ')
            plt.ylim(0,250)
            plt.legend(loc='best')
    #         plt.suptitle('Center X = {0}'.format(centerx) +
    #                      ', Center Y = {0}'.format(centery) +
    #                      ', ROI Box Size = {0}'.format(boxsize) +
    #                      ' \n startframe = {0}'.format(startframe) + 
    #                      ', endframe = {0}'.format(endframe)+
    #                      ', Bin Frame = {0}'.format(bin_frame)+
    #                      ', Interpolater = {0}'.format(interpolate) +
    #                      '\n SG_order = {0}'.format(SG_order) +
    #                      ', SG box size = {0}'.format(SG_window)+
    #                      ', SG Bin frame = {0}'.format(SG_bin_frame)+
    #                      ', Normalization Method = {0}'.format(normalize), fontsize=9)
            if save == 1:
                plt.savefig(fp+ '_tau.pdf',bbox_inches='tight')
            plt.show()

        if n == 1:
            plt.figure()
            plt.errorbar(Vals_fixed[5:-15,0], Vals_fixed[5:-15,3],yerr=STD_fixed[5:-15,2], color = 'blue', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label=r'$n(q)$ at '+ r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$, '+ r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$')
            plt.ylabel(r'$n(q)$')
            plt.xlabel(r'T ($t_1$ + $t_2$) [s] ')
            plt.legend(loc='best')
            plt.ylim(0.9,1.9)
    #         plt.suptitle('Center X = {0}'.format(centerx) +
    #                      ', Center Y = {0}'.format(centery) +
    #                      ', ROI Box Size = {0}'.format(boxsize) +
    #                      ' \n startframe = {0}'.format(startframe) + 
    #                      ', endframe = {0}'.format(endframe)+
    #                      ', Bin Frame = {0}'.format(bin_frame)+
    #                      ', Interpolater = {0}'.format(interpolate) +
    #                      '\n SG_order = {0}'.format(SG_order) +
    #                      ', SG box size = {0}'.format(SG_window)+
    #                      ', SG Bin frame = {0}'.format(SG_bin_frame)+
    #                      ', Normalization Method = {0}'.format(normalize), fontsize=9)
            if save == 1:
                plt.savefig(fp+ '_n.pdf',bbox_inches='tight')
            plt.show()

        if baseline == 1:
            plt.figure()
            plt.errorbar(Vals_fixed[:,0], Vals_fixed[:,4],yerr=STD_fixed[:,3], color = 'blue', fmt='.k', ecolor='lightgray', elinewidth=1, capsize=3, label='Baseline  at ' + r'$q_{//}$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[1],decimals=2)) + r' $nm^{-1}$, '+ r'$q_z′$ = ' + np.str_(np.round(convert_pixel_to_q(centerx,centery)[7],decimals=2)) + r' $nm^{-1}$')
            plt.ylabel(r'Baseline')
            plt.xlabel(r'T ($t_1$ + $t_2$) [s] ')
            plt.legend(loc='best')
            plt.suptitle('Center X = {0}'.format(centerx) +
                         ', Center Y = {0}'.format(centery) +
                         ', ROI Box Size = {0}'.format(boxsize) +
                         ' \n startframe = {0}'.format(startframe) + 
                         ', endframe = {0}'.format(endframe)+
                         ', Bin Frame = {0}'.format(bin_frame)+
                         ', Interpolater = {0}'.format(interpolate) +
                         '\n SG_order = {0}'.format(SG_order) +
                         ', SG box size = {0}'.format(SG_window)+
                         ', SG Bin frame = {0}'.format(SG_bin_frame)+
                         ', Normalization Method = {0}'.format(normalize), fontsize=9)
            plt.ylim(0.9,1.1)
            if save == 1:
                plt.savefig(fp+ '_baseline.pdf',bbox_inches='tight')
            plt.show()             







############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx TTC analysis end   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ##########
############################################################################################################

############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Speckle track start   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ##########
############################################################################################################

from numpy import r_

def see_speckle_along_a_frame(startY,endY,delY,startX,endX,delX,image,img2shift = 0,Xlogscale= False, normalize = 1, bin_frame = 1, SG_frame = 5, SG_order = 1, SG_box_size = 23, interpolater_on = False):
    """
    will call the functions in autocor.py. make sure to run this first. the g2 calculation in here is done by the code written
    Mark Sutton. It follows the multi-tau method
    
    skip => skips the frame of interest
    bin_num and bin_skip are for normalization 
    bin_num bin over specified number of frames for normalization and bin_skip determines if you want to do binning over the bin_num frames around the next frame or skip to another one. However you bin_num or bin_skip it. the normalzation will match with the t1 or t2 you are considering. This computational structer was conceived first before TTCF binning method. 
    
    skip and bin_skip are not in sync necessarily. 
    
    comparing to TTCF terms. 
    skip --> bin_frame
    bin_num --> SG_Bin_frame
    bin_skip ---> no equivalence in TTCF codes... it lets to skip frames while binning for SG to save computation time
    
    
    star, end del X Y represent the pixel numbers usually from 0 to 2160 in 4M eiger
    start frame and end frame should be within available frame range of the data
    Xlogscale will increase delX in log scale, the increment is determined by Xlog_num
    p0 are the initial parameters for the curefit to try while fitting for the exponents
    normalizewith_Sg = 0 usues the average frame of the entire start-end frame range
    normalizewith_Sg = 1 usues SG smootheing per frame with SG-order and SG-Box_size
    normalizewith_Sg = 2 does binning frames and averaging them and using only the localized average.
    normalizewith_Sg = 3 does binning frames and averaging them and using the smoothened frame of the average.
    
    """
#     numberofframes = endframe - startframe
    image = 0*image
    image += 1
    pixellist=np.where(image.ravel()==1)[0] #list of valid pixels
    #partition in bins using x,y (i,j) indices
    yy,xx=np.unravel_index(pixellist,image.shape)
    Xlogscale = False
    if Xlogscale:
        ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
        xlist=mkpartlist(np.round(np.geomspace(startX,endX,Xlog_num),decimals=0).astype(int)) # break points for bins in x
        print(xlist)
        
#         ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
#         xlist=mkpartlist(np.round(np.logspace(np.log10(startX), np.log10(endX), num=Xlog_num),decimals=0).astype(int)) # break points for bins in x
        print(xlist)
        
    else:
        ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
        xlist=mkpartlist(np.arange(startX,endX,delX)) # break points for bins in x
    plist,bind,xbind,ybind,noperbin,inpixellist,binlist,binxlist,binylist=\
        partition2d(pixellist,xx,xlist,yy,ylist)
    print(xx,yy, xlist, len(xlist))
    
    #how to look at partitions
    ii=0*image
    ii[np.unravel_index(plist,ii.shape)]=bind+1 #bind==0 is not a bin
    # Choose images to test correlations on.
#     img1=image
    img0 = image*ii
    
    
#     image_INTEREST = None
#     image_INTEREST = FD.rdframe(1005+img2shift)
#     for i in np.arange(1006+img2shift,1046+img2shift):
#         image_INTEREST += FD.rdframe(i)
#     img2= image_INTEREST#image - sgolay2d(img1, window_size=SG_box_size, order=SG_order)
    
#     if SG_order == 0:
#         img1 = image
#     else:
#         img1 = sgolay2d(image, window_size=SG_box_size, order=SG_order)
    
    
#     img2= sgolay2d(img1, window_size=SG_box_size, order=SG_order)
#     img1 = sgolay2d(img1, window_size=SG_box_size, order=SG_order)
#     img0 = img1*ii
    img1 = image
    CC = crosscor(img1.shape,mask=ii,normalization="symavg")
    rr=CC(img1,img1)
    
    
    idno = 0 
    #Show how this works for an example
    a=rr[idno] #speckle cross cor
    #solve by quadratic on 3 points
    # print(res[idno,:]) #fit results
    p=np.unravel_index(np.argmax(a, axis=None), a.shape) #position of max
    ax=(a[p[0]-1,p[1]]+a[p[0]+1,p[1]]-2*a[p[0],p[1]])/2.
    dx=(a[p[0]-1,p[1]]-a[p[0]+1,p[1]])/4./ax
    cx=a[p[0],p[1]]-ax*dx*dx
    ay=(a[p[0],p[1]-1]+a[p[0],p[1]+1]-2*a[p[0],p[1]])/2.
    dy=(a[p[0],p[1]-1]-a[p[0],p[1]+1])/4./ay
    cy=a[p[0],p[1]]-ay*dy*dy
    fwhmx=2*np.sqrt(-((cx+cy)/2.-1)/2./ax) 
    fwhmy=2*np.sqrt(-((cx+cy)/2.-1)/2./ay) 
    #print(ax,dx,cx,ay,dy,cy)
    dp=r_[dx,dy] #shif
    # print(p,p+dp)
    # print(a[p[0],p[1]],res[idno,1]+res[idno,6],cx,cy,(cx+cy)/2.)
    # print(res[idno,[4,5]],2*np.sqrt(-(cx-1)/ax/2.),2*np.sqrt(-(cy-1)/ay/2.),fwhmx,fwhmy)
    #Let's calulate for all of them
    res1=np.zeros((len(rr),5))
    idno=0
    for j in tqdm(range(len(rr))):
        a=rr[j]
        p=np.unravel_index(np.argmax(a, axis=None), a.shape) #position of max

    #     p=r_[5,5]+np.unravel_index(np.argmax(a[5:-5,5:-5], axis=None), (10,10))
        ax=(a[p[0]-1,p[1]]+a[p[0]+1,p[1]]-2.0*a[p[0],p[1]])/2.0
        dx=(a[p[0]-1,p[1]]-a[p[0]+1,p[1]])/4.0/ax
        cx=a[p[0],p[1]]-ax*dx*dx
        ay=(a[p[0],p[1]-1]+a[p[0],p[1]+1]-2*a[p[0],p[1]])/2.
        dy=(a[p[0],p[1]-1]-a[p[0],p[1]+1])/4./ay
        cy=a[p[0],p[1]]-ay*dy*dy
        res1[idno,:]=r_[(cx+cy)/2.,p[0]+dx,p[1]+dy,np.sqrt(-cx/2/ax),np.sqrt(-cy/2/ay)]
        idno += 1
        plt.figure()
        plt.imshow(a)
        plt.colorbar()
        plt.show()
        
        
    #the q// range
    x_pixel_range = np.linspace(startX,endX,len(res1[:,0]))
    Y_center = startY + delY//2
    q_parallel_range = []
    for i in range(len(x_pixel_range)):
        q_parallel_range.append(np.round(convert_pixel_to_q(x_pixel_range[i],Y_center)[1],decimals=3))
    q_parallel_range = np.asarray(q_parallel_range)
    
    #plot results beta,px,py,sx,sy
    #note how shift in speckle peaks occurs between these two images.
    fig,((ax0),(ax1),(ax2))=plt.subplots(3,1,figsize=(8,10))
    ax0.plot(q_parallel_range,res1[:,0]-1.0) #beta
    ax0.set_xlabel(r'$q_{||}$ [nm$^{-1}$]')
    ax0.set_ylabel(r'$\beta$')
#     ax0.set_ylim(0.1,0.5)
    
    ax1.plot(q_parallel_range,res1[:,1], label = 'x') #center in x
    ax1.set_xlabel(r'$q_{||}$ [nm$^{-1}$]')
    ax1.set_ylabel('Center')
    
    
    ax1.plot(q_parallel_range,res1[:,2], label = 'x') #center in y
   
    
    ax2.plot(q_parallel_range,res1[:,3], label = 'x') #width in x
    ax2.set_xlabel(r'$q_{||}$ [nm$^{-1}$]')
    ax2.set_ylabel('Width')
#     ax2.set_ylim(0,4)
    ax2.plot(q_parallel_range,res1[:,4], label = 'y') #width in y
    ax1.legend(loc='best')
    ax2.legend(loc='best')
    
#     ax0.plot(res2[j,:,0]-1.0,label="ROI No =%d"%j)#beta
#         ax0.set_xlabel('New binned frame number')
#         ax0.set_ylabel(r'\beta')

#         ax2.plot(res2[j,:,2],'.-',label="Center X =%d"%j)#center x
#         ax2.set_xlabel('New binned frame number')
#         ax2.set_ylabel('Center')

#         ax3.plot(res2[j,:,3]+j,label="Width X =%d"%j)#width x
#         ax3.set_xlabel('New binned frame number')
#         ax3.set_ylabel('Width')

#         ax2.plot(res2[j,:,1],label="Center Y =%d"%j)#center y


#         ax3.plot(res2[j,:,4]+j,label="Witdh Y =%d"%j)#width y

#         ax1.plot(res2[j,:,2]+j,'.-',label="center x with offset =%d"%j) #center x with offset
#         ax1.plot(res2[j,:,1]+j,'.-',label="center y with offset =%d"%j) #center y with offset
#         ax1.set_xlabel('New binned frame number')
#         ax1.set_ylabel(r'Width')

#         ax0.legend(loc='best')
#         ax2.legend(loc='best')

def track_speckle(startY,endY,delY,startX,endX,delX,startframe,endframe, refframe,skip_ROI_every___frame = 1, Xlogscale= False, normalize = 1, bin_frame = 1, SG_frame = 5, SG_order = 1, SG_box_size = 23, interpolater_on = False):
    """
    SG doesn't work here ? 
    
    ***multiple calculations of cross-correlation will only work for the same qzs. I.e. you can only have one row of ROI (multiple columns allowed). The labelling in the legend will fail, but the code should still work if you have multiple rows of ROIs ... just disable putting legends in this scenario 
    
    will call the functions in autocor.py. make sure to run this first. the g2 calculation in here is done by the code written
    Mark Sutton. It follows the multi-tau method
    
    skip => skips the frame of interest
    bin_num and bin_skip are for normalization 
    bin_num bin over specified number of frames for normalization and bin_skip determines if you want to do binning over the bin_num frames around the next frame or skip to another one. However you bin_num or bin_skip it. the normalzation will match with the t1 or t2 you are considering. This computational structer was conceived first before TTCF binning method. 
    
    skip and bin_skip are not in sync necessarily. 
    
    comparing to TTCF terms. 
    skip --> bin_frame
    bin_num --> SG_Bin_frame
    bin_skip ---> no equivalence in TTCF codes... it lets to skip frames while binning for SG to save computation time
    
    
    star, end del X Y represent the pixel numbers usually from 0 to 2160 in 4M eiger
    start frame and end frame should be within available frame range of the data
    Xlogscale will increase delX in log scale, the increment is determined by Xlog_num
    p0 are the initial parameters for the curefit to try while fitting for the exponents
    normalizewith_Sg = 0 usues the average frame of the entire start-end frame range
    normalizewith_Sg = 1 usues SG smootheing per frame with SG-order and SG-Box_size
    normalizewith_Sg = 2 does binning frames and averaging them and using only the localized average.
    normalizewith_Sg = 3 does binning frames and averaging them and using the smoothened frame of the average.
    
    """
#     numberofframes = endframe - startframe
    
    pixellist=np.where(mask.ravel()==1)[0] #list of valid pixels
    #partition in bins using x,y (i,j) indices
    yy,xx=np.unravel_index(pixellist,mask.shape)
    
    if Xlogscale:
        ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
        xlist=mkpartlist(np.round(np.geomspace(startX,endX,Xlog_num),decimals=0).astype(int)) # break points for bins in x
#         print(xlist)
        
#         ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
#         xlist=mkpartlist(np.round(np.logspace(np.log10(startX), np.log10(endX), num=Xlog_num),decimals=0).astype(int)) # break points for bins in x
#         print(xlist)
        
    else:
        ylist=mkpartlist(np.arange(startY,endY,delY)) # break points for bins in y
        xlist=mkpartlist(np.arange(startX,endX,delX)) # break points for bins in x
    plist,bind,xbind,ybind,noperbin,inpixellist,binlist,binxlist,binylist=\
        partition2d(pixellist,xx,xlist,yy,ylist)
    
    
    if refframe !=0:
        img1=FD.rdframe(refframe-bin_frame//2) #reference (in middle)
        for i in np.arange(refframe-bin_frame//2 +1, refframe+bin_frame//2+1):
            img1 += FD.rdframe(i)
        moving_reference = False
    else:
        refframe = startframe-bin_frame
        img1=FD.rdframe(refframe-bin_frame//2) #reference (in middle)
        for i in np.arange(refframe-bin_frame//2 +1, refframe+bin_frame//2+1):
            img1 += FD.rdframe(i)
        moving_reference = True
    
    
    #how to look at partitions
    image = img1
    
    
    ii=-1+0*img1 #make bind==0 stand out
    ii[np.unravel_index(plist,ii.shape)]=bind
#     ii=0*image
#     ii[np.unravel_index(plist,ii.shape)]=bind+1 #bind==0 is not a bin
    # Choose images to test correlations on.
    img1=image
    img0 = image*ii
    
    ii=-1+0*img1 #make bind==0 stand out
    ii[np.unravel_index(plist,ii.shape)]=bind
    
    beg=startframe
    no=endframe-startframe+1
    
    
    if (endX-startX) < delX*2:
        res2=np.zeros((1,(no-1)//bin_frame,5)) #len(rr)
        make_bigger_res2 = False
    else:
        make_bigger_res2 = True
    CC = crosscor(img1.shape,mask=ii+1,normalization="symavg")

    for imgno in tqdm(np.arange(beg,beg+no-1,bin_frame)):
        if moving_reference:
            img1=None
            img1=FD.rdframe(refframe-bin_frame//2) #reference (in middle)
            for i in np.arange(refframe-bin_frame//2 +1, refframe+bin_frame//2+1):
                img1 += FD.rdframe(i)
            refframe = imgno-bin_frame
        img2=FD.rdframe(imgno-bin_frame//2) #reference (in middle)
        for i in np.arange(imgno-bin_frame//2 +1, imgno+bin_frame//2+1):
            img2 += FD.rdframe(i)
            
    #     cx=a[p[0],p[1]]-ax*dx*dx
        rr=CC(img1,img2)
#         print(np.shape(rr))
        
        if make_bigger_res2 == True:
            res2 = np.zeros((len(rr),(no-1)//bin_frame,5))
            make_bigger_res2 = False
        
        if (endX-startX) < delX*2:
            for j in range(1):#len(rr)
                try:
                    a=rr#[j]
        #             p=r_[5,5]+np.unravel_index(np.argmax(a[5:-5,5:-5], axis=None), (10,10))
                    p=np.unravel_index(np.argmax(a, axis=None), a.shape) #position of max
        #             print(a.shape)
                    ax=(a[p[0]-1,p[1]]+a[p[0]+1,p[1]]-2*a[p[0],p[1]])/2.
                    dx=(a[p[0]-1,p[1]]-a[p[0]+1,p[1]])/4./ax
                    ay=(a[p[0],p[1]-1]+a[p[0],p[1]+1]-2*a[p[0],p[1]])/2.
                    dy=(a[p[0],p[1]-1]-a[p[0],p[1]+1])/4./ay
                    cy=a[p[0],p[1]]-ay*dy*dy
                    cx=a[p[0],p[1]]-ax*dx*dx
                    res2[j,(imgno-beg)//bin_frame,:]=r_[(cx+cy)/2.,p[0]+dx,p[1]+dy,np.sqrt(-cx/2/ax),np.sqrt(-cy/2/ay)]
                except:
                    pass
        else:
            for j in range(len(rr)):#len(rr)
                try:
                    a=rr[j]
        #             p=r_[5,5]+np.unravel_index(np.argmax(a[5:-5,5:-5], axis=None), (10,10))
                    p=np.unravel_index(np.argmax(a, axis=None), a.shape) #position of max
        #             print(a.shape)
                    ax=(a[p[0]-1,p[1]]+a[p[0]+1,p[1]]-2*a[p[0],p[1]])/2.
                    dx=(a[p[0]-1,p[1]]-a[p[0]+1,p[1]])/4./ax
                    ay=(a[p[0],p[1]-1]+a[p[0],p[1]+1]-2*a[p[0],p[1]])/2.
                    dy=(a[p[0],p[1]-1]-a[p[0],p[1]+1])/4./ay
                    cy=a[p[0],p[1]]-ay*dy*dy
                    cx=a[p[0],p[1]]-ax*dx*dx
                    res2[j,(imgno-beg)//bin_frame,:]=r_[(cx+cy)/2.,p[0]+dx,p[1]+dy,np.sqrt(-cx/2/ax),np.sqrt(-cy/2/ay)]
                except:
                    pass
            
            
    #the q// range
    x_pixel_range = np.linspace(startX,endX,len(res2[:,0,0]))
    Y_center = startY + delY//2
    q_parallel_range = []
    for i in range(len(x_pixel_range)):
        q_parallel_range.append(np.round(convert_pixel_to_q(x_pixel_range[i],Y_center)[1],decimals=3))
    q_parallel_range = np.asarray(q_parallel_range)
    print(q_parallel_range)
    
    #plot results beta,px,py,sx,sy
    time = np.linspace(startframe,endframe - bin_frame, num =len(res2[j,:,0]-1.0))*datatakingtime
    fig,((ax0,ax1),(ax2,ax3))=plt.subplots(2,2,figsize=(15,8))
    if (endX-startX) < delX*2:
        range_frame = range(1)
    else:
        range_frame = np.arange(0,len(rr),skip_ROI_every___frame)
    for j in range_frame:#len(rr) # here len(rr) or range_frame represents the number of ROIs
        '''
        slicing [j,:,i] means at ROI number # j, and i = 0 is beta, i = 2 is center x position, i  = 3 is speckle width x, 
        i = 1 is center y position, i =4 is with y
        '''
        ax0.plot(time,res2[j,:,0]-1.0,label=r"$q_{//}$ = %s"%q_parallel_range[j] + r' $nm^{-1}$')#beta
        ax0.set_xlabel('Time [s]')
        ax0.set_ylabel(r'$\beta$')

        ax2.plot(time,res2[j,:,2],'.-',label=r"Center X @ $q_{//}$ = %s"%q_parallel_range[j]+ r' $nm^{-1}$')#center x
        ax2.set_xlabel('Time [s]')
        ax2.set_ylabel('Center')

        ax3.plot(time, res2[j,:,3]+j,label=r"Width X with offset @ $q_{//}$ = %s"%q_parallel_range[j]+ r' $nm^{-1}$')#width x
        ax3.set_xlabel('Time [s]')
        ax3.set_ylabel('Width')

        ax2.plot(time,res2[j,:,1],label=r"Center Y @ $q_{//}$ = %s"%q_parallel_range[j]+ r' $nm^{-1}$')#center y


        ax3.plot(time,res2[j,:,4]+j,label=r"Witdh Y with offset @ $q_{//}$ = %s"%q_parallel_range[j]+ r' $nm^{-1}$')#width y

        ax1.plot(time,res2[j,:,2]+j,'.-',label=r"center X with offset @ $q_{//}$ = %s"%q_parallel_range[j]+ r' $nm^{-1}$') #center x with offset
        ax1.plot(time,res2[j,:,1]+j,'.-',label=r"center Y with offset @ $q_{//}$ = %s"%q_parallel_range[j]+ r' $nm^{-1}$') #center y with offset
        ax1.set_xlabel('Time [s]')
        ax1.set_ylabel('Center with offset')

        # Shrink current axis by 20%
        box = ax0.get_position()
        ax0.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax0.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        # Shrink current axis by 20%
        box = ax1.get_position()
        ax1.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax1.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        
        # Shrink current axis by 20%
        box = ax2.get_position()
        ax2.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax2.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        # Shrink current axis by 20%
        box = ax3.get_position()
        ax3.set_position([box.x0, box.y0, box.width * 0.9, box.height])
        ax3.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    plt.suptitle(r'Speckle track along $q_{z}^′$ = %s'%np.round(convert_pixel_to_q(x_pixel_range[i],Y_center)[7],decimals=3)+ r' $nm^{-1}$' )
    fig.tight_layout()
    plt.show()
    return(time, res2)


############################################################################################################
########   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxx Speckle track end   xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx ##########
############################################################################################################