"""
Last modified date 01/22/2019

@author: Peco
"""
'''The following code is use latex computer modern font for all plots'''
import matplotlib as mpl
mpl.rcParams.update(mpl.rcParamsDefault)
mpl.rcParams['font.family'] = "serif"
mpl.rcParams.update({'font.size': 12})
mpl.rcParams.update({'figure.figsize': (7.2,4.45)})

# from matplotlib import rc
# rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# ## for Palatino and other serif fonts use:
# #rc('font',**{'family':'serif','serif':['Palatino']})
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



# mpl.rcParams['font.serif'] = "cm"
# from matplotlib import rc
# # rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})
# # for Palatino and other serif fonts use:
# rc('font',**{'family':'serif','serif':['Palatino']})
# rc('text', usetex=True)


def long_wave_form(q,Sy,B):
    return -B*(q**4)-Sy*(q**2)


def Rwithorchard(q,Sy,A,h0):
    Q = q*h0
    return -Sy*(q**2)-A*((Q*(np.sinh(2*Q)-2*Q))/(1+2*(Q**2)+np.cosh(2*Q)))

def plotandfit_R(Rs_list,name_list = [], colors_select = [], manual_name = '', use_errors = True, long_wave=True, long_wave_restrict=False,long_wave_sy = -11.2, param_lw_bounds = ([-np.inf, 0],[np.inf, np.inf]), param_lw_restrict_bounds=([0],[np.inf]), orchard=True,orchard_restrict=False,orchard_sy= 0,orchard_h0=0, orchard_po = [300.03,0.2],orchard_bounds=([-np.inf,0,0],[np.inf,np.inf,np.inf]), param_orchard_restrict_bounds=([0,0],[np.inf,np.inf]),savefigure=0):    
    length_oflist = len(Rs_list)
    plt.figure()
    
    for j in range(length_oflist):
        Rs = Rs_list[j]
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

        R = R_copy
        
        
        #pick different colors for each error bar plot
        if len(colors_select) == 0:
            colors_select = plt.cm.rainbow(np.linspace(0, 1, len(Rs_list)))

        if len(name_list) == 0:
            plt.errorbar(R[:,0], R[:,1]*1E3,yerr=R[:,2]*1E3, fmt='.k', color='black', ecolor='lightgray', elinewidth=1, capsize=3, label= r'R($q_{||}$) at $q_z$ = ' + np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+np.absolute(ysliceend-yslicestart)//2)[2],decimals=3)) + ' ' + r'nm$^{-1}$' ) #
        else:
            plt.errorbar(R[:,0], R[:,1]*1E3,yerr=R[:,2]*1E3,color=colors_select[j], fmt='.', ecolor='lightgray', elinewidth=1, capsize=3,label= name_list[j] + r'R($q_{||}$) at $q_z$ = ' + np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+np.absolute(ysliceend-yslicestart)//2)[2],decimals=3)) + ' ' + r'nm$^{-1}$' ) #, 
            print('here')
    #     param_orchard_bounds=([-np.inf,0,0],[np.inf,np.inf,2])
    #     param_orchard_bounds_2=([0,0],[np.inf,100])

    #     param_lw_bounds = ([-np.inf, 0],[np.inf, np.inf])
    #     param_lw_bounds_2=([0],[np.inf])


        if use_errors:
            errors = R[:,2]
        else:
            errors = None
        param_orchard_bounds = orchard_bounds 





    #     param_bounds=([-np.inf,-np.inf,-np.inf],[np.inf,np.inf,np.inf])
    
    
        # now replace R[:,0] with a more finely spaced array for better plotting
        q_pll = np.linspace(R[:,0][0],R[:,0][-1],num=100)
        if long_wave:
            if long_wave_restrict:
                Sy = long_wave_sy #here setting things brute force
                params, params_covariance = optimize.curve_fit(lambda q, B: long_wave_form(q, B, Sy), R[:,0], R[:,1], sigma = errors, bounds= param_lw_restrict_bounds) 
                params_covariance_std_deviations = np.sqrt(np.diag(params_covariance))
                plt.plot(q_pll, long_wave_form(q_pll, params[0],Sy)*1E3,label=r'$- S_{y}q^{2} - Bq^{4}$ fit:' + '\n'+ r'$S_y$ = {0}'.format(np.round(Sy,decimals=2)) + r' nm$^2$/s, B = {0}'.format(np.round(params,decimals=2))+ r' nm$^4$/s.')
                text = ' \n Sy = {0}'.format(Sy) + ', B = {0}'.format(params)+ ' nm^2/s nm^4/s. 1std = {0}'.format(params_covariance_std_deviations)
            else:

                params, params_covariance = optimize.curve_fit(long_wave_form, R[:,0], R[:,1], sigma = errors, bounds = param_lw_bounds )#, p0=[1.4,21]) #lambda q, B: long_wave_form(q, B, Sy), R[:,0], R[:,1], p0=[2]) 
                params_covariance_std_deviations = np.sqrt(np.diag(params_covariance))
                plt.plot(q_pll, long_wave_form(q_pll, params[0],params[1])*1E3,label=r'$- S_{y}q^{2} - Bq^{4}$ fit:' + '\n'+ r'$S_y$ = {0} $\pm$ {1:.2f}'.format(np.round(params[0],decimals=2),params_covariance_std_deviations[0]) + r' nm$^2$/s, '+ '\n'+'B = {0} $\pm$ {1:.2f}'.format(np.round(params[1],decimals=2),params_covariance_std_deviations[1])+ r' nm$^4$/s') 
                text = ' \n [Sy, B] = {0}'.format(params)+ 'nm^2/s nm^4/s. 1std = {0}'.format(params_covariance_std_deviations)

        if orchard:
            if orchard_restrict:
                
                if orchard_h0 == 0 and orchard_sy !=0:
                    Sy = orchard_sy #here setting things brute force
                    paramsOrchard, params_covarianceOrchard = optimize.curve_fit(lambda q, A, h0: Rwithorchard(q,Sy,A,h0), R[:,0], R[:,1], sigma = errors, p0=orchard_po, bounds = param_orchard_restrict_bounds) 
                    params_covarianceOrchard_std_deviations = np.sqrt(np.diag(params_covarianceOrchard))
                    plt.plot(q_pll, Rwithorchard(q_pll, Sy, paramsOrchard[0], paramsOrchard[1])*1E3,label= r'Orchard Term fit: $S_y$ = {0}'.format(np.round(orchard_sy,decimals=2)) + r' nm$^2/s$' +', A = {0}'.format(np.round(paramsOrchard[0],decimals=2)) + r' s^${-1}$, $h_0$ =' + ' {0}'.format(np.round(paramsOrchard[1],decimals=2)) + ' nm')
                    text += '\n Sy = {0}'.format(orchard_sy) + ' , [A, h0] = {0}'.format(paramsOrchard) + 'nm^2/s  s^-1  nm. 1std = {0}'.format(params_covarianceOrchard_std_deviations)
                
                if orchard_h0 != 0 and orchard_sy ==0:
                    h0 = orchard_h0 #here setting things brute force
                    paramsOrchard, params_covarianceOrchard = optimize.curve_fit(lambda q, Sy, A: Rwithorchard(q,Sy,A,h0), R[:,0], R[:,1], sigma = errors, p0=orchard_po, bounds = param_orchard_restrict_bounds) 
                    params_covarianceOrchard_std_deviations = np.sqrt(np.diag(params_covarianceOrchard))
                    plt.plot(q_pll, Rwithorchard(q_pll, paramsOrchard[0], paramsOrchard[1],h0)*1E3,label= r'Orchard Term fit: $S_y$ = {0}'.format(np.round(paramsOrchard[0],decimals=2)) + r' $nm^2/s$' +', A = {0}'.format(np.round(paramsOrchard[1],decimals=2)) + r' $s^{-1}$, $h_0$ =' + ' {0}'.format(np.round(h0,decimals=2)) + ' nm')
                    text += '\n h0 = {0}'.format(h0) + ' , [Sy, A] = {0}'.format(paramsOrchard) + 'nm nm^2/s  s^-1. 1std = {0}'.format(params_covarianceOrchard_std_deviations)
                
                
            else:
                paramsOrchard, params_covarianceOrchard = optimize.curve_fit(Rwithorchard, R[:,0], R[:,1], sigma = errors, bounds = param_orchard_bounds) #, p0=[10.29, 10.84,0.51]
                params_covarianceOrchard_std_deviations = np.sqrt(np.diag(params_covarianceOrchard))
                plt.plot(q_pll, Rwithorchard(q_pll, paramsOrchard[0], paramsOrchard[1],paramsOrchard[2] )*1E3,label=r'Orchard Term fit: $S_y$ = {0}'.format(np.round(paramsOrchard[0],decimals=2)) + r' $nm^2/s$' +', A = {0}'.format(np.round(paramsOrchard[1],decimals=2)) + r' $s^{-1}$, $h_0$ =' + ' {0}'.format(np.round(paramsOrchard[2],decimals=2)) + ' nm')
                text += '\n [Sy, A, h0] = {0}'.format(paramsOrchard) + 'nm^2/s  s^-1  nm. 1std = {0}'.format(params_covarianceOrchard_std_deviations)

        #compute q_max
        q_max = np.sqrt(np.abs(params[0])/(2*params[1]))

        #label on the right positive region 
        plt.annotate(r'$q_{max}$ = ' + np.str_(np.round(q_max,decimals=2)) + r' nm$^{-1}$' +'\nR = ' + np.str_(np.round(long_wave_form(q_max,params[0],params[1])*1E3,decimals=2)) + r' $\times$ $10^{-3}$ s$^{-1}$',xy=(q_max, long_wave_form(q_max,params[0],params[1])*1E3), xycoords='data',xytext=(-80, -90), textcoords='offset points', fontsize=13,arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.1"))
        #label on the left -- negative region 
        plt.annotate(r'$q_{max}$ = ' + np.str_(np.round(-1*q_max,decimals=2)) + r' nm$^{-1}$' +'\nR = ' + np.str_(np.round(long_wave_form(q_max,params[0],params[1])*1E3,decimals=2)) + r' $\times$ $10^{-3}$ s$^{-1}$',xy=(-1*q_max, long_wave_form(q_max,params[0],params[1])*1E3), xycoords='data',xytext=(-50, -90), textcoords='offset points', fontsize=13,arrowprops=dict(arrowstyle="->", connectionstyle="arc3,rad=.1")) # = \sqrt{\frac{\left| S_y \right|}{2B}} 
        plt.ylabel(r'R($q_{||}$) $10^{-3}$ [s$^{-1}$]')
        plt.xlabel(r'$q_{||}$ [nm$^{-1}$]')
        
        if length_oflist > 1:
            plt.legend(loc='center left', bbox_to_anchor=(1, 0.5))
        else:
            plt.legend(loc='best')
        
        if length_oflist == 1:
            if long_wave==False and orchard==False:
                text = None
#             if chi_end_frame != 0:
#                 more_label = '\n chi_end_frame = {0}'.format(chi_end_frame) + ', increment = {0}'.format(increment)
#                 plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
#                              ', ysliceend = {0}'.format(ysliceend) + 
#                              ' xslicestart = {0}'.format(xslicestart) + 
#                              ', xsliceend = {0}'.format(xsliceend) + 
#                              ' delx = {0}'.format(delx) + 
#                              ' startframe = {0}'.format(startframe) + 
#                              ', endframe = {0}'.format(endframe) + text + more_label, fontsize=8, y = 1.1)
#             else:
#                 plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
#                              ', ysliceend = {0}'.format(ysliceend) + 
#                              ' xslicestart = {0}'.format(xslicestart) + 
#                              ', xsliceend = {0}'.format(xsliceend) + 
#                              ' delx = {0}'.format(delx) + 
#                              ' startframe = {0}'.format(startframe) + 
#                              ', endframe = {0}'.format(endframe) + text, fontsize=8, y = 1.1)
        from datetime import datetime
        # current date and time
        now = datetime.now()
        # timestamp = datetime.timestamp(now)
        now = np.str_(now)
        now = now.replace('-','')
        now = now.replace(':','')
        now = now.replace(' ','')
        index = now.index('.')
        now = now.replace('.','')
        now = now[:index+2]
    
    if len(name_list) >1:
        print("well, uid WASN'T defined! saving in primary folder")
        fp = '/home/pmyint/fits_plots/'+ manual_name +'_Rs_fitted_at_qz_'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3)) + '_' + now + '.pdf'
    else:
        try:
            uid
        except NameError:
            print("well, uid WASN'T defined! saving in primary folder")
            fp = '/home/pmyint/fits_plots/'+ manual_name +'_Rs_fitted_at_qz_'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3)) + '_' + now + '.pdf'
        else:
            print("uid was defined. Saving in respective folder")
            directory_exists = os.path.isdir('/home/pmyint/' + '%s_'%(uid) +'/')
            if directory_exists == False:
                os.makedirs('/home/pmyint/' + '%s_'%(uid) +'/')
            ion_angle = md['ion angle'].replace(' ', '_')
            name = md['sample'] +'_'+ ion_angle +'_' + md['gas'] + '_'+'Energy'+  md['beam voltage']  +'_%s_'%(uid)
            fp = '/home/pmyint/' + '%s_'%(uid) +'/' + name +'_Rs_fitted_at_qz_'+ np.str_(np.round(convert_pixel_to_q(xslicestart,yslicestart+(ysliceend-yslicestart)//2)[2],decimals=3))  + '_' + now + '.pdf'
    if savefigure == 1:
        plt.savefig(fp,bbox_inches='tight')
    plt.show()

def plotandfit_R_constraint(R,yslicestart,ysliceend,xslicestart,xsliceend,delx,startframe,endframe):
    """
    here we constarint one parameter with brute force while fitting Rs
    for example, we can also set h0 as a certain value 
    """
    Sy = 0 #here setting something as bruteforce
    params, params_covariance = optimize.curve_fit(lambda q, B: long_wave_form(q, B, Sy), R[:,0], R[:,1], p0=[2]) 
    paramsOrchard, params_covarianceOrchard = optimize.curve_fit(lambda q, A, h0: Rwithorchard(q,Sy,A,h0), R[:,0], R[:,1], p0=[0.203,3.4]) 
    plt.figure()
    plt.scatter(R[:,0], R[:,1]*1E3, color = 'black', marker = '.', label='data') #multiplied by 0.1*1e3/1e3 which mean it is per 1000 seconds
    plt.plot(R[:,0], Rwithorchard(R[:,0], Sy, paramsOrchard[0],paramsOrchard[1] )*1E3,label='Orchard Term fit')
    plt.plot(R[:,0], long_wave_form(R[:,0], params[0],Sy)*1E3,label=r'$- S_{y}q^{2} - Bq^{4}$ fit') 
    plt.ylabel(r'R($q_{y}$) ($10^{-3} s^{-1}$)')
    plt.xlabel(r'q ($nm^{-1}$)')
    plt.legend(loc='best')
    plt.suptitle('yslicestart = {0}'.format(yslicestart) + 
                 ' ysliceend = {0}'.format(ysliceend) + 
                 ' \n xslicestart = {0}'.format(xslicestart) + 
                 ' xsliceend = {0}'.format(xsliceend) + 
                 ' delx = {0}'.format(delx) + 
                 ' \n startframe = {0}'.format(startframe) + 
                 ' endframe = {0}'.format(endframe) + 
                 ' \n Sy = 0, B = {0}'.format(params)+ 
                 ' Sy = 0, A, h0 = {0}'.format(paramsOrchard), fontsize=8)

    plt.show()
    plt.savefig(DD.DD['filename'][-4] + '_R_fitted')