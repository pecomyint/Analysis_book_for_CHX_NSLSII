"""
@author: Mark

edited by peco 08/01/2019 and interpolation between modular gaps is implemented 
"""

#multiqtau in python.
import numpy as np
from tqdm import tqdm_notebook as tqdm
# def insertimg(FD,n,plist,bind,nobind,nobuf,nolev,norm=1, to_interpolate = False):
def insertimg(FD,n,plist,bind,nobind,nobuf,nolev,normalizewith_SG, norm_order, pixel_boxsize, to_interpolate = False):
    """insertimg(FD,n,plist,bind,nobind,nobuf,nolev,norm=1)
       read and insert image, n, into accumulator buffers,
       plist is list of pixels, bind their bins and
       nobind is the number in each bind.
       then increments appropriate correlators for new image.
       Use function rdframe(n) to interface to the image reading routines.
       norm is 1 or SG used to normalize (from smoothed image).
    """
    global buf,cur,cts, norm_avg, norm_bins, norm_bins_SG, startframe, skip_to_bin
    
    skip = skip_to_bin
    
    cur[0]=(1+cur[0])%nobuf
#     img=eliminate_module_beamstop(FD.rdframe(n))[1] #eliminate_module_beamsstop added by Peco
#     img=FD.rdframe(n)
    img = FD.rdframe(n-skip//2)
    for i in np.arange(n-skip//2+1, n+skip//2):
        img+=FD.rdframe(i)
        
    
    #interpolate the gaps -- Peco
    if to_interpolate == True:
        import pandas
        df = pandas.DataFrame(img)
        mask_df = pandas.DataFrame(mask)
        mask_df = mask_df.replace(0,np.nan)
        df = df*mask_df
        df = df.interpolate(method='linear', limit_direction='forward', limit = 15, axis=1) #along x
        df = df.interpolate(method='linear', limit_direction='backward', axis=0) #along y
        df = df.replace(np.nan,0)
        img = df.as_matrix()
    #interpolate the gaps end -- Peco    
    
    
    if normalizewith_SG == 0:
        norm_avg = norm_avg* (skip+1)
        norm = np.square(norm_avg.ravel()[plist]) #flatten and then pick the right pixels and then square them to normalize
        
#         norm = 1
    elif normalizewith_SG == 1:
        boxsize = pixel_boxsize
        smoothened = sgolay2d(img, window_size=pixel_boxsize, order=norm_order)
#         print('smoothened with ' + np.str_(boxsize)+ ' '+ np.str_(norm_order))
        smooth_img = np.asarray(smoothened)* (skip+1)
        norm = np.square(smooth_img.ravel()[plist]) #flatten and then pick the right pixels and then square them to normalize
    elif normalizewith_SG == 2:
        mod_img = norm_bins[n-startframe]* (skip+1)
        norm = np.square(mod_img.ravel()[plist]) #flatten and then pick the right pixels and then square them to normalize
    elif normalizewith_SG == 3:
        smooth_img = norm_bins_SG[n-startframe]* (skip+1)
        norm = np.square(smooth_img.ravel()[plist]) #flatten and then pick the right pixels and then square them to normalize
        
    buf[0,cur[0],:]= img.ravel()[plist]/norm
    process(0,cur[0],bind,nobind,nobuf)
    processing=1
    lev=1
    while(processing):
        if(cts[lev]):
            prev= (cur[lev-1]+nobuf-1)%nobuf
            cur[lev]= (cur[lev]+1)%nobuf
            buf[lev,cur[lev],:]=(buf[lev-1,prev,:]+buf[lev-1,cur[lev-1],:])/2.0
            cts[lev]=0
            process(lev,cur[lev],bind,nobind,nobuf)
            lev += 1
            processing= 1 if (lev<nolev) else 0
        else:
            cts[lev]=1
            processing=0

def process(lev,bufno,bind,nobind,nobuf):
    """ process(lev,bufno,bind,nobind,nobuf)
        The buffer bufno at level lev has just been added so update
        the correlation and intensity averages it affects. bind are
        bin indexes and nobind is histogram(bind)
    """
    global num,G,IAP,IAF
    num[lev]+=1
    imin=0 if lev==0 else nobuf//2

    for i in range(imin, min(num[lev],nobuf)):
        ptr=lev*nobuf//2+i
        delayno=(bufno-(i)+nobuf)%nobuf
        IP=buf[lev,delayno, :]
        IF=buf[lev,bufno,:]
        G[ptr,:]+= (np.bincount(bind,IF*IP)/nobind-  G[ptr,:])/(num[lev]-i)
        IAP[ptr,:] += (np.bincount(bind,IP)/nobind-IAP[ptr,:])/(num[lev]-i)
        IAF[ptr,:] += (np.bincount(bind,IF)/nobind-IAF[ptr,:])/(num[lev]-i)
    #print('proc: ', G[lev*nobuf//2+imin,:])

def autocor(FD,begframe,noframes,plist,bind, nobuf=8,nolev=8,skip=1,norm=1,bin_frames = 1,bin_frame_skip =  5, interpolater = False, norm_order = 1, pixel_boxsize = 23):
    """autocor(FD,begframe,noframes,plist,bind,nobuf=8,nolev=8,skip=1,norm=1)
        Function to drive the acquisition of data and to autocorrelate
        the data using multiple tau method. Uses variables in .des.i for info.
        begframe is frame number of first pattern to use.
        noframes is number of frames to process.
        plist is list of pixels to analyze (pixellist).
        bind is binning index for pixels.
        KEYWORD: norm=1 or SG  to normalize by Smooth image SG
                 nobuf=8 number of images at each level
                 nolev=8 number of levels (8 goes to 896 delays)
                 skip=1 analyze skip'th frame
    """
    global buf,cts,cur,G,IAP,IAF,num, norm_avg, norm_bins, norm_bins_SG, startframe, skip_to_bin
    
    skip_to_bin = skip
    
    startframe = begframe
    nobs=1+int(max(bind))
    buf=np.zeros((nolev,nobuf,plist.size),dtype=float)
    cts=np.zeros(nolev,dtype=int)
    cur=(nobuf-1)*np.ones(nolev,dtype=int)
    G=np.zeros([(nolev+1)*nobuf//2,nobs])
    IAP=np.zeros([(nolev+1)*nobuf//2,nobs])
    IAF=np.zeros([(nolev+1)*nobuf//2,nobs])
    num=np.zeros(nolev,dtype=int)
    nobind=np.bincount(bind)
    w=np.flatnonzero(nobind==0)
    if(w.size > 0):
        print("{:d} bins have no pixels in them.".format(w.size))
        nobind[w]=1 #prevent divide by 0s
#     print(begframe)
#     print(noframes)
    
    if norm == 0:
        temp_img =FD.rdframe(begframe)
        print('calculating average')
        for i in tqdm(np.arange(1,noframes)):
            temp_img += FD.rdframe(begframe+i)
        norm_avg = temp_img/noframes
    
    if norm == 2:
        norm_bins = []
        print('binning and averaging. binning every {0}'.format(bin_frame_skip))
        loop = np.arange(startframe,noframes+startframe,bin_frame_skip)
#         print(loop)
        for i in tqdm(loop):
            img = None
            img=FD.rdframe(i)
            for j in range(bin_frames):
                img+=FD.rdframe(i+j-bin_frames//2) #shifting the center so it bins equal amount of frames before and after i
            img = img - FD.rdframe(i) # to prevent double counting
            img = img/bin_frames # get the average
            for k in range(bin_frame_skip): # need to have len(norm_bins) == number of frames that are being considered
                norm_bins.append(img)
    if norm == 3:
        norm_bins_SG = []
        print('binning and averaging  and applying SG. binning every {0}'.format(bin_frame_skip))
        loop = np.arange(startframe,noframes+startframe,bin_frame_skip)
        for i in tqdm(loop):
            img = None
            img=FD.rdframe(i)
            for j in range(bin_frames):
                img+=FD.rdframe(i+j-bin_frames//2) #shifting the center so it bins equal amount of frames before and after i
            img = img - FD.rdframe(i) # to prevent double counting
            img = img/bin_frames # get the average
            smoothened = sgolay2d(img, window_size=pixel_boxsize, order=norm_order) # smoothen it
            for k in range(bin_frame_skip): # need to have len(norm_bins) == number of frames that are being considered
                norm_bins_SG.append(smoothened)

        
    print('calculating g2')            
    for n in tqdm(np.arange(0,noframes,skip)):
        insertimg(FD,begframe+n,plist,bind,nobind,nobuf,nolev, norm,norm_order, pixel_boxsize, to_interpolate = interpolater)
#         insertimg(FD,begframe+n,plist,bind,nobind,nobuf,nolev,norm=norm, to_interpolate = interpolater)

    gmax=np.flatnonzero(IAP[:,0]==0)[0]
    if(gmax.size==0):gmax=IAP[:,0].size
    else: gmax=gmax.min()
    a=np.arange(nobuf//2)
    #generate times for G.
    tt=np.append(a,np.multiply.outer(\
            2**np.arange(1+(gmax-nobuf//2)//(nobuf//2)),a+nobuf//2).ravel())
    return(tt[:gmax],G[:gmax,:]/(IAP[:gmax,:]*IAF[:gmax,:]))
