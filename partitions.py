"""
@author: Mark
"""
import numpy as np

def mkpartlist(lst,width=0):
    '''
    list = mkpartlist(values,width=0)
      Return a list suitable for partitioning and selecting bins.
      A bin is all pixels such that values[i]<=x<values[i+1] allowing
      for empty gaps between values.

      If width <=0 returns [values[0],values[1],values[1],values[2]...]
      else returns [values[0]-width/2,values[0]+width/2,...].

      Beware if width > spacing as this returns a valid list but probably
      not what you expect.

      NOTE: upper limit of values should not equal lower limit of next bin.
         This program offsets them by 1e-10.

         print(list.reshape((len(list)//2,2)).T) is a good way to see list.
    '''
    lst=np.asarray(lst)
    if(width <= 0.0):
        list=np.zeros(2*lst.size-2)
        list[0::2]=lst[:-1]
        list[1::2]=lst[1:]-1e-10 #ensure monotonic
    else:
        list=np.zeros(2*lst.size)
        list[::2]=lst-width/2
        list[1::2]=lst+width/2
        #make monotonic
        w=1+2*np.where(list[1:-1:2]>=list[2::2])[0]
        list[w] = list[1+w]-1e-10
    return(list)

def partition1d(pxlist,values,list,thres=1):
    '''
    newpx,bind,noperbin,w,binlist = partition1d(pixellist,values,list,thres=1)
       select and assign bins for averaging over.
    inputs:
       pixellist is list of pixels to consider.
       values are the values assigned to each pixel (like qs,phis,intensities).
       list are the break points (see mkpartlist).
       thres is minimum number of pixels in a bin. Default is 1).
    outputs:
       newpx is the selected subset of pixellist.
       bind is the bin assignment for each selected pixel (has all values 0 to n-1).
       noperbin is number of pixels in each bin.
       w is index into pixellist to generate newpx (ie newpx=pixellist[w]).
          can be used to select from other lists (ie if qs used then phis[w] are phis).
       binlist are bins left from list. (IE list[2*binlist[i]]<=x<list[1+2*binlist[i]]).
    '''
    bind=np.digitize(values,list)
    w=np.where(bind%2)[0]
    newpx=pxlist[w]
    bind=bind[w]//2
    noperbin=np.bincount(bind)
    binlist=np.arange(len(noperbin))
    #check have enough pixels in each bin
    w0=noperbin>=thres
    if(np.sum(w0)!=len(noperbin)):
        n=len(noperbin)-np.sum(w0)
        msg="Warning: {:d} bin{:s} less than threshold.".format(n,"s are" if n>1 else " is")
        print(msg,"bins:",binlist[~w0])
    #select pixels left
    ww=np.where(w0[bind])[0]
    bind=bind[ww]
    newpx=newpx[ww]
    w=w[ww]
    #make bin list
    noperbin=noperbin[w0]
    binlist[w0]=np.arange(len(noperbin))
    #make bins fill 0 to len(noperbin)-1
    bind=binlist[bind]
    binlist=np.where(w0)[0]
    return(newpx,bind,noperbin,w,binlist)

def partition2d(pxlist,values0,list0,values1,list1,thres=1):
    '''
    newpx,bind,bind0,bind1,noperbin,w,binlist0,binlist1 =
            partition2d(pixellist,values0,list0,values1,list1,thres=1)
       select and assign bins for averaging over.
    inputs:
       pixellist is list of pixels to consider.
       values[0,1] are the values assigned to each pixel (like qs,phis,intensities).
       list[0,1] are the break points (see mkpartlist).
       thres is minimum number of pixels in a bin. Default is 1).
    outputs:
       newpx is the selected subset of pixellist.
       bind is the bin assignment for each selected pixel (has all values
           0 to len(noperbin)-1).
       bind[0,1] are the assignments into list[0,1] of each pixel.
       noperbin is number of pixels in each bin.
       w is index into pixellist to generate newpx (ie newpx=pixellist[w]).
          can be used to select from other lists (ie if qs used then phis[w] are phis).
       binlist[0,1] are resutling bins from list[0,1]
    '''
    bind0=np.digitize(values0,list0)
    bind1=np.digitize(values1,list1)
    w=np.where((bind0%2)*(bind1%2))[0]
    newpx=pxlist[w]
    bind0=bind0[w]//2
    bind1=bind1[w]//2
    mn0=np.min(bind0)
    mx1=np.max(bind1)
    mn1=np.min(bind1)
    bind=(bind0-mn0)*(mx1-mn1+1)+(bind1-mn1)#give unique bind numbers
    noperbin=np.bincount(bind)
    binlist=np.arange(len(noperbin))
    #check enough pixels in each bin
    w0=noperbin>=thres
    if(np.sum(w0)!=len(noperbin)):
        n=len(noperbin)-np.sum(w0)
        msg="Warning: {:d} bin{:s} less than threshold.".format(n,"s are" if n>1 else " is")
        print(msg,"bins:",binlist[~w0])
    #select pixels left
    ww=np.where(w0[bind])[0]
    bind0=bind0[ww]
    bind1=bind1[ww]
    bind=bind[ww]
    newpx=newpx[ww]
    w=w[ww]
    # make bin lists
    noperbin=noperbin[w0]
    btmp=binlist[w0]
    binlist0=btmp//(mx1-mn1+1)+mn0
    binlist1=btmp%(mx1-mn1+1)+mn1
    #make bins fill from 0 to len(noperbin)-1
    binlist[w0]=np.arange(np.sum(w0))
    bind=binlist[bind]
    return(newpx,bind,bind0,bind1,noperbin,w,binlist,binlist0,binlist1)
