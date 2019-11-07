#!/usr/bin/env python
# coding: utf-8
import sys
import numpy as np
import matplotlib.pyplot as plt
#get_ipython().run_line_magic('matplotlib', 'inline')

#pfam = 'PF00504'
#s = np.loadtxt('pfam_2_40k.txt',dtype='str')
#pfam_list = s[:,0]
pfam_list = ['PF00200']

#=========================================================================================
def di_top(d,top):
    # find value of top biggest
    d1 = d.copy()
    np.fill_diagonal(d1, 0)
    #print(d1)
    
    a = d1.reshape((-1,))
    #print(a)    
    a = np.sort(a)[::-1] # descreasing sort
    #print(a)

    top_value = a[top]
    #print(top_value)
       
    # fill the top largest to be 1, other 0
    top_pos = d1 > top_value
    #print(top_pos)
    d1[top_pos] = 1.
    d1[~top_pos] = 0.
    #print(d1)
    
    xy = np.argwhere(d1==1)  
    return xy

#=========================================================================================
top_list = [40,60,80,100]
    
for pfam in pfam_list:
    ct = np.loadtxt('../pfam_50_80pos/%s_ct.txt'%pfam)
    di = np.loadtxt('%s/di.dat'%pfam)

    nx,ny = 4,5
    nfig = nx*ny
    fig, ax = plt.subplots(ny,nx,figsize=(nx*3.,ny*2.8))

    for j,cutoff in enumerate([3,4,5,6,7]):
        ct_top = np.argwhere(ct < cutoff)

        for i,top in enumerate(top_list):
            xy_di = di_top(di,top)

            #ax[j,i].plot(ct_top[:,0],ct_top[:,1],'ko',markersize=2)
            #ax[j,i].plot(xy_di[:,0],xy_di[:,1],'ro',markersize=4)
            ax[j,i].plot(ct_top[:,0],ct_top[:,1],'co',markersize=5,mfc='none',label='contact map')
            ax[j,i].plot(xy_di[:,0],xy_di[:,1],'r*',markersize=6,label='direct information')

    for i in range(4):
        ax[0,i].set_title('top: %i'%top_list[i])

    plt.tight_layout(h_pad=1, w_pad=1.5)
    
    plt.savefig('%s/ct_di_top.pdf'%pfam,format='pdf', dpi=50)
    plt.close()
       
    #-------------------------------------------------------------------   
    # map
    plt.figure(figsize=(8,3.2))

    plt.subplot2grid((1,2),(0,0))
    plt.title('contact map')
    plt.imshow(ct,cmap='rainbow_r',origin='lower')
    plt.xlabel('i')
    plt.ylabel('j')
    plt.clim(0,10)
    plt.colorbar(fraction=0.045, pad=0.05)

    plt.subplot2grid((1,2),(0,1))
    plt.title('direct info')
    plt.imshow(di,cmap='rainbow',origin='lower')
    plt.xlabel('i')
    plt.ylabel('j')
    plt.clim(0,0.01)
    plt.colorbar(fraction=0.045, pad=0.03)

    plt.tight_layout(h_pad=1, w_pad=1.5)
    plt.savefig('%s/ct_di.pdf'%pfam,format='pdf', dpi=100)  
    
    plt.close()
