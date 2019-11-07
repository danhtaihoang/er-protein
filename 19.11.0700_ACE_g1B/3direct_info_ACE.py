#!/usr/bin/env python
# coding: utf-8

import sys
import numpy as np
#import pandas as pd
from direct_info import direct_info

np.random.seed(1)

#pfam_id = 'PF00200'
pfam_id = sys.argv[1]

#s0 = np.loadtxt('../pfam_50_80pos/%s_s0.txt'%(pfam_id)).astype(int)
s0 = np.loadtxt('../pfam_2_100k/%s_s0.txt'%(pfam_id)).astype(int)

n_var = s0.shape[1]
mx = np.array([len(np.unique(s0[:,i])) for i in range(n_var)])
mx_sum = mx.sum()

ia_tab = np.loadtxt('%s/ia_tab.dat'%pfam_id).astype(int)
jb_tab = np.loadtxt('%s/jb_tab.dat'%pfam_id).astype(int)

print(ia_tab.shape)
print(jb_tab.shape)

w = np.zeros((mx_sum,mx_sum))

with open('%s/w_pred.j'%pfam_id) as file:    
    for irow,line in enumerate(file):
        
        if irow >= n_var: # skip the local field part        
            irow1 = irow - n_var
            a = [float(x) for x in line.split('\t')]

            for icol in range(len(a)):
                ia = ia_tab[irow1,icol]
                jb = jb_tab[irow1,icol]

                w[ia,jb] = a[icol]
                
w = w + w.T                

# direct information
di = direct_info(s0,w)
np.savetxt('%s/di.dat'%pfam_id,di,fmt='% 3.8f')
