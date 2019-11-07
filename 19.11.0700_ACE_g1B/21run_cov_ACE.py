import os
import numpy as np

#s = np.loadtxt('pfam_2_40k.txt',dtype='str')
#pfam_list = s[:,0]

pfam_list = ['PF00091','PF00092','PF00125','PF00174','PF00190','PF00207','PF00239','PF00339',\
             'PF00388','PF00394','PF00479','PF00640','PF00743','PF00763','PF00925']

n = len(pfam_list)
#--------------------------------------------------------------
# create pfam folder
for i in range(n):
    os.system('rm -r %s'%(pfam_list[i]))
    os.system('mkdir %s'%(pfam_list[i]))

#--------------------------------------------------------------
# run covariance
for pfam in pfam_list:
    print(pfam)
    os.system('python 1cov_for_ACE.py %s &'%pfam)
