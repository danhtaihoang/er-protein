import os
import numpy as np

pfam_list = ['PF00091','PF00092','PF00125','PF00174','PF00190','PF00207','PF00239','PF00339',\
             'PF00388','PF00394','PF00479','PF00640','PF00743','PF00763','PF00925']

#n = len(pfam_list)
#--------------------------------------------------------------
# run ACE
for pfam in pfam_list:
    print(pfam)
    s0 = np.loadtxt('../pfam_2_100k/%s_s0.txt'%(pfam)).astype(int)
    n_seq = s0.shape[0]
    
    g2 = 1./n_seq

    #os.system('python 1cov_for_ACE.py %s &'%pfam)
    os.system('./bin/ace -d %s -i cov -o w_pred -g2 %f -b %i &'%(pfam,g2,n_seq))

