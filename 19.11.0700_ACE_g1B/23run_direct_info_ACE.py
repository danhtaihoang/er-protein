import os
import numpy as np

pfam_list = ['PF00091','PF00092','PF00125','PF00174','PF00190','PF00207','PF00239','PF00339',\
             'PF00388','PF00394','PF00479','PF00640','PF00743','PF00763','PF00925']

#n = len(pfam_list)
#--------------------------------------------------------------
# run ACE
for pfam in pfam_list:
    print(pfam)

    os.system('python 3direct_info_ACE.py %s &'%pfam)

