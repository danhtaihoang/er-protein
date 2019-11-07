import os
import numpy as np

s = np.loadtxt('pfam_2_40k.txt',dtype='str')
pfam_list = s[:,0]

for pfam in pfam_list:
    print(pfam)
    os.system('python 1main_LR.py %s'%pfam)
