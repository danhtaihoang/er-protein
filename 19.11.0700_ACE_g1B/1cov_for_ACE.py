#!/usr/bin/env python
# coding: utf-8
import sys
import numpy as np
from sklearn.preprocessing import OneHotEncoder

pfam_id = sys.argv[1]
#pfam_id = 'PF00200'

#s0 = np.loadtxt('../pfam_50_80pos/%s_s0.txt'%(pfam_id)).astype(int)
s0 = np.loadtxt('../pfam_2_100k/%s_s0.txt'%(pfam_id)).astype(int)
print(s0.shape)

n_var = s0.shape[1]
mx = np.array([len(np.unique(s0[:,i])) for i in range(n_var)])
mx_cumsum = np.insert(mx.cumsum(),0,0)
i1i2 = np.stack([mx_cumsum[:-1],mx_cumsum[1:]]).T 

mx_sum = mx.sum()

onehot_encoder = OneHotEncoder(sparse=False,categories='auto')
#onehot_encoder = OneHotEncoder(sparse=False)

s = onehot_encoder.fit_transform(s0)
print(s.shape)

#=========================================================================================
def hamming(seq1,seq2):
    """
    Return the Hamming distance between two sequences (of any kind).
    """

    d = np.sum(np.array(seq1)!=np.array(seq2))

    return d

#=========================================================================================
def seqreweight(msa,threshold=1.0):
    """
    Return weighting for a set of sequences based on similarity.
    """

    # Sequence reweighting (similarity)
    
    Beff=1.0
    weight=np.array([1.0 for i in range(len(msa))])

    thresh=int((1.0 - threshold) * len(msa[0]))
    for i in range(len(msa)):
        for j in range(i+1,len(msa)):
            if hamming(msa[i],msa[j])<thresh:
                weight[i]+=1.0
                weight[j]+=1.0
    weight=1.0/weight
    Beff=weight.sum(0)
        
    return Beff, weight
#=========================================================================================
msa = s

theta = 0.2
threshold = 1 - theta

Beff,weight = seqreweight(msa, threshold=threshold)

# Remove sites with no variation
p1      = np.sum(weight * msa.T, axis=1) / Beff
nonsing = (p1>0) * (p1<1)

#if removeSingular:
#    msa = msa[:,nonsing]
#    N = len(msa[0])

# Compute correlations
p12 = np.einsum('i,ij,ik->jk', weight, msa, msa) / Beff

p_ia = p1      # for convenience
p_iajb = p12   # for convenience

#=========================================================================================
# print frequency and correlation
with open('%s/cov.p'%pfam_id,'a') as f:
    # p_ia:
    for i0 in range(n_var):
        i1,i2 = i1i2[i0,0],i1i2[i0,1]

        # (exclude the last state)
        f.write(" ".join([str(p_ia[ia]) for ia in range(i1,i2-1)])+"\n")

    # p_iajb:
    for i0 in range(n_var-1):
        i1,i2 = i1i2[i0,0],i1i2[i0,1]

        for j0 in range(i0+1,n_var):
            j1,j2 = i1i2[j0,0],i1i2[j0,1]

            f.write(" ".join([str(p_iajb[ia,jb]) for ia in range(i1,i2-1) \
                                     for jb in range(j1,j2-1)])+"\n") 

#=========================================================================================
# find ia_tab, jb_tab
nrow = int(n_var*(n_var-1)/2.) 
ncol = int(20*20)

ia_tab = np.zeros((nrow,ncol)).astype(int)
jb_tab = np.zeros((nrow,ncol)).astype(int)

# row
irow = 0
for i0 in range(n_var-1):
    i1,i2 = i1i2[i0,0],i1i2[i0,1]

    for j0 in range(i0+1,n_var):
        j1,j2 = i1i2[j0,0],i1i2[j0,1]

        # column
        icol = 0
        for ia in range(i1,i2-1):
            for jb in range(j1,j2-1):
                
                ia_tab[irow,icol] = ia
                jb_tab[irow,icol] = jb
                
                icol += 1
        irow += 1
                
np.savetxt('%s/ia_tab.dat'%pfam_id,ia_tab,fmt='%i')                
np.savetxt('%s/jb_tab.dat'%pfam_id,jb_tab,fmt='%i')

