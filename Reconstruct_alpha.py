# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 15:30:25 2020

@author: Can Kizilkale
"""

import os
import time
import copy
import numpy as np
import pandas as pd
from datetime import datetime
from argparse import ArgumentParser
import pygraphviz as pyg
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'






def GPR_NA_3(M_input):             # takes noisy matrix where NA are marked as 3
    
    NA_position=np.argwhere(M_input==3)
    for j in NA_position:
        M_input[j[0],j[1]]=0
    
    M_copy=M_input.astype(bool)
    approximate_order=greedy_row_rec(M_copy)[1]
    M_out=greedyPtreeNADyn(M_copy,sum(approximate_order),M_copy,0,21,0)[2]
    
    return M_out








###############################################################################
#############    Helper Functions  ############################################


def greedy_row_rec(M_input):     # reconstructing wrt. rows instead of columns
    sorted_rows=np.argsort(sum(M_input.T))
    ISet=[]
    M_copy=M_input.astype(bool)
    M_remaining=M_copy.copy()
    for i in range(M_copy.shape[0]):
        ISet.append(sorted_rows[i])
 
    S_rem=[]
    bret=[]                      #Location of detected false positives
        
    while len(ISet)>1:
        
        pivot_index=ISet[0]                   # index of the pivot vector 
        Sremaining=ISet.copy()
        pivot_vector=M_copy.T[:,pivot_index]  # vector used for pivoting the current iteration
        cum_vector=np.copy(pivot_vector)      # holds the union vector
        while_cont=1
        
        while while_cont==1:                  # main loop
            while_cont=0
            
            for j in Sremaining:
                cap_j=cum_vector*M_copy.T[:,j] # intersection of the pivot and the jth column
                
                
                if np.sum(cap_j)>0:            # continue as long as there is a column having non-empty intersection
                    M_copy.T[:,j]=M_copy.T[:,j]+pivot_vector
                    while_cont=1
                    Sremaining.remove(j)
               

        M_copy.T[:,pivot_index]=cum_vector
        ISet.remove(pivot_index)
            
    bret=np.argwhere(M_copy.astype(int)<M_input.astype(int))
    return [bret,M_copy]

def greedyPtreeNADyn(M_input,approx_order,M_grnd,overlap_coeff,hist_coeff,use_cht): # Modified Greedy algorithm for NA and false positive values.
    # Matrix M_input has to be BOOLEAN for the code to work right
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    #overlap_coeff=0
    #hist_coeff= 15                       # AROUND 21 for NRM, 22 for CHT seems optimum.
    M_copy=M_input.copy()
    ISet1=np.argsort(sum(M_copy))         # Set of indices of columns sorted wrt. the number of ones they contain
    if use_cht==1:
        approx_order=sum(M_grnd)
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

    bret=[]    #Location of detected false negatives
    pret=[]    #Location of detected false positives
    
    def n_of_cuts(S_i,c_vec,M_current):
        n_cut=0
        max_index=S_i[-1]
        for j in S_i:                                    # clean up the false positives wrt. the established pivot
            capj=c_vec*M_current[:,j]                   # intersection of union with column j
            difj=M_current[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                n_cut=n_cut+np.sum(M_copy[:,j])-np.sum(capj)
                if approx_order[max_index]<approx_order[j]:
                    max_index=j
                #M_copy[:,j]=capj
            else:
                n_cut=n_cut+np.sum(M_copy[:,j])-np.sum(difj)
                M_copy[:,j]=difj
        return n_cut,max_index
    
    def clean_fp(S_i,c_vec,M_current):
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=c_vec*M_current[:,j]                   # intersection of union with column j
            difj=M_current[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                M_current[:,j]=capj
            else:
                M_current[:,j]=difj
    S_rcn=[]
    while len(ISet)>1:

        pivot_index=ISet[-1]               # index of the pivot vector 
        Sremaining=ISet.copy()             # set of indices that are not included in the union
        pivot_vector=M_copy[:,pivot_index] # vector used for pivoting the current iteration
        cum_vector=np.copy(pivot_vector)   # holds the union vector
        while_cont=1
        cum_hist=np.zeros(M_input.shape[0]) # holds the histogram for the union
       
        
        while while_cont==1:               # Continue uniting vectors until no candidate remains
            while_cont=0
            for j in Sremaining:
                cap_i=pivot_vector*M_copy[:,j]
                min_vec_size=np.min([np.sum(M_copy[:,j]),np.sum(pivot_vector)])
                cap_size=np.sum(cap_i)
#               print(np.sum(M_copy[:,j]),np.sum(pivot_vector),np.sum(cap_i),cap_size, min_vec_size,np.floor(cap_size/min_vec_size))
                
                if  cap_size/min_vec_size > overlap_coeff: # we check if the columns have a meaningful overlap
                    cum_hist=cum_hist+M_copy[:,j].astype(int)
                    while_cont=1                          # found an overlapping vector so we keep on going
                    Sremaining.remove(j)                  # united vector is removed from the search set
                    #cnumT=np.floor(cum_hist.max()/hist_coeff)
                    #pivot_vector=cum_hist>cnumT
                
        cnumT=np.floor(cum_hist.max()/hist_coeff)                 # the elements that repeated few times are considered to be false positives
        cum_vector=cum_hist>cnumT
        [n_1_0,n_max_index]=n_of_cuts(ISet,cum_vector,M_copy)
        #print(n_max_index,pivot_index)
        
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=cum_vector*M_copy[:,j]                   # intersection of union with column j
            difj=M_copy[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                M_copy[:,j]=capj
            else:
                M_copy[:,j]=difj
        
        #if np.sum(cum_vector)==0:
        #            print("error pivot column is zero !!!")
        pivot_index=n_max_index                           # guessed pivot
        M_copy[:,pivot_index]=cum_vector                  # correcting the false negatives in the pivot
        ISet.remove(pivot_index)                          # removing the pivot from the search space 
        S_rcn.append(pivot_index)

    bret=np.argwhere(M_copy.astype(int)>M_input.astype(int))
    pret=np.argwhere(np.argwhere(M_copy.astype(int)<M_input.astype(int)))
    return [bret,pret,M_copy]


#
#file_name='simNo_1-s_30-m_4430-h_1-minVAF_0.02-ISAV_0-n_94-fp_0.001-fn_0.05-na_0-d_0-l_10000.SC'
#file_ground='simNo_1-s_30-m_4430-h_1-minVAF_0.02-ISAV_0-n_94-fp_0.001-fn_0.05-na_0-d_0-l_10000.SC.before_FP_FN_NA'
#m_g_sim=ReadFfile(file_ground)
#m_n_sim=ReadFileNA(file_name)
#m_rec=GPR_NA_3(m_n_sim)
#a=compareAD(m_g_sim,m_rec)
#compareDF(m_g_sim,m_rec)

