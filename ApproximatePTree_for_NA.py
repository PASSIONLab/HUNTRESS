# -*- coding: utf-8 -*-
"""
Created on Thu Feb 13 19:18:26 2020

@author: Can KIZILKALE
"""

#Approximate Error Correction for phylogenetic tree matrices.
#Can Kizilkale

import os
import time
import copy
import numpy as np
import pandas as pd
from datetime import datetime
from argparse import ArgumentParser
import pygraphviz as pyg
os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

# Both algorithms take inut matrices of BOOL type.
def greedyPtreeNew(M_input): #very greedy algorithm that constructs a ptree matrix from M_inputs by adding 1's
    # M_input has to be of type BOOLEAN !
    # Returns 
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
        
    M_copy=M_input.copy()
    M_copy=M_copy.astype(bool)
    ISet1=np.argsort(sum(M_copy))
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])
    
    bret=[]    #Location of detected false negatives
        
    while len(ISet)>1:
        #print("pivoting column", i,bret)
        pivot_index=ISet[-1]                # index of the pivot vector 
        Sremaining=ISet.copy()
        pivot_vector=M_copy[:,pivot_index]  # vector used for pivoting the current iteration
        cum_vector=np.copy(pivot_vector)    # holds the union vector
        while_cont=1
        
        while while_cont==1:                # main loop
            while_cont=0
            
            for j in Sremaining:
                cap_j=cum_vector*M_copy[:,j] # intersection of the pivot and the jth column
                
                if np.any(cap_j):            # continue as long as there is a column having non-empty intersection
                    cum_vector=cum_vector+M_copy[:,j]
                    while_cont=1
                    Sremaining.remove(j)
#                    print("pivot ", pivot_index, j, sum(cap_j),sum(cum_vector),sum(M_copy[:,j]))
#        print(len(Si))
#        print(sum(pivot_vector),sum(cum_vector))
        M_copy[:,pivot_index]=cum_vector
        ISet.remove(pivot_index)        
        bret=np.argwhere(M_copy.astype(int)>M_input.astype(int))
    return [bret,M_copy]

def greedyPtreeFPnew(M_input): # Modified Greedy algorithm for false positives.
    # Matrix M_input has to be BOOLEAN for the code to work right
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    overlap_coeff=0.5                     # Two columns having less than 50 percent overlap will not be connected  
    hist_coeff=20                         # The number elements less than 1/hist_coeff of the maximum element will be marked as false positives 
    
    M_copy=M_input.copy()
    M_copy=M_copy.astype(bool)
    ISet1=np.argsort(sum(M_copy))         # Set of indices of columns sorted wrt. the number of ones they contain
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

    bret=[]    #Location of detected false negatives
    pret=[]    #Location of detected false positives

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
        
        cnumT=np.floor(cum_hist.max()/hist_coeff)                 # the elements that repeated few times are considered to be false positives
        cum_vector=cum_hist>cnumT

        
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=cum_vector*M_copy[:,j]                   # intersection of union with column j
            difj=M_copy[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                M_copy[:,j]=capj
            else:
                M_copy[:,j]=difj
        
        
        M_copy[:,pivot_index]=cum_vector                  # correcting the false negatives in the pivot
        ISet.remove(pivot_index)                          # removing the pivot from the search space       
#        print(len(ISet),np.sum(M_copy[:,pivot_index]),np.sum(M_input[:,pivot_index]))
    
    bret=np.argwhere(M_copy.astype(int)>M_input.astype(int))
    pret=np.argwhere(np.argwhere(M_copy.astype(int)<M_input.astype(int)))
    return [bret,pret,M_copy]

def greedyPtreeNA(M_input): # Modified Greedy algorithm for NA and false positive values.
    # M_input is the noisy INTEGER matrix where NAs are represented as 3.
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    overlap_coeff=0.1
    hist_coeff=20

    [M_i,S_order]=Noisy_Read_Order(M_input)
    M_copy=M_i.copy()
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(S_order[i])

    bret=[]    #Location of detected false negatives
    pret=[]    #Location of detected false positives

    while len(ISet)>1:
#        pivot_index=ISet[np.argmax(sum((M_copy[:,ISet].T.dot(M_copy[:,ISet])>0).T))]
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
        
        cnumT=np.floor(cum_hist.max()/hist_coeff)                 # the elements that repeated few times are considered to be false positives
        cum_vector=cum_hist>cnumT

        
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=cum_vector*M_copy[:,j]                   # intersection of union with column j
            difj=M_copy[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                M_copy[:,j]=capj
            else:
                M_copy[:,j]=difj
        
        
        M_copy[:,pivot_index]=cum_vector                  # correcting the false negatives in the pivot
        ISet.remove(pivot_index)                          # removing the pivot from the search space       
#        print(len(ISet),np.sum(M_copy[:,pivot_index]),np.sum(M_input[:,pivot_index]))
    
    bret=np.argwhere(M_copy.astype(int)>M_input.astype(int))
    pret=np.argwhere(np.argwhere(M_copy.astype(int)<M_input.astype(int)))
    return [bret,pret,M_copy]

def greedyPtreeNAmod(M_input): # Modified Greedy algorithm for NA and false positive values, BETA.
    # M_input is the noisy INTEGER matrix where NAs are represented as 3.
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    overlap_coeff=0
    hist_coeff=7

    [M_i,S_order]=Noisy_Read_Order(M_input)
    M_copy=M_i.copy()
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(S_order[i])

    bret=[]    #Location of detected false negatives
    pret=[]    #Location of detected false positives
    
    while len(ISet)>1:
#        pivot_index=ISet[np.argmax(sum((M_copy[:,ISet].T.dot(M_copy[:,ISet])>0).T))]
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
                    pivot_vector=cum_hist>np.floor(cum_hist.max()/hist_coeff)  # this time we keep an updated pivot vector at each step.

                    Sremaining.remove(j)                  # united vector is removed from the search set
        
        cnumT=np.floor(cum_hist.max()/hist_coeff)                 # the elements that repeated few times are considered to be false positives
        cum_vector=cum_hist>cnumT

        
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=cum_vector*M_copy[:,j]                   # intersection of union with column j
            difj=M_copy[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                M_copy[:,j]=capj
            else:
                M_copy[:,j]=difj
        
        
        M_copy[:,pivot_index]=cum_vector                  # correcting the false negatives in the pivot
        ISet.remove(pivot_index)                          # removing the pivot from the search space       
#        print(len(ISet),np.sum(M_copy[:,pivot_index]),np.sum(M_input[:,pivot_index]))
    
    bret=np.argwhere(M_copy.astype(int)>M_input.astype(int))
    pret=np.argwhere(np.argwhere(M_copy.astype(int)<M_input.astype(int)))
    return [bret,pret,M_copy]


def Noisy_Read_Order(matrix_input): # Takes a noisy integer matrix where NA's are marked as 3
    # Returns a matrix where NA's are zeroed and an approximate ordering wrt the number of ones and NAs in the matrix
    M=matrix_input.astype(float)
    M_o=M.copy()
    NA_position=np.argwhere(M==3)
    print(len(NA_position))
    for j in NA_position:
        #print(M[j[0],j[1]])
        M_o[j[0],j[1]]=0
        c=np.sum(M[:,j[1]]==1)/(np.sum(M[:,j[1]]==1)+np.sum(M[:,j[1]]==0))
        if (np.sum(M[:,j[1]]==1)<0):
            c=0
        M[j[0],j[1]]=c
    print(sum(sum(M)))
    return [M_o.astype(bool),np.argsort(sum(M))]
    

def ReadFile(filename): # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep='\t', index_col=0)
    return df.values

def WriteTfile(filename,matrix,filename2): # writes matrix output as an integer matrix
    df_input = pd.read_csv(filename2, sep='\t', index_col=0)
    matrix_output=matrix.astype(int)
    df_output = pd.DataFrame(matrix_output)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"
    df_output.to_csv(filename, sep="\t")


def isPtree(matrix_in):   # brute force check if matrix_in is a pTree
    M=matrix_in.astype(bool)
    for j in range(M.shape[1]):
        for i in range(j,M.shape[1]):
            cap=M[:,i]*M[:,j]
            cap_size=np.sum(cap)
            Mi_size=np.sum(M[:,i])
            Mj_size=np.sum(M[:,j])
            if (cap_size != 0):
                if (cap_size != Mi_size):
                    if (cap_size != Mj_size):
                        return False
                
    print("Seems to be a PTree ...")
    return True

def compareAD(M1,M2): #M1 is the ground truth
    error_pairs=[]
    n_adpairs=0
    for i in range(M1.shape[1]):
#        print(i)
        for j in range(i,M1.shape[1]):
            cap1=M1[:,i]*M1[:,j]
            cap2=M2[:,i]*M2[:,j]
            if (np.sum(cap1)>0 and np.sum(M1[:,i]) != np.sum(M1[:,j])):
                n_adpairs=n_adpairs+1
                if (np.sum(cap2)==0):
                    error_pairs.append([i,j])
                else:
                    if (np.sum(M1[:,j])>np.sum(M1[:,i]) and np.sum(M2[:,j])<=np.sum(M2[:,i])):
                        error_pairs.append([i,j])
                    else:
                        if (np.sum(M1[:,i])>np.sum(M1[:,j]) and np.sum(M2[:,i])<=np.sum(M2[:,j])):
                            error_pairs.append([i,j])
                        #print(i,j,sum(M1[:,i]),sum(M1[:,j]),sum(M2[:,i]),sum(M2[:,j]))
    print('Number of AD pairs = ',n_adpairs,"errors : ",len(error_pairs), "error rate = ", 1 - len(error_pairs)/n_adpairs)
    return error_pairs                

def compareDF(M_orj,M_rec):
    error_pairs=[]
    d_pairs=0
    for i in range(M_orj.shape[1]):
        for j in range(i,M_orj.shape[1]):
            cap1=M_orj[:,i]*M_orj[:,j]
            cap2=M_rec[:,i]*M_rec[:,j]
            if np.sum(cap1)==0:
                d_pairs=d_pairs + 1
                if np.sum(cap2)>0:
                    error_pairs.append([i,j])
    print("Number of Diff pairs = ",d_pairs, "errors :",len(error_pairs), "score :", 1-len(error_pairs)/d_pairs)
    return 




###################################################################################################################
#                  TESTING 
###################################################################################################################    

file_name="simNo_10-s_15-m_1000-h_1-minVAF_0.04-ISAV_0-n_1000-fp_0.005-fn_0.1-na_0-d_0-l_1000000-seed_4.SC"
file_name_ground="simNo_10-s_15-m_1000-h_1-minVAF_0.04-ISAV_0-n_1000-fp_0.005-fn_0.1-na_0-d_0-l_1000000.SC.before_FP_FN_NA"

m_noisy=ReadFile(file_name)
m_ground=ReadFfile(file_name_ground)

print('Computing scores for :',file_name)
m_NRM=greedyPtreeNA(m_noisy)[2]
ADerrors=compareAD3(m_ground,m_NRM)
compareDF(m_ground,m_NRM)
#file_out_nrm=file_name + "NRM_20_1"
#WriteTfile(file_out_nrm,m_NRM,file_name)
#print("Now cheating :")
#m_CHT=greedyPtreeNA(m_noisy,m_ground)[2]
#ADerrors=compareAD3(m_ground,m_CHT)
#compareDF(m_ground,m_CHT)
#file_out_cht=file_name + "CHT_20_1"
#WriteTfile(file_out_cht,m_CHT,file_name)
