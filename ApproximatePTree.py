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
#import pygraphviz as pyg
#os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'

# Both algorithms take inut matrices of BOOL type.
def greedyPtreeNew(M_input): #very greedy algorithm that constructs a ptree matrix from M_inputs by adding 1's
    # M_input has to be of type BOOLEAN !
    # Returns 
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
        
    M_copy=M_input.copy()
    ISet1=np.argsort(sum(M_copy))
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])
    
    bret=[]    #Location of detected false negatives (not necessary)
        
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
    M_copy=M_input.copy()
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
                
                if  cap_size/min_vec_size > 0.5: # we check if the columns have a meaningful overlap
                    cum_hist=cum_hist+M_copy[:,j].astype(int)
                    while_cont=1                          # found an overlapping vector so we keep on going
                    Sremaining.remove(j)                  # united vector is removed from the search set
        
        cnumT=np.floor(cum_hist.max()/20)                 # the elements that repeated few times are considered to be false positives
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




def ReadFfile(filename): # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values.astype(bool)
    return M

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



#matrix_input=ReadFfile('simNo_1-s_10-m_100-n_100-fn_0.1-k_144.SC.noisy')
#
#start=time.time()
#val=greedyPtreeFPnew(matrix_input)
#end=time.time()
#print(end-start)

#start=time.time()
#val=greedyPtreeNew(matrix_input)
#end=time.time()
#print(end-start)
#matrix_out=val[1]
#print(matrix_out)
##matrix_oldout=ReadFfile("simNo_1-s_15-m_1000-h_1-minVAF_0.04-ISAV_0-n_100-fp_0.005-fn_0.20-na_0-d_0-l_1000000.SCFpApprox")
##matrix_out=correctitFp(matrix_input,val)
##print(sum(matrix_input)-sum(matrix_out))
#isPtree(matrix_out)
## =============================================================================
#fnames=os.listdir()
#fnames.pop(0)
#sancheck=[]
#
#
#start=time.time()
#for finputname in fnames:
#    #finputname=fnames[2]
#    matrix_input=ReadFfile(finputname)
##    print("Approximation started for :",finputname)
#    val=greedyPtreeFPnew(matrix_input)
#    #M= correctitFp (matrix_input,val)
#    M=val[2]
#    isitP=isPtree(M)
##    print('approximation completed',sum(sum(M)),sum(sum(matrix_input)))
#    #print(finputname)
#    
#    sancheck.append(isitP)    
#    
#    fapproximated=finputname+"Approx-m_500"
#    WriteTfile(fapproximated,M,finputname)
#
#end=time.time()
#print("Approximation ended, elapsed time :",end-start) 







