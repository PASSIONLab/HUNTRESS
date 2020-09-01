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

# ---------------------READ ME ---------------------------
# Use Reconstruct as the main function.
# Reconstruct(input_file,output_file,Algchoice,S_or_Alpha=0,auto_tune=0,overlapp_coeff=0,hist_coeff=80):
# The path and name of the the noisy matrix is given in input_file
# The reconstructed error free matrix is written in output_file
# Algchoice defines the version of the algorithm to be used.
# Algchoice = "FN" for matrices that only have false negatives
#           = "FPFN" for matrices that have both false positives and false negatives
#           = "NAFPFN" for matrices that have false positives , false negatives and NA (entries that could not be read) entries  that are marked as 3
# S_or_Alpha = 0 uses standard version of the algorithms(Default)
#            = 1 uses an experimental version that seems to give improved results for extreme NA cases such as very fat matrices, highly non-uniform errors.
# auto_tune = automatically adjusts hist_coeff to minimize 1-0 switches from the noisy matrix to reconstructed one
# auto_tune = 0 OFF (Default)
#           = 1 ON
# overlap_coeff,hist_coeff = reconstruction parameters if one wants to adjust them. USe Autotune instead
#
# EXAMPLE: Reconstruct("simNo_1-s_10-m_100-h_1-minVAF_0.03-ISAV_0-n_100-fp_0.0001-fn_0.1-na_0-d_0-l_1000000.SC", "myoutput", "NAFPFN",0,1) will reconstruct the given matrix for a matrix with NAs, uses auto tune(alpha mode off) and write it to file myoutput.sc
#------------------- END README ---------------------------------------


def Reconstruct(input_file,output_file,Algchoice,S_or_Alpha=0,auto_tune=0,overlapp_coeff=0,hist_coeff=80):
    
    
    
                
    if Algchoice == "FN":
        matrix_input=ReadFile(input_file)
        matrix_recons=greedyPtreeNew(matrix_input.astype(bool))[1]
        WriteTfile(output_file,matrix_recons,input_file)
        return
    
    if Algchoice == "FPFN" and auto_tune==0:
        matrix_input=ReadFile(input_file)
        matrix_recons=greedyPtreeFPnew(matrix_input.astype(bool),overlapp_coeff,hist_coeff)[2]
        WriteTfile(output_file,matrix_recons,input_file)
        return
    
    if Algchoice == "NAFPFN" and auto_tune==0:
        matrix_input=ReadFileNA(input_file)
        apprx_ordr=Approx_csize(input_file)
        if S_or_Alpha == 1:
            apprx_ordr=sum(greedy_row_rec(matrix_input.astype(bool))[1])
            
        matrix_recons=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,overlapp_coeff,hist_coeff)[2]
        WriteTfile(output_file,matrix_recons,input_file)
        return
    
    if Algchoice == "FPFN" and auto_tune==1:
        print("auto tuning coeffcicient for minimum 1-0 switches")
        matrix_input=ReadFile(input_file)
        matrix_recons=greedyPtreeFPnew(matrix_input.astype(bool),0,1)[2]
        for i in range(2,100):
            matrix_rec_i=greedyPtreeFPnew(matrix_input.astype(bool),0,i)[2]
            n_10=np.sum(matrix_rec_i<matrix_input)
            print("H coefficient = ",i,"10 switches = ", n_10 , "best :",np.sum(matrix_recons<matrix_input))
            if n_10>2*np.sum(matrix_recons<matrix_input):
                break
            if n_10<np.sum(matrix_recons<matrix_input):
                matrix_recons=matrix_rec_i
            
            
        #matrix_recons=greedyPtreeFPnew(matrix_input.astype(bool),overlapp_coeff,hist_coeff)[2]
        WriteTfile(output_file,matrix_recons,input_file)
        return
    
    if Algchoice == "NAFPFN" and auto_tune==1:
        
        matrix_input=ReadFileNA(input_file)
        apprx_ordr=Approx_csize(input_file)
        if S_or_Alpha == 1:
            apprx_ordr=sum(greedy_row_rec(matrix_input.astype(bool))[1])
            print("Alpha mode ON!")
        print("auto tuning coeffcicient for minimum 1-0 switches for NA")    
        matrix_recons=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,0,1)[2]
        for i in range(2,100):
            matrix_rec_i=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,0,i)[2]
            n_10=np.sum(matrix_rec_i<matrix_input)
            print("H coefficient = ",i,"10 switches = ", n_10 , "best :",np.sum(matrix_recons<matrix_input))
            if n_10>2*np.sum(matrix_recons<matrix_input):
                break
            if n_10<np.sum(matrix_recons<matrix_input):
                matrix_recons=matrix_rec_i
        WriteTfile(output_file,matrix_recons,input_file)
        return    
    


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

def greedyPtreeFPnew(M_input,overlap_coeff=0.5,hist_coeff=20): # Modified Greedy algorithm for false positives.
    # Matrix M_input has to be BOOLEAN for the code to work right
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    #overlap_coeff=0.5
    #hist_coeff=20
    
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

def greedyPtreeNA(M_input,approx_order,oc,hc): # Modified Greedy algorithm for NA and false positive values.
    # Matrix M_input has to be BOOLEAN for the code to work right
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    overlap_coeff=oc  # 0.1
    hist_coeff=hc     # 25
    
    M_copy=M_input.copy()
#    ISet1=np.argsort(sum((M_copy.T.dot(M_copy)>0).T))    # Set of indices of columns sorted wrt. number of overlapping vectors
#    ISet1=np.argsort(sum(M_copy))         # Set of indices of columns sorted wrt. the number of ones they contain
    #ISet1=np.argsort(sum(M_NA))
    ISet1=np.argsort(approx_order)
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

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
        
#        if(np.sum(pivot_vector) > 10):
#            overlap_coeff=0.2
#        else:
#            overlap_coeff=0
            
        
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

def greedyPtreeNADyn_old(M_input,approx_order,M_grnd,hist_coeff): # Modified Greedy algorithm for NA and false positive values.
    # Matrix M_input has to be BOOLEAN for the code to work right
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    overlap_coeff=0
    #hist_coeff=80 
    M_copy=M_input.copy()
#    ISet1=np.argsort(sum((M_copy.T.dot(M_copy)>0).T))    # Set of indices of columns sorted wrt. number of overlapping vectors
#    ISet1=np.argsort(sum(M_copy))         # Set of indices of columns sorted wrt. the number of ones they contain
    #ISet1=np.argsort(sum(M_grnd))
    ISet1=np.argsort(approx_order)
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])

    bret=[]    #Location of detected false negatives
    pret=[]    #Location of detected false positives
    
    def n_of_cuts(S_i,c_vec,M_current):
        n_cut=0
        for j in S_i:                                    # clean up the false positives wrt. the established pivot
            capj=c_vec*M_current[:,j]                   # intersection of union with column j
            difj=M_current[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                n_cut=n_cut+np.sum(M_copy[:,j])-np.sum(capj)
                #M_copy[:,j]=capj
            else:
                n_cut=n_cut+np.sum(M_copy[:,j])-np.sum(difj)
                M_copy[:,j]=difj
        return n_cut
    
    def clean_fp(S_i,c_vec,M_current):
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=c_vec*M_current[:,j]                   # intersection of union with column j
            difj=M_current[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                M_current[:,j]=capj
            else:
                M_current[:,j]=difj
    while len(ISet)>1:
#        pivot_index=ISet[np.argmax(sum((M_copy[:,ISet].T.dot(M_copy[:,ISet])>0).T))]
        pivot_index=ISet[-1]               # index of the pivot vector 
        Sremaining=ISet.copy()             # set of indices that are not included in the union
        pivot_vector=M_copy[:,pivot_index] # vector used for pivoting the current iteration
        cum_vector=np.copy(pivot_vector)   # holds the union vector
        while_cont=1
        cum_hist=np.zeros(M_input.shape[0]) # holds the histogram for the union
        
#        if(np.sum(pivot_vector) > 10):
#            overlap_coeff=0.2
#        else:
#            overlap_coeff=0
            
        
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
        n_1_0=n_of_cuts(ISet,cum_vector,M_copy)
        
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=cum_vector*M_copy[:,j]                   # intersection of union with column j
            difj=M_copy[:,j]!=capj                        # difference of column j from the union 
            if np.sum(capj) > np.sum(difj):
                M_copy[:,j]=capj
                
#                if np.sum(capj) > 4*np.sum(pivot_vector):
#                    M_copy[:,j]=cum_vector
                    
            else:
                M_copy[:,j]=difj
        
        #print(n_1_0)
#        if n_1_0 > 50:
#            cum_vector=np.ones(M_copy.shape[0])
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

def ReadFile(filename): # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values
    return M

def ReadFileNA(filename): # reads the file and fills the NA with 0's.
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values
    
    NA_position=np.argwhere(M==3)
    #print("number of NA : ",len(NA_position))
    for j in NA_position:
        #print(M[j[0],j[1]])
        M[j[0],j[1]]=0
        
    #print(sum(sum(M)))    
    return M.astype(bool)

def ReadFileNA2(filename): # reads the file and fills the NA with 0's.
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values.astype(float)
    n_ones=len(np.argwhere(M==1))
    n_zeros=len(np.argwhere(M==0))
    o_rate=n_ones/(n_ones+n_zeros)
    t=0
    NA_position=np.argwhere(M==3)
    #print(len(NA_position))
    for j in NA_position:
        #print(M[j[0],j[1]])
        n1=len(np.argwhere(M[:,j[1]]==1))
        n0=len(np.argwhere(M[:,j[1]]==0))
        #c=len(np.argwhere(M[:,j[1]]==1))/(np.sum(M[:,j[1]]==1)+np.sum(M[:,j[1]]==0))
        if n1+n0 == 0:
            c=0
        else:
            c=n1/(n1+n0)
        if (n1+n0<t):
            c=o_rate
        M[j[0],j[1]]=c
        
    #print(sum(sum(M)))    
    return M

def Approx_csize(filename): # reads the file and fills the NA with 0's.
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values.astype(float)
    approx_lengths=np.zeros(M.shape[1])
    f_n=0.2
    f_p=0.005
    tresh=8
    
    n_ones=len(np.argwhere(M==1))
    n_zeros=len(np.argwhere(M==0))
    gen_or=n_ones/(n_ones+n_zeros)
    
    for j in range(M.shape[1]):
        if len(np.argwhere(M[:,j]<3)) < tresh:
            approx_lengths[j]=(1+f_n)*(gen_or*(len(np.argwhere(M[:,j]==3)))+len(np.argwhere(M[:,j]==1)))
        else:
            ones_rate=len(np.argwhere(M[:,j]==1))/(len(np.argwhere(M[:,j]==1))+len(np.argwhere(M[:,j]==0)))
            approx_lengths[j]=(1+f_n)*(ones_rate*(len(np.argwhere(M[:,j]==3)))+len(np.argwhere(M[:,j]==1)))
            
            
    #print("Expected n of ones: ",sum(approx_lengths))    
    return approx_lengths

def WriteTfile(filename,matrix,filename2): # writes matrix output as an integer matrix
    df_input = pd.read_csv(filename2, sep='\t', index_col=0)
    matrix_output=matrix.astype(int)
    df_output = pd.DataFrame(matrix_output)
    df_output.columns = df_input.columns
    df_output.index = df_input.index
    df_output.index.name = "cellIDxmutID"
    df_output.to_csv(filename, sep="\t")

##########################################################################3
############# TESTING FUNCTIONS ##########################################
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
    print('Number of AD pairs = ',n_adpairs,"errors : ",len(error_pairs), "AD score = ", 1 - len(error_pairs)/n_adpairs)
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
            
def ReadFasis(filename): # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values
    return M        

def compute_fnfp(M_n,M_r):
    n_01=0
    n_10=0
    n_1=0
    n_0=0
    for x in np.argwhere(M_n<3):
        if M_r[x[0],x[1]]==0:
            n_0=n_0 + 1
        if M_r[x[0],x[1]]==1:
            n_1=n_1 + 1    
        if M_n[x[0],x[1]] > M_r[x[0],x[1]]:
            
            n_10=n_10 + 1
            
            
        if M_n[x[0],x[1]] < M_r[x[0],x[1]]:
            n_01=n_01 + 1
            n_1=n_1 + 1
    print("computed fn :",n_01/n_1," fp : ",n_10/n_0)
    return [n_01/n_1,n_10/n_0]
        
# =====================================================================
        
#######################################################################
######################## ROW based ptree reconstructor ################
       
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
            
        
            
