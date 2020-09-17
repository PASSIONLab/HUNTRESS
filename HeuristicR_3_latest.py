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
import networkx as nx
from datetime import datetime
from argparse import ArgumentParser
import pygraphviz as pyg

os.environ["PATH"] += os.pathsep + 'C:/Program Files (x86)/Graphviz2.38/bin/'
np.set_printoptions(threshold=np.inf)

# ---------------------READ ME ---------------------------
# Use Reconstruct as the main function.
# Reconstruct(input_file,output_file,Algchoice,S_or_Alpha=0,auto_tune=0,overlapp_coeff=0,hist_coeff=80):
# The path and name of the the noisy matrix is given in input_file
# The reconstructed error free matrix is written in output_file. Output file name must end with ".CFMatrix" for proper operation
# Algchoice defines the version of the algorithm to be used.
# Algchoice = "FN" for matrices that only have false negatives
#           = "FPFN" for matrices that have both false positives and false negatives
#           = "NAFPFN" for matrices that have false positives , false negatives and NA (entries that could not be read) entries  that are marked as 3
# approx_method defines which approximation will be used for number wise column ordering
# approx_method = 0 columns are ordered wrt. how many other columns they overlap
#               = 1 columns are ordered wrt. the number of ones they contain.
#            
# auto_tune = automatically adjusts hist_coeff to minimize 1-0 switches from the noisy matrix to reconstructed one
# auto_tune = 0 OFF 
#           = 1 ON  (Default)
# overlap_coeff,hist_coeff = reconstruction parameters if one wants to adjust them. USe Autotune instead
#
# postprocessing minimizes the total number of 1->0,0->1while keeping the topology of the reconstructed graph. The files will have extension ".processed"
# postprocessing = 0 postprocessing OFF 
#                = 1 postprocessing ON  (default)    it still finds and writes the unprocessed matrix so keeping this always on is a good idea
# EXAMPLE: Reconstruct("simNo_1-s_10-m_100-h_1-minVAF_0.03-ISAV_0-n_100-fp_0.0001-fn_0.1-na_0-d_0-l_1000000.SC", "myoutput.CFMatrix", "NAFPFN",approx_method=0,auto_tune=1,postprocessing=1) will reconstruct the given matrix for a matrix with NAs,using approximation wrt overlaps, output to myoutput.CFMATrix
# and write the postprocessed version to "myoutput.processed.CFMatrix"
#------------------- END README ---------------------------------------


def Reconstruct(input_file,output_file,Algchoice="NAFPFN",approx_method=0,auto_tune=1,overlapp_coeff=0,hist_coeff=80,postprocessing=1):
    
    
    matrix_input=ReadFileNA(input_file)
    if approx_method == 0:
        print("Ordering by number of overlaps...")
        M=matrix_input.astype(int)
        apprx_ordr=sum((M.T.dot(M))>0)
    if approx_method==1:
        print("Ordering by number of ones...")
        apprx_ordr=sum(matrix_input)
        
    tune_range=85   
#    apprx_ordr=Approx_csize(input_file)
#    M=matrix_input.astype(int)
#    apprx_ordr=sum((M.T.dot(M))>0)
#    M=ReadFfile("yardena.scistree.CFMatrix")
#    apprx_ordr=sum(M_rowdene)
    
#    if S_or_Alpha == 1:
#        apprx_ordr=sum(greedy_row_rec(matrix_input.astype(bool))[1])
#        print("Alpha mode ON!")
#    if S_or_Alpha == 2:
##        matrix_input=greedy_row_rec2(matrix_input.astype(bool))
#        apprx_ordr=sum(greedy_row_rec2(matrix_input.astype(bool)))
#        print("Super Alpha mode ON!")
#    if S_or_Alpha == 3:
#        M=matrix_input.astype(int)
#        apprx_ordr=sum((M.T.dot(M))>0)
#        print("Approximating order by number of overlaps...")
        
#    matrix_cht=ReadFfile("Yardena.CFMatrix")
#    apprx_ordr=sum(matrix_cht)
    
                
    if Algchoice == "FN":

        matrix_recons=greedyPtreeNew(matrix_input.astype(bool))[1]
        WriteTfile(output_file,matrix_recons,input_file)
        return
    
    if Algchoice == "FPFN" and auto_tune==0:

        matrix_recons=greedyPtreeFPnew(matrix_input.astype(bool),overlapp_coeff,hist_coeff)[2]
        WriteTfile(output_file,matrix_recons,input_file)
        return
    
    if Algchoice == "NAFPFN" and auto_tune==0:
           
        matrix_recons=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,overlapp_coeff,hist_coeff)[2]
        WriteTfile(output_file,matrix_recons,input_file)
        print(" difference ",np.sum(matrix_recons!=matrix_input))
#        return
    
    if Algchoice == "FPFN" and auto_tune==1:
        print("auto tuning coeffcicient for minimum 1-0 switches")
        matrix_recons=greedyPtreeFPnew(matrix_input.astype(bool),0,1)[2]
        h_current=2
        for i in range(2,tune_range):
            matrix_rec_i=greedyPtreeFPnew(matrix_input.astype(bool),0,i)[2]
            n_10=np.sum(matrix_rec_i<matrix_input)
            print("H coefficient = ",i,"10 switches = ", n_10 , "best :",np.sum(matrix_recons<matrix_input))
#            if n_10>10*np.sum(matrix_recons<matrix_input):
#                break
            if n_10<np.sum(matrix_recons<matrix_input):
                matrix_recons=matrix_rec_i
                WriteTfile(output_file,matrix_recons,input_file)
                h_current=i
            
        #matrix_recons=greedyPtreeFPnew(matrix_input.astype(bool),overlapp_coeff,hist_coeff)[2]
        print(" Optimum h is ",h_current)
        WriteTfile(output_file,matrix_recons,input_file)
        return
    
    if Algchoice == "NAFPFN" and auto_tune==1:
        h_current=2
        print("auto tuning coeffcicient for minimum 1-0 switches for NA")    
        matrix_recons=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,0,1)[2]
        for i in range(3,tune_range):
            matrix_rec_i=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,0,i)[2]
            n_10=np.sum(matrix_rec_i<matrix_input)
            n_01=np.sum(matrix_rec_i>matrix_input)
        
            print("Opt h =",h_current,"H coefficient = ",i," 01 switches : ",n_01,"10 switches = ", n_10 , "best :",np.sum(matrix_recons!=matrix_input))
#            if n_10>10*np.sum(matrix_recons<matrix_input):
#                break
            if n_10+n_01<np.sum(matrix_recons!=matrix_input):
                matrix_recons=matrix_rec_i
                WriteTfile(output_file,matrix_recons,input_file)
                h_current=i
            if n_10+n_01>2*np.sum(matrix_recons!=matrix_input):
                break
            
        WriteTfile(output_file,matrix_recons,input_file)
        print(" difference ",np.sum(matrix_recons!=matrix_input)," at h=",h_current)
    if postprocessing==1:
        postprocess(input_file,output_file)
        
        
    

def deletemutations(M_in,M_out): # Finds the rows which are too far from the input, just replicates them to their closest neighbour.
    x=M_in.shape
    M_return=M_out.copy()
    treshold_cut=0.5
    dt=np.divide(sum(M_in),sum(M_out))<=treshold_cut
    for i in range(x[1]):
        if dt[i]==1:
            M_return[:,i]=np.zeros(x[0])
    
    return M_return

def precombination(M_in):
    x=M_in.shape
    trsh=0.22
    M_return=M_in.copy()
    for i in range(x[0]):
        for j in range(i+1,x[0]):
            rij=np.sum(M_in[i,:]!=M_in[j,:])/np.sum(M_in[i,:]+M_in[j,:])
            print(rij)
            if rij<trsh:
                cupi=M_in[i,:]+M_in[j,:]
                M_return[i,:]=cupi
                M_return[j,:]=cupi
                print("combined ",i,j)
    return M_return

def findclosestinter(i,M_input):
    vec_int=np.zeros(M_input.shape[1]).astype(bool)
    for j in range(M_input.shape[0]):
        if i!=j and np.sum(M_input[j,:]>M[i,:])>0 and np.sum(M_input[j,:]*M[i,:])>0:
            if np.sum(vec_int):
                vec_int=M_input[j,:]
            else:
                vec_int=vec_int*M_input[j,:]
    return vec_int

            
    
            
        
        
        
                            


# Both algorithms take inut matrices of BOOL type.
def greedyPtreeNew(M_input): #very greedy algorithm that constructs a ptree matrix from M_inputs by adding 1's
    # M_input has to be of type BOOLEAN !
    # Returns 
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
#    M_apprx=greedy_row_rec(M_input)[1]    
    M_copy=M_input.copy()

    ISet1=np.argsort(sum(M_copy))
#    ISet1=np.argsort(sum(M_apprx))
#    ISet1=np.argsort(sum(M_rowdene))
    ISet=[]
    for i in range(M_copy.shape[1]):
        ISet.append(ISet1[i])
    
    bret=[]    #Location of detected false negatives
    print(M_copy.shape,len(ISet))    
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
#            difj=M_copy[:,j]!=capj                        # difference of column j from the union 
            difj=M_copy[:,j]>capj                        # difference of column j from the union
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

        pivot_est=pivot_index
        
        ncap=np.sum(cum_vector*M_copy[:,pivot_est])
#        ndel=np.sum(M_copy[:,j]!=capj)
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=cum_vector*M_copy[:,j]                   # intersection of union with column j
#            difj=M_copy[:,j]!=capj                        # difference of column j from the union
            difj=M_copy[:,j]>capj                        # difference of column j from the union
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

def greedyPtreeNA2(M_input,approx_order,oc,hc): # Modified Greedy algorithm for NA and false positive values.
    # Matrix M_input has to be BOOLEAN for the code to work right
    # Returns:
    # M_copy(reconstructed matrix),
    # bret (list of positions of assumed false negatives)
    # pret (list of approximate false positives)
    overlap_coeff=oc  # 0.1
    hist_coeff=hc     # 25
    def numberofdifference(IS,c_vec):
        n_of_changes=0
        for j in IS:                                    # clean up the false positives wrt. the established pivot
            capj=c_vec*M_copy[:,j]                   # intersection of union with column j
#            difj=M_copy[:,j]!=capj                        # difference of column j from the union
            difj=M_copy[:,j]>capj                        # difference of column j from the union
            
            if np.sum(capj) > np.sum(difj):
                n_of_changes=n_of_changes+np.sum(difj)
            
            else:
                n_of_changes=n_of_changes+np.sum(capj)
#            print(n_of_changes)
        return n_of_changes
                
    def find_hc(IS,c_hist,p_vec):
        n_current=10000
        h_current=1
        for i in range(100):
            cnT=np.floor(c_hist.max()/i)                 # the elements that repeated few times are considered to be false positives
            c_vector=c_hist>cnT
            n_temp=numberofdifference(IS,c_vector) + 0.1*np.sum(p_vec!=c_vector)
            if n_temp<n_current:
                n_current=n_temp
                h_current=i
        print("found ", n_current, " h_coefficient ",h_current)
        return h_current
            
            
        
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
    it_count=0
    while len(ISet)>1:
#        pivot_index=ISet[np.argmax(sum((M_copy[:,ISet].T.dot(M_copy[:,ISet])>0).T))]
        pivot_index=ISet[-1]               # index of the pivot vector 
        Sremaining=ISet.copy()             # set of indices that are not included in the union
        pivot_vector=M_copy[:,pivot_index] # vector used for pivoting the current iteration
        cum_vector=np.copy(pivot_vector)   # holds the union vector
        while_cont=1
        cum_hist=np.zeros(M_input.shape[0]) # holds the histogram for the union
        it_count=it_count+1
        
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
        
        hist_coeff=find_hc(ISet,cum_hist,pivot_vector)                 # dynamically deterines h_coeff
        cnumT=np.floor(cum_hist.max()/hist_coeff)                 # the elements that repeated few times are considered to be false positives
        cum_vector=cum_hist>cnumT

        pivot_est=pivot_index
        
        ncap=np.sum(cum_vector*M_copy[:,pivot_est])
#        ndel=np.sum(M_copy[:,j]!=capj)
        for j in ISet:                                    # clean up the false positives wrt. the established pivot
            capj=cum_vector*M_copy[:,j]                   # intersection of union with column j
#            difj=M_copy[:,j]!=capj                        # difference of column j from the union
            difj=M_copy[:,j]>capj                        # difference of column j from the union
            if np.sum(capj) > np.sum(difj):
                M_copy[:,j]=capj

            else:
                M_copy[:,j]=difj
        
        
        M_copy[:,pivot_index]=cum_vector                  # correcting the false negatives in the pivot
        ISet.remove(pivot_index) 
        print(np.sum(M_input!=M_copy),"Step : ",it_count)                         # removing the pivot from the search space       
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

def greedy_row_rec2(M_input):     # reconstructing wrt. rows instead of columns
    sorted_rows=np.argsort(sum(M_input.T))
    ISet=[]
    M_copy=M_input.astype(bool)
    M_remaining=M_copy.copy()
    for i in range(M_copy.shape[0]):
        ISet.append(sorted_rows[i])
 
    S_rem=[]
    bret=[]                      #Location of detected false positives
        
    while len(ISet)>1:
        pivot_index=ISet[np.argmax(sum(M_copy[ISet,:].T))]
#        pivot_index=ISet[-1]                   # index of the pivot vector 
        Sremaining=ISet.copy()
        pivot_vector=M_copy[pivot_index,:]  # vector used for pivoting the current iteration
#        cap_hist=pivot_vector.astype(int)      # holds the union vector
        capj=pivot_vector*M_copy[ISet[0],:]
        cap_hist=capj
        for j in Sremaining:
            capjT=pivot_vector*M_copy[j,:]
            if j!=pivot_index:
                cap_hist=cap_hist+capjT
            if np.sum(capj) < np.sum(capjT) and j!=pivot_index:
                capj=capjT
                jt=j
#                cap_hist=cap_hist+capj.astype(int)
#        hist_T=2
#        cap_vec=cap_hist>hist_T
#        cap_vec=capj
        cap_vec=cap_hist
        print("Intersect : ",np.sum(cap_vec), " pivot length ",np.sum(pivot_vector),np.sum(cap_hist))
        if pivot_index==jt:
            print("Error index conflict")
#        print("max intersection :",np.sum(cap_vec), " pivot vector ", np.sum(pivot_vector))
        
        for j in Sremaining:
#            capij=M_copy[pivot_index,:]*M_copy[j,:]
            if np.sum(cap_vec>pivot_vector) >0:
                print("Error pivot does not contain the intersection")
#            print(np.sum(cap_vec>pivot_vector))
            capij=cap_vec*M_copy[j,:]
            difj=M_copy[j,:]>pivot_vector
#            print(np.sum(difj),np.sum(capij),np.sum(M_copy[j,:]))
#            print(np.sum(capij),np.sum(M_copy[pivot_index,:]*M_copy[j,:]))
            
            if j != pivot_index:
#                print(np.sum(capij),np.sum(M_copy[pivot_index,:]*M_copy[j,:]))
                print(np.sum(M_copy[j,:]),np.sum(capij),np.sum(difj))
                M_copy[j,:]=difj + capij
                print(np.sum(M_copy[j,:]))            
            
        print(np.sum(M_copy!=M_input))    
        ISet.remove(pivot_index)
            
    
    return M_copy
        
            
# Cenk's group code for printing a phylogeny graph from a Matrix.
#----------------------------------------------------------------    
def draw_tree(filename):
    add_cells = True

    from collections import Counter
    import pygraphviz as pyg
    import networkx as nx
    from networkx.drawing.nx_agraph import graphviz_layout, to_agraph

    def contains(col1, col2):
        for i in range(len(col1)):
            if not col1[i] >= col2[i]:
                return False
        return True

    df = pd.read_csv(filename, sep="\t", index_col=0)
    splitter_mut = "\n"
    matrix = df.values
    names_mut = list(df.columns)

    i = 0
    while i < matrix.shape[1]:
        j = i + 1
        while j < matrix.shape[1]:
            if np.array_equal(matrix[:, i], matrix[:, j]):
                matrix = np.delete(matrix, j, 1)
                x = names_mut.pop(j)
                names_mut[i] += splitter_mut + x
                j -= 1
            j += 1
        i += 1

    rows = matrix.shape[0]
    cols = matrix.shape[1]
    dimensions = np.sum(matrix, axis=0)
    indices = np.argsort(dimensions)
    dimensions = np.sort(dimensions)
    names_mut = [names_mut[indices[i]] for i in range(cols)]

    G = nx.DiGraph(dpi=300)
    G.add_node(cols)
    G.add_node(cols - 1)
    G.add_edge(cols, cols - 1, label=names_mut[cols - 1])
    node_mud = {}
    node_mud[names_mut[cols - 1]] = cols - 1

    i = cols - 2
    while i >= 0:
        if dimensions[i] == 0:
            break
        attached = False
        for j in range(i + 1, cols):
            if contains(matrix[:, indices[j]], matrix[:, indices[i]]):
                G.add_node(i)
                G.add_edge(node_mud[names_mut[j]], i, label=names_mut[i])
                node_mud[names_mut[i]] = i
                attached = True
                break
        if not attached:
            G.add_node(i)
            G.add_edge(cols, i, label=names_mut[i])
            node_mud[names_mut[i]] = i
        i -= 1

    clusters = {}
    for node in G:
        if node == cols:
            # G._node[node]['label'] = '<<b>germ<br/>cells</b>>'
            G._node[node]["fontname"] = "Helvetica"
            G._node[node]["width"] = 0.4
            G._node[node]["style"] = "filled"
            G._node[node]["penwidth"] = 3
            G._node[node]["fillcolor"] = "gray60"
            continue
        untilnow_mut = []
        sp = nx.shortest_path(G, cols, node)
        for i in range(len(sp) - 1):
            untilnow_mut += G.get_edge_data(sp[i], sp[i + 1])["label"].split(splitter_mut)
        untilnow_cell = df.loc[
            (df[untilnow_mut] == 1).all(axis=1)
            & (df[[x for x in df.columns if x not in untilnow_mut]] == 0).all(axis=1)
        ].index
        if len(untilnow_cell) > 0:
            clusters[node] = f'<b>{", ".join(untilnow_cell)}</b>'
        else:
            clusters[node] = "––"

        if add_cells:
            if "––" not in clusters[node]:
                G._node[node]["fillcolor"] = "#80C4DF"
            else:
                G._node[node]["fillcolor"] = "gray90"
            G._node[node]["label"] = clusters[node]
        else:
            G._node[node]["label"] = ""
            G._node[node]["shape"] = "circle"
        G._node[node]["fontname"] = "Helvetica"
        G._node[node]["width"] = 0.4
        G._node[node]["style"] = "filled"
        G._node[node]["penwidth"] = 2.5
    i = 1
    for k, v in clusters.items():
        if v == "––":
            clusters[k] = i * "––"
            i += 1

    for node in G:
        if node != cols:
            num = 0
            paths = nx.shortest_path(G, source=cols, target=node)
            for i in range(len(paths) - 1):
                x = paths[i]
                y = paths[i + 1]
                num += len(G[x][y]["label"].split(splitter_mut))
            G._node[node]["label"] = f"<[{node}]  " + G._node[node]["label"] + f"  ({num})>"
        else:
            G._node[node]["label"] = f"<[{node}]  germ cells>"

    data = G.edges.data("label")
    outputpath = filename[: -len(".CFMatrix")]
    for u, v, l in data:
        ll = l.split(splitter_mut)
        genes = [x.split(".")[0] for x in ll]
        a = Counter(genes)
        a = a.most_common()
        lll = list(set([x.split(".")[0] for x in ll]))
        G.add_edge(u, v, label=splitter_mut.join(lll))
        print(f"[{u}]->[{v}]: {' '.join(ll)}", file=open(f"{outputpath}.mutsAtEdges", "a"))
        G.add_edge(u, v, label=f" {len(ll)}")

    header = ""
    temp = df.columns[(df == 0).all(axis=0)]
    if len(temp) > 0:
        header += f"Became Germline: {len(temp)}<br/>"

    H = nx.relabel_nodes(G, clusters)
    html = """<{}>""".format(header)
    H.graph["graph"] = {
        "label": html,
        "labelloc": "t",
        "resolution": 300,
        "fontname": "Helvetica",
        "fontsize": 8,
    }
    H.graph["node"] = {"fontname": "Helvetica", "fontsize": 12}
    H.graph["edge"] = {"fontname": "Helvetica", "fontsize": 12, "penwidth": 2}

    mygraph = to_agraph(H)
    mygraph.layout(prog="dot")
    mygraph.draw(f"{outputpath}.png")


#-----------------------------------------------------------------------------
#-End of graph draw..
    
def estimated_graph(mutAtEdgesFile,noisySC):
    
    N = pd.read_table(noisySC, index_col=0)
    data = []
    with open(mutAtEdgesFile) as fin:
        for line in fin:
            line = line.strip()
            a = line.split('->')[0].replace('[','').replace(']','')
            b = line.split('->')[1].split(': ')[0].replace('[','').replace(']','')
            c = line.split('->')[1].split(': ')[1].split(' ')
            d = {}
            d['from'] = a
            d['to'] = b
            for x in c:
                d[x] = 1
            data.append(d)
    CAN = pd.DataFrame(data, columns=['from','to']+list(N.columns)).fillna(0).astype(int)
    M = CAN.drop(['from','to'], axis=1).values
    root = max(CAN['from'].values)
    G = nx.DiGraph()
    for x, y in CAN[['from','to']].itertuples(index=False):
        G.add_edge(x, y)
    def return_path(to):
        l = list(nx.shortest_path(G, root, to))
        o = []
        for z in [(l[i],l[i+1]) for i in range(0,len(l)-1,1)]:
            o.append(CAN[(CAN['from'] == z[0]) & (CAN['to'] == z[1])].index[0])
        return o
    Mp = []
    for x in CAN['to'].values:
        y = return_path(x)
        Mp.append(np.sum(M[y,:], axis=0))
    Mp = np.array(Mp)
    Mp[Mp > 0] = 1
    out = pd.DataFrame(Mp, columns=N.columns, index=CAN.to)
    
 
    Mt=CAN.to_numpy()[:,:2]
    root_node=np.max(Mt)
    print("root node num ", root_node)
    M_nodes=np.zeros((root_node+1,Mp.shape[1]))
    
    for i in range(root_node):
        node_ind=np.argwhere(Mt[:,1]==i)[0][0]    # initial node index
        M_nodes[i,:]=Mp[node_ind,:]
    #    c_vec=Mp[node_ind,:].copy()
    #    node_an=Mt[node_ind,0]
    #    print("Computing ",i,np.sum(c_vec))
    #    while node_an!=root_node:
    #        
    #        node_ind=np.argwhere(Mt[:,1]==node_an)[0][0]
    #        c_vec=c_vec+Mp[node_ind,:]
    #        print(node_an,np.sum(c_vec))
    #        node_an=Mt[node_ind,0]
    #    print(node_an,np.sum(c_vec))
    #    M_nodes[i,:]=c_vec 
    return M_nodes

def find_dist(node_piv,M_nodes):
    distances=np.zeros(M_nodes.shape[0])
    for i in range(M_nodes.shape[0]):
        distances[i]=np.sum(M_nodes[i,:]!=node_piv)
    return distances

def closest_matrix(M_input,M_nodes):
    M_out=M_input.copy()
    for i in range(M_input.shape[0]):
        pivot_v=M_input[i,:]
        distance_i=find_dist(pivot_v,M_nodes)
        min_index=np.argmin(distance_i)
        M_out[i,:]=M_nodes[min_index,:]
        print("M input :", sum(M_input[i,:]),sum(M_out[i,:]))
    return M_out
    
def postprocess(input_file,out_file):
    M_noisy=ReadFfile(input_file)
    
    edge_file=out_file[:-8]+"mutsAtEdges"
    if os.path.exists(edge_file):
        print("Edge file exists, deleting previous one")
        os.remove(edge_file)
    draw_tree(out_file)
    M_nds=estimated_graph(edge_file,input_file)
    M_postprocessed=closest_matrix(M_noisy,M_nds)
    processed_file=out_file[:-8] + "processed.CFMatrix"
    print("Writing to file ",processed_file)
    WriteTfile(processed_file,M_postprocessed,input_file)
    print("complete, difference  ",np.sum(M_postprocessed!=M_noisy))
    draw_tree(processed_file)
    
    
    