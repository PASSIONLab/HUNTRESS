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
# Reconstruct(input_file,output_file):    
# The path and name of the the noisy matrix is given in input_file
# The reconstructed error free matrix is written in output_file. The extension ".CFMatrix" will be added to output_file
# The optional inputs to the function Reconstruct are as follows:
# Algchoice defines the version of the algorithm to be used.
# Algchoice = "FN" for matrices that only have false negatives
#           = "FPNA" for matrices that have false positives , false negatives and NA (entries that could not be read) entries  that are marked as 3
#            
# auto_tune = automatically adjusts hist_coeff to minimize 1-0 switches from the noisy matrix to reconstructed one
# auto_tune = 0 OFF 
#           = 1 ON  (Default)
# overlap_coeff,hist_coeff = reconstruction parameters if one wants to adjust them. Use Autotune instead
#
# postprocessing minimizes the total number of 1->0,0->1 by modyfying rows. The files will have extension ".processed"
# postprocessing = 0 postprocessing OFF (default)
#                = 1 postprocessing ON   it still finds and writes the unprocessed matrix
# EXAMPLE: Reconstruct("simNo_1-s_10-m_100-h_1-minVAF_0.03-ISAV_0-n_100-fp_0.0001-fn_0.1-na_0-d_0-l_1000000.SC", "myoutput.CFMatrix", "NA",auto_tune=1,postprocessing=0) will reconstruct the given matrix for a matrix with NAs,using approximation wrt overlaps, output to myoutput_optH{the history coefficient used}.CFMATrix
# and write the postprocessed version to "myoutput.processed.CFMatrix"
#------------------- END README ---------------------------------------


def Reconstruct(input_file,output_file,Algchoice="FPNA",auto_tune=1,overlapp_coeff=0.1,hist_coeff=20,postprocessing=0,fnfp=1,fnc=1):
    Flog = open(output_file + ".LOG", "a+")
    matrix_input=ReadFileNA(input_file)
    matrix_input_raw=ReadFasis(input_file)
    matrix_NA_est=Estimated_Matrix(input_file)
    print(np.sum(matrix_input),np.sum(matrix_NA_est),file=Flog)

    tune_range=30   
    
    running_time = 0
    if Algchoice == "FN":
        s_time = time.time()
        matrix_recons=greedyPtreeNew(matrix_input.astype(bool))[1]
        e_time = time.time()
        running_time = e_time - s_time
        WriteTfile(output_file,matrix_recons,input_file)
 
    if Algchoice == "FPNA" and auto_tune==0:

        s_time = time.time()
        apprx_ordr=sum(matrix_NA_est)
        matrix_recons=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,overlapp_coeff,hist_coeff)[2]
        e_time = time.time()
        running_time = e_time - s_time
        output_file=output_file+"_optH_{}TEMP.CFMatrix".format(hist_coeff)
        WriteTfile(output_file,matrix_recons,input_file)
#        print(" difference ",np.sum(matrix_recons!=matrix_input))
    
    if Algchoice == "FPNA" and auto_tune==1:

        s_time = time.time()
        apprx_ordr=sum(matrix_NA_est)
        h_current=1
        print("auto tuning coeffcicient for minimum 1-0 switches for FP and/or NA",file=Flog)    
        matrix_recons=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,0,1)[2]
        matrix_rec_Temp=deleteNas(matrix_input_raw,matrix_recons)
        n_10=np.sum(matrix_rec_Temp<matrix_input)
        n_01=np.sum(matrix_rec_Temp>matrix_input)
        
#        distance=np.square(n_10)+n_01
        distance=fnfp*n_10+fnc*n_01
        oc_current=0
        for oc_j in range(0,6):
            overlapp_coeff=oc_j/10
            
            for i in range(h_current + 1,tune_range):
                
                matrix_rec_i=greedyPtreeNA(matrix_input.astype(bool),apprx_ordr,overlapp_coeff,i)[2]
                matrix_rec_Temp=deleteNas(matrix_input_raw,matrix_rec_i)
                n_10=np.sum(matrix_rec_Temp<matrix_input,dtype='int64')
                n_01=np.sum(matrix_rec_Temp>matrix_input,dtype='int64')
                print("Opt h = ",h_current,"Opt_O = ",oc_current,"H,OC coefficient = ",i,overlapp_coeff," 01 switches : ",n_01,"10 switches = ", n_10 , "best :",distance,file=Flog)
                print("Opt h = ",h_current,"Opt_O = ",oc_current,"H,OC coefficient = ",i,overlapp_coeff," 01 switches : ",n_01,"10 switches = ", n_10 , "best :",distance)
    #            distance_i=np.square(n_10)+n_01
                distance_i=fnfp*n_10+fnc*n_01
    #            if n_10<np.sum(matrix_recons<matrix_input):
                if distance_i<distance:
                    matrix_recons=matrix_rec_i.copy()
                    distance=distance_i
    #                WriteTfile(output_file,matrix_recons,input_file)
                    h_current=i
                    oc_current=overlapp_coeff
        e_time = time.time()
        running_time = e_time - s_time

        output_file=output_file+"_optH_TEMP.CFMatrix"    
        WriteTfile(output_file,matrix_recons,input_file)
        print(" difference ",np.sum(matrix_recons!=matrix_input)," at h=",h_current,file=Flog)
    
    Flog.close()
    postprocess_col(input_file,output_file)   #Column based post processing, seem to give the best result.
    if postprocessing==1:
        postprocess(input_file,output_file[:-13]+".CFMatrix")   #Column based post processing, seem to give the best result.
        
    return running_time
        
    
def deleteNas(M_in,M_out):
    M_o=M_out.copy()
    NA_position=np.argwhere(M_in==3)
    #print("number of NA : ",len(NA_position))
    for j in NA_position:
        #print(M[j[0],j[1]])
        M_o[j[0],j[1]]=0
        
    return M_o


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


def Estimated_Matrix(filename):            # Creates an estimate of the matrix such that each element is given the expectation wrt the column 1/0 frequencies.
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values.astype(float)
    
    for i in range(M.shape[1]):
        if np.sum(M[:,i]!=3)==0:
            one_ratio=0
        else:
            one_ratio=np.sum(M[:,i]==1)/np.sum(M[:,i]!=3)
        for j in range(M.shape[0]):
            if M[j,i]==3:
                M[j,i]=one_ratio


    return M



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
        

            
# Cenk's group code for printing a phylogeny graph from a Matrix.
#----------------------------------------------------------------    
def draw_tree(filename):
    add_cells = False
    change_edges_to_number = False
    combine = False

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
            if not add_cells:
                G._node[node]["label"] = "<<b>zygote</b>>"
            G._node[node]["fontname"] = "Helvetica"
            G._node[node]["style"] = "rounded"
            G._node[node]["shape"] = "box"
            G._node[node]["margin"] = 0.05
            G._node[node]["pad"] = 0
            G._node[node]["width"] = 0
            G._node[node]["height"] = 0
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
                G._node[node]["color"] = "#0000FF"
            else:
                G._node[node]["color"] = "gray70"
                G._node[node]["fontcolor"] = "gray70"
            G._node[node]["label"] = clusters[node]
        else:
            G._node[node]["label"] = ""
        G._node[node]["shape"] = "box"
        G._node[node]["fontname"] = "Helvetica"
        G._node[node]["style"] = "rounded"
        G._node[node]["margin"] = 0.05
        G._node[node]["pad"] = 0
        G._node[node]["width"] = 0
        G._node[node]["height"] = 0

    i = 1
    for k, v in clusters.items():
        if v == "––":
            clusters[k] = i * "––"
            i += 1

    outputpath = filename[: -len(".CFMatrix")]
    if add_cells:
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
                G._node[node]["label"] = f"<[{node}]  <b>zygote</b>  (0)>"

    if change_edges_to_number:
        data = G.edges.data("label")
        # with open(f"{outputpath}.mutsAtEdges", "w") as fout:
        for u, v, l in data:
            ll = l.split(splitter_mut)
            # fout.write(f"[{u}]->[{v}]: {' '.join(ll)}\n")
            if "––" in G._node[v]["label"]:
                G.add_edge(u, v, label=f"  {len(ll)}  ", color="gray70", fontcolor="gray70")
            else:
                G.add_edge(u, v, label=f"  {len(ll)}  ")


    if combine:
        H = G.copy()

        for _ in range(len(G.nodes)):
            d_in = H.in_degree(H)
            d_out = H.out_degree(H)
            for node in H.nodes():
                if d_out[node] == 1 and d_in[node] == 1:
                    parent = [x for x in H.predecessors(node)][0]
                    child = [x for x in H.successors(node)][0]
                    if d_out[parent] < 2 and d_in[parent] == 1:
                        new_node = f'{parent}{node}'
                        new_edge = f"{int(H[parent][node]['label']) + int(H[node][child]['label'])}"

                        H = nx.contracted_nodes(H, parent, node, self_loops=False)
                        mapping = {parent: new_node}
                        H = nx.relabel_nodes(H, mapping)
                        H[new_node][child]['label'] = new_edge
                        break

        d_in = H.in_degree(H)
        d_out = H.out_degree(H)
        nodes = []
        for node in H.nodes():
            if d_out[node] == 0:
                nodes.append(node)

        for node in nodes:
            parent = [x for x in H.predecessors(node)][0]
            if d_out[parent] == 1 and d_in[parent] == 1:
                grandparent = [x for x in H.predecessors(parent)][0]

                new_node = f'{parent}{node}'
                new_edge = f"{int(H[grandparent][parent]['label']) + int(H[parent][node]['label'])}"

                H = nx.contracted_nodes(H, parent, node, self_loops=False)
                mapping = {parent: new_node}
                H = nx.relabel_nodes(H, mapping)
                H[grandparent][new_node]['label'] = new_edge

        mygraph = nx.drawing.nx_agraph.to_agraph(H)
    else:
        mygraph = nx.drawing.nx_agraph.to_agraph(G)

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
 
    return Mp

def find_dist(node_piv,M_samples):
#    distances=np.zeros(M_nodes.shape[0])
#    for i in range(M_nodes.shape[0]):
##        distances[i]=np.sum(M_nodes[i,:]!=node_piv)
##        distances[i]=np.sum(M_nodes[i,:]<node_piv)
##        distances[i]=np.linalg.norm(M_nodes[i,:]-node_piv)
#        distances[i]=0
#        d_10=0
#        d_01=0
#        for j in range(M_nodes.shape[1]):
#            
#            if node_piv[j]!=3:
#                if node_piv[j]>M_nodes[i,j]:
#                    d_10=d_10+1
#                if node_piv[j]<M_nodes[i,j]:
#                    d_01=d_01+1
#        distances[i]= np.square(d_10) +d_01
##            else:
##                distances[i]=1-M_nodes[i,j]+distances[i]
    M_nodes=M_samples.copy()
    for j in range(M_nodes.shape[1]):
        if node_piv[j]==3:
            node_piv[j]=0
            M_nodes[:,j]=0*M_nodes[:,j]
    distances=np.zeros(M_nodes.shape[0])
    for i in range(M_nodes.shape[0]):
        d_10=np.sum(node_piv>M_nodes[i,:],dtype='int64')
        d_01=np.sum(node_piv<M_nodes[i,:],dtype='int64')
        distances[i]=np.square(d_10) + d_01            
                        
                
    return distances

def find_dist_col(node_piv,M_samples):
    M_nodes=M_samples.copy()
    for j in range(M_nodes.shape[0]):
        if node_piv[j]==3:
            node_piv[j]=0
            M_nodes[j,:]=0*M_nodes[j,:]
    distances=np.zeros(M_nodes.shape[1])
    for i in range(M_nodes.shape[1]):
        d_10=np.sum(node_piv>M_nodes[:,i],dtype='int64')
        d_01=np.sum(node_piv<M_nodes[:,i],dtype='int64')
#        distances[i]=np.square(d_10) + d_01
        distances[i]=d_10 + d_01
    return distances

def closest_matrix(M_input,M_nodes,M_rec):
    M_out=M_input.copy()
    for i in range(M_input.shape[0]):
        pivot_v=M_input[i,:]
        distance_i=find_dist(pivot_v,M_nodes)
        
        min_index=np.argmin(distance_i)
#        print("Current distance ",np.sum(M_input[i,:]>M_rec[i,:]),"minimum possible ",np.min(distance_i),min_index)
        M_out[i,:]=M_nodes[min_index,:]
#        if np.sum(M_input[i,:]>M_rec[i,:])>2*np.min(distance_i):
#            min_index=np.argmin(distance_i)
#            M_out[i,:]=M_nodes[min_index,:]
#        else:
#            M_out[i,:]=M_rec[i,:]
#        print("M input :", sum(M_input[i,:]),sum(M_out[i,:]))
    return M_out
 
def closest_matrix_col(M_input,M_nodes,M_rec):
    M_out=M_input.copy()
    for i in range(M_input.shape[1]):
        pivot_v=M_input[:,i]
        distance_i=find_dist_col(pivot_v,M_nodes)
        
        min_index=np.argmin(distance_i)
#        print("Old New difference ",np.sum(M_nodes[:,i]!=M_nodes[:,min_index]),i)
        M_out[:,i]=M_nodes[:,min_index]

    return M_out
def postprocess(input_file,out_file):
    M_noisy=ReadFasis(input_file)
    
    edge_file=out_file[:-9]+".mutsAtEdges"
    if os.path.exists(edge_file):
        print("Edge file exists, deleting previous one")
        os.remove(edge_file)
    draw_tree(out_file)
    M_nds=estimated_graph(edge_file,input_file)
    M_postprocessed=closest_matrix(M_noisy,M_nds,ReadFfile(out_file))
    processed_file=out_file[:-8] + "processed.CFMatrix"
#    print("Writing to file ",processed_file)
    WriteTfile(processed_file,M_postprocessed,input_file)

    draw_tree(processed_file)
    
def postprocess_col(input_file,out_file):
    M_noisy=ReadFasis(input_file)

    M_nds=ReadFfile(out_file)
    M_postprocessed=closest_matrix_col(M_noisy,M_nds,ReadFfile(out_file))
    processed_file=out_file[:-13] + ".CFMatrix"
#    print("Writing to file ",processed_file,file=Flog)
    WriteTfile(processed_file,M_postprocessed,input_file)

    draw_tree(processed_file)
    