Heuristic_r2

Use Reconstruct as the main function.
``` Reconstruct(input_file,output_file,Algchoice,S_or_Alpha=0,auto_tune=0,overlapp_coeff=0,hist_coeff=80):``` 

The path and name of the the noisy matrix is given in input_file
The reconstructed error free matrix is written in output_file

Algchoice defines the version of the algorithm to be used.
* = "FN" for matrices that only have false negatives
* = "FPFN" for matrices that have both false positives and false negatives
* = "NAFPFN" for matrices that have false positives , false negatives and NA (entries that could not be read) entries  that are marked as 3

S_or_Alpha 
* = 0 uses standard version of the algorithms(Default), you should use this.
* = 1 uses an experimental version that seems to give improved results for extreme NA cases such as very fat matrices, highly non-uniform errors.
            
auto_tune = automatically adjusts hist_coeff to minimize 1-0 switches from the noisy matrix to reconstructed one
* = 0 OFF (Default)
* = 1 ON

overlap_coeff,hist_coeff = reconstruction parameters if one wants to adjust them. USe Autotune instead

EXAMPLE: ```Reconstruct("simNo_1-s_10-m_100-h_1-minVAF_0.03-ISAV_0-n_100-fp_0.0001-fn_0.1-na_0-d_0-l_1000000.SC", "myoutput", "NAFPFN",0,1)``` will reconstruct the given matrix for a matrix with NAs, uses auto tune(alpha mode off) and write it to file myoutput.sc

