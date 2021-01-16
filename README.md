#Needs to be updated

HeuristicR_3_latest.py

Use Reconstruct as the main function.
```
Reconstruct(input_file,output_file,Algchoice,S_or_Alpha=0,auto_tune=0,overlapp_coeff=0,hist_coeff=80):
```
The path and name of the the noisy matrix is given in input_file

The reconstructed error free matrix is written in output_file. Output file name must end with ".CFMatrix" for proper operation

Algchoice defines the version of the algorithm to be used. Algchoice=
- "FN" for matrices that only have false negatives
- "FPFN" for matrices that have both false positives and false negatives
- "NAFPFN" for matrices that have false positives , false negatives and NA (entries that could not be read) entries  that are marked as 3

approx_method defines which approximation will be used for number wise column ordering
approx_method = 
- 0 columns are ordered wrt. how many other columns they overlap
- 1 columns are ordered wrt. the number of ones they contain.
            
auto_tune automatically adjusts hist_coeff to minimize 1-0 switches from the noisy matrix to reconstructed one. auto_tune = 
- 0 OFF 
- 1 ON  (Default)

overlap_coeff,hist_coeff arfe reconstruction parameters if one wants to adjust them. Use Autotune instead

postprocessing minimizes the total number of 1->0,0->1 while keeping the topology of the reconstructed graph. The files will have extension ".processed". postprocessing = 
- 0 postprocessing OFF 
- 1 postprocessing ON  (default)    it still finds and writes the unprocessed matrix so keeping this always on is a good idea
                
 EXAMPLE: 
 ```
 Reconstruct("simNo_1-s_10-m_100-h_1-minVAF_0.03-ISAV_0-n_100-fp_0.0001-fn_0.1-na_0-d_0-l_1000000.SC", "myoutput.CFMatrix", "NAFPFN",approx_method=0,auto_tune=1,postprocessing=1) 
 ```
 will reconstruct the given matrix for a matrix with NAs,using approximation wrt overlaps, output to myoutput.CFMATrix
 and write the postprocessed version to "myoutput.processed.CFMatrix"



