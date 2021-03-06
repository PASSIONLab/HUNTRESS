# HUNTRESS

HUNTRESS is a fast heuristic for reconstructing phylogenetic trees of tumor evolution.

The code requires Anaconda installed. 

Call this function like:

python huntress.py "Input_matrix_filename" "Output_filename" --nofcpus 8 --algorithmchoice "FPNA" --fn_fpratio 100 --fp_coeff 0.0001 --fn_coeff 0.1

Input Matrix: The path and name of the the noisy matrix is given here. 
Output Matrix: The reconstructed error free matrix is written in Output_filename with the extension ".CFMatrix"

The format of the input matrix is as follows.

Single-cell input is assumed to be represented in the form of ternary, __tab-delimited__, matrix with rows corresponding to single-cells and columns corresponding to mutations. We assume that this file contains headers and that matrix is ternary matrix with 0 denoting the absence and 1 denoting the presence of mutation in a given cell, whereas 3 represents the lack of information about presence/absence of mutation in a given cell (i.e. missing entry). __In order to simplify parsing of the matrix, we also assume that upper left corner equals to string `cellID/mutID`__.

Below is an example of single-cell data matrix. Note that mutation and cell names are arbitrary strings not containing tabs or spaces, however they must be unique.
```
cellID/mutID  mut0  mut1  mut2  mut3  mut4  mut5  mut6  mut7
cell0         0     0     3     0     0     0     0     0
cell1         0     3     1     0     0     0     1     1
cell2         0     0     1     0     0     0     1     1
cell3         1     1     0     0     0     0     0     0
cell4         0     0     1     0     0     0     0     0
cell5         1     0     0     0     0     0     0     0
cell6         0     0     1     0     0     0     1     1
cell7         0     0     1     0     0     0     0     0
cell8         3     0     0     0     3     0     3     1
cell9         0     1     0     0     0     0     0     0
```


The optional inputs are as follows:

--nofcpus defines the number of cpus to be used for tuning in parallel. Default is 7

--algorithmchoice defines the version of the algorithm to be used.
           = "FN" for matrices that only have false negatives
           = "FPNA" for matrices that have false positives , false negatives and NA (entries that could not be read) entries. These entries must be given as 3 in the input matrix

--fn_fpratio is the ratio of the weights of 0->1 switches over 1->0 switches that is used by the algorithm to tune the parameters.
 Default=100            

--fp_coeff false positive probability coefficient used for postprocessing.
 Default: 0.0001

--fn_coeff false negative probability coefficient used for postprocessing.
 Default: 0.1 

This draw_tree function can be used to visualize the Reconstructed matrix.


# Drawing tree

For drawing the tree please visit [https://github.com/algo-cancer/PhISCS-BnB](https://github.com/algo-cancer/PhISCS-BnB)
