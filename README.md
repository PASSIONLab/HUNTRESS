# HUNTRESS

HUNTRESS is a fast algorithm for reconstructing phylogenetic trees of tumor evolution.

Huntress is written is python. 

## Package Dependencies
Tested on:

- Python   (ver 3.6) 
- Numpy    (version 1.20.3)
- Pandas   (version 1.3.4)

> We encourage using Anaconda for package management(not necessary). 

## Quick Start
If python and dependencies are already installed HUNTRESS can quickly be used as follows. 

Download "HUNTRESS.py" from [here]( https://github.com/PASSIONLab/HUNTRESS/blob/master/HUNTRESS.py).

Open Terminal/powershell in the folder containing "HUNTRESS.py" and type:

`python HUNTRESS.py --i "Input_filename" --o "Output_filename" `

- 'Input_filename' the path to the file that contains the noisy matrix.
- 'Output_filename' the path to the output file. Reconstructed matrix is written in "Output_filename.CFMATRIX". 


## Installation

Installation is very easy. Just follow these steps.   

Make sure that python version>=3.6 is installed.

Download the repository. From the terminal type the following

`$ git clone https://github.com/PASSIONLab/HUNTRESS  `

Go to the repo directory.

`$ cd HUNTRESS`

Inside the repo type the following to install dependencies. 

`$ pip install -r requirements.txt`

Finally test the script:

`$ python Test.py`

At the end of the run you should see "Test PASSED ...". And you are good to go.


## Usage

Open a terminal/powershell and go to the directory containing "HUNTRESS.py" and type.

`python HUNTRESS.py --i "Input_filename" --o "Output_filename `

- 'Input_filename' the path to the file that contains the noisy matrix.
- 'Output_filename' the path to the output file. Reconstructed matrix is written in "Output_filename.CFMATRIX". 

> If using Anaconda use Anaconda's own powershell/terminal.  

### List of arguments/tuneable parameters and their default values

- '--i' Path to the input file

- '--o' Path to the output file 

- '--t' defines the number of threads to be used for tuning in parallel. Default is number of threads of the computers it is run. 

- '--algorithmchoice' defines the version of the algorithm to be used.
           = "FN" for matrices that only have false negatives
           = "FPNA" for matrices that have false positives , false negatives and NA (entries that could not be read) entries. These entries must be given as 3 in the input matrix

- '--fn_fpratio' is the ratio of the weights of 0->1 switches over 1->0 switches that is used by the algorithm to tune the parameters.
 Default=51            

- '--fp_coeff' false positive probability coefficient (used for postprocessing).
 Default: 0.0001

- '--fn_coeff' false negative probability coefficient (used for postprocessing).
 Default: 0.1   

## Input Output Format:
The format of the input matrix is as follows.

Single-cell input is assumed to be represented in the form of ternary, __tab-delimited__, matrix with rows corresponding to single-cells and columns corresponding to mutations. We assume that this file contains headers and that matrix is ternary matrix with 0 denoting the absence and 1 denoting the presence of mutation in a given cell, whereas 3 represents the lack of information about presence/absence of mutation in a given cell (i.e. missing entry). __In order to simplify parsing of the matrix, we also assume that upper left corner equals to string `cellID/mutID`__.

Below is an example of single-cell data matrix. Note that mutation and cell names are arbitrary strings not containing tabs or spaces, however they must be unique.
```
cellID/mutID  mut0  mut1  mut2  mut3  mut4  mut5  mut6  mut7
cell0         0     0     3     0     0     0     0     0
cell1         0     3     1     0     0     0     1     1
cell2         0     0     1     0     0     0     1     1
cell3         1     1     0     0     0     0     0     0
cell4         0     0     1     0     3     0     0     0
cell5         1     0     0     0     0     0     0     0
cell6         0     0     1     0     0     0     1     1
cell7         0     0     1     0     0     0     0     0
cell8         3     0     0     0     3     0     3     1
cell9         0     1     0     0     0     0     0     0
```

The output is of the same format. Since it represents the reconstructed matrix its elements are only 1's and 0's. 

### EXAMPLE and Testing 

As stated in installation a quick test can be performed by running the following while inside the folder containing repo files. 
`python Test.py`
If you see "Test PASSED ..." at the end of the execution it is done.

For manual testing do the following.

Download Demo.Sc from https://github.com/PASSIONLab/HUNTRESS/blob/master/Demo.SC 

`python HUNTRESS.py --i "Demo.SC" --o "Demo_reconstructed" `

The reconstructed matrix will be written in "Demo_reconstructed.CFMatrix"

For testing purposes check if the matrices in files "Demo_out.CFMatrix" and "Demo_reconstructed.CFMatrix" are identical.


### Additional functions and outputs
The code includes additional testing functions CompareAD(Ground Truth Matrix,Reconstructed Matrix) computes ancestor descendent score of the reconstructed matrix,CompareDF(Ground Truth Matrix,Reconstructed Matrix) which computes ancestor descendent score of the reconstructed matrix. The code also generates a log file as "Output_filename.LOG". 



