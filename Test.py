# Simple test for functionality of HUNTRESS
# Author: Can Kizilkale

import numpy as np
import pandas as pd
import os



def ReadFasis(filename):     # reads a matrix from file and returns it in BOOL type
    df = pd.read_csv(filename, sep='\t', index_col=0)
    M=df.values
    return M        

inputfile="Demo.SC"
outputfile="Demo_test"
reconstructedfile=outputfile + ".CFMatrix"
testfile="Demo_out.CFMatrix"
mycommand="python HUNTRESS.py "+ " --i "+ inputfile + " --o " + outputfile# used default values
os.system(mycommand)
M_r=ReadFasis(reconstructedfile)
M_test=ReadFasis(testfile)

if np.array_equal(M_r,M_test):
    print("Test PASSED ...")
else:
    print("Test Failed !!!")







# files=os.listdir()
# n_cpu=8
# for f in files:
#     Scorelog=open("Doub003.scores", "a+")
#     inputfile=f
#     outputfile=inputfile
#     groundfile=inputfile + ".before_FP_FN_NA"
#     if f.endswith(".SC") and f.startswith("simNo") and os.path.exists(groundfile):
#         Pfn=float(inputfile.split("fn_")[1].split("-")[0])
#         Pfp=0.03
#         reconstructedfile=inputfile+".CFMatrix"
#         #mycommand="python huntress.py "+ inputfile + " " + outputfile + " " + "--nofcpus " + str(n_cpu)  +" --fp_coeff "+str(Pfp) +" --fn_coeff "+str(Pfn) + " --fn_fpratio " + str(10)# used default values
#         mycommand="python huntress_cleaned.py "+ inputfile + " " + outputfile + " --fn_fpratio " + str(2)+" --nofcpus " + str(n_cpu)# used default values              
#         #+ " --fp_coeff " + str(0.00001) + " --fn_fpratio " + str(25)+ " --fn_fpratio " + str(25)
#         t_start=time.time()
#         #Reconstruct(args.inputfile,args.outputfile,Algchoice="FPNA",n_proc=mp.cpu_count(),fnfp=51,post_fn=fn_conorm,post_fp=fp_conorm)
#         os.system(mycommand)
#         t_end=time.time()
#         M_g=ReadFasis(groundfile)
#         if os.path.exists(reconstructedfile):
#             M_r=ReadFasis(reconstructedfile)
#             print(inputfile, " \ ",compareAD(M_g,M_r), " \ ", compareDF(M_g,M_r)," \ ", "Elapsed time: " , t_end-t_start," \ ", "n of cpus: ", n_cpu, file=Scorelog )
#             #print(inputfile, " \ ",compareAD(M_g,M_r), " \ ", compareDF(M_g,M_r), "Completed !...")
#         else:
#             print(inputfile, " Reconstruction NOT found ! ", file=Scorelog )

