## define enviroment
import sys,os
currentdir = os.getcwd()
core_dir = currentdir+'/core'
conf_dir = currentdir+"/conf"
sys.path.insert(0, core_dir)


## import all
from filters import *


run_file = "indata/chembl25_DAT.tsv"
buffer = read_data(run_file, Verbose = True)
buffer = filter_confidence(buffer, Verbose = True)
buffer = filter_assay_type(buffer, Verbose = True)
buffer = filter_affinity(buffer, Verbose = True, keepIC50=False, keepKi=True)
buffer = filter_units(buffer, Verbose = True)
buffer = filter_exact(buffer, Verbose = True)


print("inhibitor:")
A = filter_assaydefinition_DAT(buffer, Verbose = True)
B = filter_assaydefinition_DAT2(buffer, Verbose = False)

A[['description','pubmed_id','doi']].to_csv("DAT_data_inhibitor.dat",sep='\t', index=False)
A.to_csv("DAT_data_inhibitor.tsv", sep='\t', index=False)




print("uptake:")
C = filter_assaydefinition_DAT_uptake(B, Verbose = True)
print("out:")
D = filter_assaydefinition_DAT_uptake2(B, Verbose = True)


C[['description','pubmed_id','doi']].to_csv("DAT_data_uptake.dat",sep='\t', index=False)
C.to_csv("DAT_data_uptake.tsv", sep='\t', index=False)
D[['description','pubmed_id','doi']].to_csv("DAT_data_out.dat",sep='\t', index=False)
D.to_csv("DAT_data_out.tsv", sep='\t', index=False)


#buffer[['description','pubmed_id','doi']].to_csv("DAT_data_keep.dat",sep='\t', index=False)
