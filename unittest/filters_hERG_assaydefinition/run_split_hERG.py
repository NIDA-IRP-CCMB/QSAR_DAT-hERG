## define enviroment
import sys,os
currentdir = os.getcwd()
core_dir = currentdir+'/core'
conf_dir = currentdir+"/conf"
sys.path.insert(0, core_dir)


## import all
from filters import *
target = "hERG"


run_file = "indata/chembl25_hERG.tsv"
buffer = read_data(run_file, Verbose = True)
buffer = filter_confidence(buffer, Verbose = True)
buffer = filter_assay_type(buffer, Verbose = True)
buffer = filter_affinity(buffer, Verbose = True, keepIC50=True, keepKi=False)
buffer = filter_units(buffer, Verbose = True)
filtered_out = filter_exact(buffer, Verbose = True)

keys=["clamp", "binding", "other"]
for key in keys:
    print(key)
    filterfile=conf_dir+'/assaydefinition_hERG_'+key+'.txt'
    filtered_in, filtered_out = filter_assaydefinition(filtered_out, target, key, Verbose = True)

filtered_out[['description','pubmed_id','doi']].to_csv("hERG_data_out.dat",sep='\t', index=False)
filtered_out.to_csv("hERG_data_out.tsv",sep='\t', index=False)
