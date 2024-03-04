#!/usr/bin/python

## define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
core_dir = home+'/repositories/QSAR_DAT-hERG'
conf_dir = core_dir+"/conf"
sys.path.insert(0, core_dir)
sys.path.insert(0, conf_dir)

## import all
from filters import *
from buildmodel import *
from misc import *

import getopt





def main(argv):


    try:
        opts, args = getopt.getopt(argv,"hp:a:t:o:b:s:",["protein=","tsv=","outdir=","basename=","standard_type"])

    except getopt.GetoptError:
        print(sys.argv[0]+' -t <tsv_file>  -p <target_protein> -a <assaydefinition> -o <output_dir> -b <output_file_basename> -s <standard_type>')
        sys.exit(2)

    for opt, arg in opts:
        if opt == '-h':
            print(sys.argv[0]+' -t <tsv_file>  -p <target_protein> -a <assaydefinition> -o <output_dir> -b <output_file_basename> -s <standard_type>')
            sys.exit()
            
        elif opt in ("-p", "--protein"):
            target = arg
            
        elif opt in ("-a", "--assaydefinition"):
            assaydefinition = arg

        elif opt in ("-t", "--tsv"):
            chembl_tsv_file = arg
            
        elif opt in ("-o", "--outdir"):
            output_dir = arg
            
        elif opt in ("-b", "--basename"):
            base_name = arg
            
        elif opt in ("-s", "-standard_type"):
            standard_type = arg

    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d%H%M%S")
    
    
    ## check target is defined
    try: target
    except NameError: target = ""
    if len(target) == 0:
        print("error: Need to specify target protein DAT or hERG")
        print(sys.argv[0]+' -t <tsv_file>  -p <target_protein> -a <assaydefinition> -o <output_dir> -b <output_file_basename> -s <standard_type>')
        sys.exit()
    if target != "DAT":
        if target != "hERG":
            print("error: Need to specify target protein DAT or hERG")
            print(sys.argv[0]+' -t <tsv_file>  -p <target_protein> -a <assaydefinition> -o <output_dir> -b <output_file_basename> -s <standard_type>')
            sys.exit()

    ## check assaydefinition is defined
    try: assaydefinition
    except NameError: assaydefinition = ""
    if len(assaydefinition) == 0:
        print("error: Need to specify assaydefinition DAT (uptake/inhibitor) or hERG(clamp/binding/others)")
        print(sys.argv[0]+' -t <tsv_file>  -p <target_protein> -a <assaydefinition> -o <output_dir> -b <output_file_basename> -s <standard_type>')
        sys.exit()
    if assaydefinition != "uptake":
        if assaydefinition != "inhibitor":
            if assaydefinition != "clamp":
                if assaydefinition != "binding":
                    if assaydefinition != "others":
                        print("error: Need to specify assaydefinition DAT (uptake/inhibitor) or hERG(clamp/binding/others)")
                        sys.exit()            
            
    ## check standard type Ki or IC50 is definded
    try: standard_type
    except NameError: standard_type = ""
    if len(standard_type) == 0:
        print("error: Need to specify standard_type Ki or IC50")
        sys.exit()
    if standard_type != "Ki":
        if standard_type != "IC50":
            print("error: Need to specify standard_type Ki or IC50")
            sys.exit()

            
            
    ## check tsv source exist
    output_dir=output_dir+"_"+standard_type
    if not os.path.isfile(chembl_tsv_file):
        #print("error: tsv file doesn't exist")
        sys.exit()

    ## check output directory exists
    if not os.path.exists(output_dir):
        #print("creating "+output_dir)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
    else:
        os.rename(output_dir, output_dir+"_old_"+timestamp)
        Path(output_dir).mkdir(parents=True, exist_ok=True)
        #print("move "+output_dir+" to "+output_dir+"_old_"+timestamp)
        #print("creating "+output_dir)            

    return (target, assaydefinition, chembl_tsv_file, output_dir, base_name, standard_type)





if __name__ == "__main__":
    
    target,assaydefinition,chembl_tsv_file,output_dir,base_name,standard_type = main(sys.argv[1:])
    
    buffer = read_data(chembl_tsv_file, Verbose = True)
    buffer = filter_confidence(buffer, Verbose = True)
    buffer = filter_assay_type(buffer, Verbose = True)
    if standard_type == "Ki":
        buffer = filter_affinity(buffer, Verbose = True, keepIC50=False, keepKi=True)
    if standard_type == "IC50":
        buffer = filter_affinity(buffer, Verbose = True, keepIC50=True, keepKi=False)
    buffer = filter_units(buffer, Verbose = True)
    filtered_out = filter_exact(buffer, Verbose = True)

        
    
    
    if target == "DAT":
        ## for DAT:using filter_assaydefinition_DAT 
        ## filter_secondary_test_set is only for hERG nad it is not applied to DAT. 
        if assaydefinition == "inhibitor":
            keys=["inhibitor"]
        elif assaydefinition == "uptake":
            keys=["uptake" ]
            
        for key in keys:
            filtered_in, filtered_out = filter_assaydefinition(filtered_out, target, key, Verbose = False)
            buffer = filtered_in
        print("Number of compounds after Displacement Assay filter:", len(buffer))
        
        #buffer = filter_secondary_test_set(buffer, Verbose = True) # only do this for hERG/clamp
        print("Number of compounds after removing testset 2 compounds:  n/a")
        
    if target == "hERG":
        ## for hERG: using filter_assaydefinition_hERG
        #buffer = filter_assaydefinition_hERG(buffer, Verbose = True)

        if assaydefinition == "clamp":
            keys=["clamp"]
        elif assaydefinition == "binding":
            keys=["clamp", "binding" ]
        elif assaydefinition == "others":
            keys=["clamp", "binding" ]

        for key in keys:
            filtered_in, filtered_out = filter_assaydefinition(filtered_out, target, key, Verbose = False)

        if assaydefinition == "clamp":
            buffer = filtered_in
        elif assaydefinition == "binding":
            buffer = filtered_in
        elif  assaydefinition == "others":
            buffer = filtered_out

        if assaydefinition == "clamp":
            print("Number of compounds after Displacement Assay filter:", len(buffer))
            buffer = filter_secondary_test_set(buffer, Verbose = True)
        else:
            print("Number of compounds after Displacement Assay filter:", len(buffer))
            #buffer = filter_secondary_test_set(buffer, Verbose = True) # only do this for hERG/clamp
            print("Number of compounds after removing testset 2 compounds:  n/a")
            
    buffer = filter_small_sets(buffer, Verbose = True, threshold=4)
    buffer = filter_salts(buffer, Verbose = True)
    buffer = filter_elements(buffer, Verbose = True)
    buffer = filter_size(buffer, Verbose = True)
    buffer = filter_pchembl_values(buffer, Verbose = True, replace=True)
    buffer = filter_weirdos(buffer, Verbose = True)
    buffer = deduplicate_mols(buffer, Verbose = True)
    write_smi_act_reg(buffer, base_name, output_dir, add_extra=False)
    write_smi_act_class(buffer, base_name, output_dir, inact_val = 5.0, act_val = 6.0, Verbose = True)

