## 1. Datasets (data_sql) provided in this example

    # dataset provide in data_sql fodler

        data_sql
        ├── DAT
        │   ├── chembl25_homo.tsv			<-------- dataset for human DAT
        │   ├── chembl25_rat.tsv			<-------- dataset for rat DAT
        │   ├── chembl25_raw.tsv
        │   ├── chembl25.sql			<-------- sql
        │   └── chembl25.tsv			<-------- dataset for all DAT
        └── hERG
            ├── chembl25_raw.tsv
            ├── chembl25.sql                        <-------- sql
            └── chembl25.tsv			<-------- dataset for hERG


### Steps to generate initial ChEMBL dataset from our inhouse sql server

Use psql and chembl25_Classifier.sql to generate chembl25_Classifier_raw.tsv
   
    # You would need the chembl password for running our inhouse sql server. Please ask the admin for the password.   
    /usr/pgsql-9.6/bin/psql -h 127.0.0.1 chembl_25 chembl  --field-separator=$'\t'  --no-align -f chembl25_Classifier.sql -o chembl25_Classifier_raw.tsv

Clean up the tsv file. remove the header and the last line of the file

    sed '1d;$d' chembl25_Classifier_raw.tsv > chembl25_Classifier.tsv

Extract human/rat bioactivity for DAT

    cat chembl25_Classifier.tsv | awk '$3=="Homo" {print}' > chembl25_Classifier_homo.tsv
    cat chembl25_Classifier.tsv | awk '$3=="Rattus" {print}' > chembl25_Classifier_rat.tsv





## 2. Examples of running filters and buildmodels

All output files in the example can be found in output. 

### load python modules
Note that you would need to have necessary modules installed in your python environment.  

    module load python

### defince COREPATH

    COREPATH="/path/to/core/dir"

### Usage of run_filters.py

    python ${COREPATH}/run_filters.py -help

    run_filters.py -t <tsv_file>  -p <target_protein> -a <assaydefinition> -o <output_dir> -b <output_file_basename> -s <standard_type>

### An example of running filters

    python ${COREPATH}/run_filters.py -p DAT -t data_sql/DAT/chembl25.tsv -a inhibitor -o dataset_all_DAT_inhibitor -b pubdata -s Ki | tee dataset_all_DAT_inhibitor_Ki.log


    touch dataset_all_DAT_inhibitor/AmyCompounds.act
    touch dataset_all_DAT_inhibitor/AmyCompounds.smi
    touch dataset_all_DAT_inhibitor/to_change.txt
    touch dataset_all_DAT_inhibitor/to_remove.txt


### Usage of run_buildmodel.py

    python ${COREPATH}/run_buildmodel.py -help
    
    run_buildmodel.py 
            
            -s <stage>              buildmodel, prediction, or both
            -m <model>              regression or classification model
            -x <method>             xgb or rf (random forest) machine learning methods
            -t <testing_percent>    proportion of dataset to reserve for testing; decimal range [0,1)
            -r <random_state>       number of random states/seeds to use to build the model; int range [1,inf)
            -n <num_splits>         number of times to randomly split your data into a testing set; int range [1,inf)
            -i <input_basename>     location of your input data with its basename (without ext ex. .act or .smi)
            -d <validation_data>    dataset used to benchmark the model prediction
    
            If -t is set to 0, then -n must be set to 1 as there is only one way to split the data into 100% training
            and 0% testing data. If -n is greater than 1, then -t must be in range (0,1).
            
            If stage is set to do prediction then -d <validation_data> must also be given. Can input nothing for -d if build
            model.


### An example of running buildmodels

    python ${COREPATH}/run_buildmodel.py -s buildmodel -m reg -x xgb   -t 0.15 -r 1 -n 1 -i dataset_all_DAT_inhibitor_Ki/pubdata

### An example of making prediction

    # generate testing compounds for example prediction
    head -n 10 dataset_all_DAT_inhibitor_Ki/pubdata.smi > testcoumpounds.smi

    # run prediction
    python $COREPATH/run_buildmodel.py -s prediction -m reg -x xgb   -t 0.15 -r 1 -n 1 -d ./testcoumpounds


