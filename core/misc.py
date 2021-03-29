import numpy as np
from pathlib import Path
import datetime


## define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
core_dir = home+'/repositories/herg/core'
conf_dir = core_dir+"/conf"
sys.path.insert(0, core_dir)
sys.path.insert(0, conf_dir)


# generate random_splits
def gen_random_splits(control_seed=2020, num_splits=50, Verbose = False):
    random_splits = []
    np.random.seed(control_seed)

    for a in range(0,num_splits):
        new_random = np.random.randint(0,1000)
        #while (new_random in random_splits_old) or (new_random in random_splits):
        while (new_random in random_splits):
            #ensure no duplicate 
            new_random = np.random.randint(0,1000)
        random_splits.append(new_random)

    #random_splits
    if Verbose:
        for random_split in random_splits:
            print(random_split)
    return random_splits


# check required files exist before doing anything
def check_required(input_basename,output_model_dir):
    now = datetime.datetime.now()
    timestamp = now.strftime("%Y%m%d%H%M%S")

    if not os.path.exists(output_model_dir):
        Path(output_model_dir).mkdir(parents=True, exist_ok=True)
    else:
        os.rename(output_model_dir, output_model_dir+"_old_"+timestamp)
        Path(output_model_dir).mkdir(parents=True, exist_ok=True)

    if not os.path.isfile(input_basename+".smi"):
        print("error: smi doesn't exist")

    if not os.path.isfile(input_basename+".act"):
        print("error: act doesn't exist")
    
    return


class MyWriter:

    def __init__(self, stdout, filename):
        self.stdout = stdout
        self.logfile = open(filename, 'a')

    def write(self, text):
        self.stdout.write(text)
        self.logfile.write(text)

    def close(self):
        self.stdout.close()
        self.logfile.close()

    def flush(self):
        pass


def get_dir(input_base):
    splits = input_base.split('/')
    output = ''
    for i in range(len(splits)-1):
        output += splits[i] + '/'
    return output


# Function used to check directory for required misc files
# We define AmyCompounds.* and to_* in conf file
def check_misc(input_dir):
    if os.path.exists(input_dir):
        with open(core_dir+'/conf/required_files.txt', "r") as f:
            for line in f:
                line = line.rstrip("\n")
                if not os.path.isfile(input_dir+line):
                    print('Missing %s data in data directory' % line)
                    sys.exit(0)
    else:
        print('Error: input directory does not exist')

