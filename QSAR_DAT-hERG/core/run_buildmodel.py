#!/usr/bin/python

## import all
from filters import *
from buildmodel import *
from misc import *

import getopt


## define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
core_dir = home+'/repositories/QSAR_DAT-hERG'
conf_dir = core_dir+"/conf"
sys.path.insert(0, core_dir)
sys.path.insert(0, conf_dir)


# initialize global variables
stages = ['buildmodel', 'prediction', 'both']


def show_msg(arg):
    msg = ''' 
        
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

        '''

    print(arg + msg)
    sys.exit()


def main(argv):

    in_validation_set = None
    in_filename = None

    try:
        opts, args = getopt.getopt(argv,"hs:m:x:t:r:n:i:d:",
                                   ["stage=", "model=", "method=", "testing_percent=", "n_rand_stat", "num_splits=",
                                    "in_filename=", "in_validation_set="])

    except getopt.GetoptError:
        show_msg(sys.argv[0])

    if len(opts) == 0:
        show_msg(sys.argv[0])

    for opt, arg in opts:
        if opt == '-h':
            show_msg(sys.argv[0])
        elif opt in ("-s", "--stage"):
            stage = arg.lower()
            if stage not in stages:
                show_msg("Invalid stage option")
        elif opt in ("-m", "--model"):
            mode = arg.lower()
        elif opt in ("-x", "--method"):
            method = arg.lower()
        elif opt in ('-t', '--tp'):      # training dataset ratio
            tp = float(arg)
            if tp < 0 or tp >= 1:
                show_msg("-t option must be set as a decimal between 0 inclusive and 1")
        elif opt in ('-r', '--n_rand_stat'):  # number of random_stat for machine learning XGB/RF
            n_rand_stat = int(arg)
            if n_rand_stat < 1:
                show_msg("-r option must be at least 1")
        elif opt in ('-n', '--num_splits'):   # number of random_splits for splitting dataset into training/testing
            num_splits = int(arg)
            if num_splits < 1:
                show_msg("-n option must be at least 1")
        elif opt in ("-i", "--in_filename"):   # input filename, eg. dir/filename
            in_filename = arg
            if not os.path.exists(get_dir(arg)):
                show_msg("Input filename path does not exist")
        elif opt in ("-d", "--in_validation_set"):
            in_validation_set = arg

    if in_validation_set is None and stage.startswith('p'):
        show_msg('Need -d option')
        sys.exit()
    if in_filename is None and stage.startswith('b'):
        show_msg('Need -i option')
        sys.exit()

    return stage, mode, method, tp, num_splits, n_rand_stat, in_filename, in_validation_set


if __name__ == "__main__":
    # return argvs 
    stage, mode, method, tp, num_splits, n_rand_stat, in_filename, in_validation_set = main(sys.argv[1:])

    do_bm, do_pred = True, True
    if stage == stages[0]:
        do_pred = False
    elif stage == stages[1]:
        do_bm = False

    if do_bm:

        # define output_dir for output & log file
        output_dir = get_output_dir(mode, method, tp)
        writer = MyWriter(sys.stdout, output_dir+'.log')
        sys.stdout = writer

        # check required files/directories
        in_dataset_dir = get_dir(in_filename)
        check_misc(in_dataset_dir)
        check_required(in_filename, output_dir)

        # number of rand_splits/rand_states
        rand_splits = gen_random_splits(control_seed=2020, num_splits=num_splits)
        if tp == 0.0:
            rand_splits = [rand_splits[0]]
        rand_states = gen_random_splits(control_seed=501, num_splits=n_rand_stat)

        # loop do build models
        for rand_state in rand_states:
            for random_split in rand_splits:
                print(random_split)
                output_ext = get_output_ext(mode, method, tp, rand_state, random_split)

                mols, acts, deletes, changes = read_data4buildmodel(in_filename, mode)
                mols, acts = curate_mols(mols, acts, deletes, changes)

                # using split_data (tp!=0.00) or all_data (tp=0.00)
                if tp > 0.00:
                    train_mols, train_names, train_acts, \
                    test_mols, test_names, test_acts = split_data(mols, acts, tp, random_split)
                else:
                    train_mols, train_names, train_acts = all_data(mols, acts)

                ad_fps, ad_rad = calc_appdom(train_mols, output_dir)

                if tp > 0:
                    test_mols, test_acts, test_names, test_mols_reject, test_acts_reject, test_names_reject \
                        = check_appdom(ad_fps, ad_rad, test_mols, test_names, test_acts, step=stage)

                # feature_topo
                train_topo_descs = calc_topo_descs(train_mols)
                train_topo_descs, topo_index, topo_names = prune_topo_descs(mode, train_topo_descs, train_acts, output_dir)

                if method == 'xgb' and mode.startswith('reg'):
                    train_phore_descs = calc_phore_descs(train_mols)
                    train_phore_descs, phore_sigbits, phore_names = prune_phore_descs(train_phore_descs, output_dir)
                    train_descs = np.concatenate((train_topo_descs, train_phore_descs), axis=1)
                    if tp > 0:
                        test_topo_descs = calc_topo_descs(test_mols, topo_index)
                        test_phore_descriptors = calc_phore_descs(test_mols, phore_sigbits)
                        test_descs = np.concatenate((test_topo_descs, test_phore_descriptors), axis=1)
                else:
                    train_descs = train_topo_descs
                    if tp > 0:
                        test_topo_descs = calc_topo_descs(test_mols, topo_index)
                        test_descs = test_topo_descs

                # build_model
                model, model_score, best_params = build_model(mode, method, rand_state, train_descs, train_acts, output_dir)

                if tp > 0:
                    # test_mols, test_acts, test_names, rej_molecules, rej_activities, rej_names = check_appdom(ad_fps, ad_rad, test_mols, test_acts, test_names, step=stage)
                    # test_descriptors = calc_topo_descs(test_molecules, topo_index)
                    if mode.startswith('reg'):
                        test_r2, test_rmse, test_mse = predict_model(model, test_descs, test_acts, train_acts, random_split, output_dir, mode, method, rand_state)
                    elif mode.startswith('class'):
                        predictions, acc = predict_model(model, test_descs, test_acts, train_acts, random_split, output_dir, mode, method, rand_state)
                    # TODO: do we really need those return?

    if do_pred:

        print("prediction")

        # hacking ARGV to old format
        in_model_dir = get_output_dir(mode, method, tp)
        splits = in_validation_set.split('/')
        in_validation_dir = get_dir(in_validation_set)
        in_validation_filename = splits[-1]

        #  input_basename== testset_pairs
        #  data_dir

        output_dir = 'pred_' + get_output_dir(mode, method, tp)
        result_list = []
        # number of rand_state

        rand_states = gen_random_splits(control_seed=501, num_splits=n_rand_stat)
        rand_splits = gen_random_splits(num_splits=num_splits)
        for random_state in rand_states:
            for random_split in rand_splits:
                output_ext = get_output_ext(mode, method, tp, random_state, random_split)

                input_data = read_mols(mode, method, in_validation_filename, datadir=in_validation_dir, modeldir=in_model_dir)
                molnames = input_data['molnames']
                mols = input_data['molecules']
                model = input_data['model']
                inds = input_data['inds']

                if mode.startswith('reg') and method == 'xgb':
                    sigbits = input_data['sigbits']
                ad_fps = input_data['ad_fps']
                ad_radius = input_data['ad_radius']
                # Check Applicability Domain
                appdom_results = check_appdom(ad_fps, ad_radius, mols, molnames, step=stage)
                mols = appdom_results['test_mols']
                molnames = appdom_results['test_names']
                molecules_rej = appdom_results['rej_mols']
                molnames_rej = appdom_results['rej_names']
                descriptors = calc_topo_descs(mols, inds)
                if mode.startswith('reg'):
                    if method == 'xgb':
                        phore_descriptors = calc_phore_descs(mols, sigbits)
                        descriptors = np.concatenate((descriptors, phore_descriptors), axis=1)
                    pred_results = make_preds(molnames, descriptors, model, random_split, mode=mode)
                    result_list.append(pred_results['predictions'])
                elif mode.startswith('class'):
                    pred_results = make_preds(molnames, descriptors, model, mode=mode)
                    if result_list == []:
                        result_list = pred_results['predictions']
                    else:
                        result_list += pred_results['predictions']

            if mode.startswith('reg'):
                compound, pred_mean, pred_error = summarize_preds(molnames, result_list)
                data = {"compound": compound, "mean": pred_mean, "stdev": pred_error}
            elif mode.startswith('class'):
                activity = []
                for prediction in result_list:
                    if prediction > 0.5 * len(rand_splits):
                        activity.append(1)
                    else:
                        activity.append(0)
                data = {"Compound": molnames, "Sum_activity": list(result_list), "Pred_activity": activity}

            df = pd.DataFrame(data=data)

            # makes the output directory if it does not already exist
            if not os.path.exists(output_dir):
                os.makedirs(output_dir)

            ## TODO, the following will case issue. if you have input directory as tetst/, it will return nothing as output_filename
            output_filename = in_model_dir.split('/')[-1]

            df.to_csv(output_dir + '/' + output_filename, header=True, index=False, sep='\t')

