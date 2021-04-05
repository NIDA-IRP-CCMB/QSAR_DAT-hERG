"""
Program: buildmodel.py

Written: Andrew Fant, Joslyn Jung

This is a python program to generate a QSAR regression model using rdkit and scipy. The molecular descriptor data is
    read and processed with extraneous descriptors removed, and the XGBoost model determines the best-performing
    parameters to be used.

input: four files.  1) named training.smi. This is the Molecules training set.
                       SMILES formatted holding the structure of the molecules in the model training set.

                    2) named training.act. This is the Activities training set.  each line has the name of one molecule
                       and its activity as a pIC50.

                    3) named testset-2.smi.
                    4) named testset-2.act.

                   N.B. molecules and activities need to be in the same order in both files.

output: generates and saves XGBoost regression model with optimized parameters in build.pickle.dat
"""

## define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
base_dir = home+'/repositories/herg/hERGvDAT/'
core_dir = base_dir+'/core'
conf_dir = core_dir+'/conf'
sys.path.insert(0, conf_dir)
sys.path.insert(0, core_dir)



import os
import numpy as np
from rdkit.Chem import AllChem as Chem
from rdkit.Chem.Pharm2D import Gobbi_Pharm2D, Generate
from rdkit import DataStructs
from rdkit.Chem.Fingerprints import FingerprintMols
from molvs import Standardizer
from scipy.stats.stats import pearsonr
from sklearn.model_selection import train_test_split, GridSearchCV
from sklearn.metrics import roc_auc_score, r2_score, mean_squared_error, confusion_matrix, recall_score, accuracy_score
from sklearn.ensemble import RandomForestRegressor, RandomForestClassifier
import operator
import copy
from descriptor_setup import dnames, dlist
from xgboost import XGBRegressor, XGBClassifier
from math import sqrt
import pickle
import matplotlib.pyplot as plt
import matplotlib
from misc import *

matplotlib.use('Agg')

# Initialize global variables
output_ext = ''


#def read_data4buildmodel(data_file, datadir, mode='regression'):
def read_data4buildmodel(input_basename, mode='regression'):

    # Read in molecules and activities
    molfh = open(input_basename+'.smi')
    actfh = open(input_basename+'.act')

    datadir = get_dir(input_basename)

    input_molecules = []
    input_activities = []

    delete_list = []
    change_dict = {}

    act = []

    for molline in molfh:
        line = molline[:-1].split()
        mol = Chem.MolFromSmiles(line[0])
        molname = line[1]
        input_molecules.append((mol, molname))

    for actline in actfh:
        line = actline[:-1].split()
        if mode.startswith('reg'):
            act = float(line[1])
        elif mode.startswith('class'):
            act = int(line[1])
        actname = line[0]
        input_activities.append((act, actname))

    molfh.close()
    actfh.close()

    if mode.startswith('reg'):
        molfh2 = open(datadir + "/AmyCompounds.smi")
        actfh2 = open(datadir + "/AmyCompounds.act")

        moldeletefh = open(datadir + "/to_remove.txt")
        molchangefh = open(datadir + "/to_change.txt")

        for molline in molfh2:
            line = molline.split()
            mol = Chem.MolFromSmiles(line[0])
            molname = line[1]
            input_molecules.append((mol, molname))

        for actline in actfh2:
            line = actline.split()
            act = float(line[1])
            actname = line[0]
            input_activities.append((act, actname))

        molfh2.close()
        actfh2.close()

        # Molecules to remove or to change the affinities of
        for delline in moldeletefh:
            if '#' in delline:
                line = delline.split('#')[0]
            else:
                line = delline
            delitem = line[:-1]
            delete_list.append(delitem)

        for changeline in molchangefh:
            if '#' in changeline:
                line = changeline.split('#')
            else:
                line = changeline
            workline = line.split()
            newvalue = float(workline[1])
            actname = workline[0]
            change_dict[actname] = newvalue

        moldeletefh.close()
        molchangefh.close()

    return input_molecules, input_activities, delete_list, change_dict


def curate_mols(input_mols, input_acts, delete_list, change_dict):
    output_mols = []
    output_acts = []
    if len(input_mols) != len(input_acts):
        print("Error: Number of molecules not equal to number of activities")
        exit(1)
    for item in range(len(input_mols)):
        if input_mols[item][1] in delete_list:
            continue
        curr_name = input_acts[item][1]
        curr_mol = input_mols[item][0]
        curr_act = input_acts[item][0]
        if input_acts[item][1] in change_dict.keys():
            curr_act = change_dict[curr_name]
        output_mols.append((curr_mol, curr_name))
        output_acts.append((curr_act, curr_name))

    if len(output_mols) != len(output_acts):
        print("Error: curated molecules length not equal to curated activities length")
        exit(1)
    for i in range(len(output_mols)):
        if output_mols[i][1] != output_acts[i][1]:
            print("Error: curated molecules order is not the same as curated activities order")
            exit(1)
    return output_mols, output_acts


def split_data(mols, acts, test_percent, split):
    mols_train = []
    mols_test = []
    molnames_train = []
    molnames_test = []
    acts_train = []
    acts_test = []
    actnames_train = []
    actnames_test = []

    # Split molecules and activities training set into training and test sets
    m_train, m_test, a_train, a_test = train_test_split(mols, acts, test_size=test_percent, random_state=split)
    # Make a list of the names of all the molecules in the training list
    names_train = []

    for mol in m_train:
        names_train.append(mol[1])

    # Iterate over all the molecules we have read in
    for i in range(len(mols)):
        # assert mols[i][1] == acts[i][1]
        if mols[i][1] in names_train:  # is the molecule in the training set?
            mols_train.append(mols[i][0])
            molnames_train.append(mols[i][1])
            acts_train.append(acts[i][0])
            actnames_train.append(acts[i][1])
        else:  # the molecule is in the test set if it isn't in the the training set
            mols_test.append(mols[i][0])
            molnames_test.append(mols[i][1])
            acts_test.append(acts[i][0])
            actnames_test.append(acts[i][1])

    assert molnames_train == actnames_train
    assert molnames_test == actnames_test

    # Standardize structures of the training set and test set
    s = Standardizer()
    standard_mols_train = []

    for mol in mols_train:
        standard_mols_train.append(s.standardize(mol))

    standard_mols_test = []

    for mol in mols_test:
        standard_mols_test.append(s.standardize(mol))

    return standard_mols_train, molnames_train, acts_train, standard_mols_test, molnames_test, acts_test


def all_data(molecules, activities):
    training_molecules = []
    training_names = []
    training_activities = []
    for i in range(len(molecules)):
        training_molecules.append(molecules[i][0])
        training_names.append(molecules[i][1])
        training_activities.append(activities[i][0])
    return training_molecules, training_names, training_activities
       
    
def get_output_ext(mode, ml, tp, rand_state, rand_split):
    global output_ext
    output_ext = "%s_%s_%.2f_%d_%d" % (mode, ml, float(tp), int(rand_state), int(rand_split))
    return output_ext


def get_output_dir(mode, ml, tp):
    global output_dir
    output_dir = "%s_%s_%.2f" % (mode, ml, float(tp))
    return output_dir


def calc_appdom(training_set, out_model_dir):
    appdom_fps = []
    # output_ext = "%s_%s_%d_%d" % (mode, method, int(rand_split), int(rand_state))

    for mol in training_set:
        # fingerprint = Chem.GetMorganFingerprintAsBitVect(mol, 3, nBits=1024)
        fingerprint = FingerprintMols.FingerprintMol(mol)
        appdom_fps.append(fingerprint)

    distances = []

    for i in range(0, (len(appdom_fps) - 1)):
        for j in range(i + 1, len(appdom_fps)):
            # dist = 1.0 - (DataStructs.TanimotoSimilarity(appdom_fps[i], appdom_fps[j]))
            dist = 1.0 - (DataStructs.FingerprintSimilarity(appdom_fps[i], appdom_fps[j]))
            distances.append(dist)

    distances = np.array(distances)
    mean_distance = np.mean(distances)
    dev_distance = np.std(distances)

    appdom_radius = mean_distance + dev_distance

    # Write fingerprints of training set and AD radius to pickle files for later prediction runs
    with open(out_model_dir + "/training-FPs_%s.dat" % output_ext, 'wb') as f:
        pickle.dump(appdom_fps, f)

    with open(out_model_dir + "/AD-radius_%s.dat" % output_ext, 'wb') as f:
        pickle.dump(appdom_radius, f)

    return appdom_fps, appdom_radius


# args: fps, rad, mols, names, acts
def check_appdom(*args, step):
    if len(args) < 4 or len(args) > 5:
        print("Error: incorrect number of arguments passed to check_appdom()")
        exit(1)

    appdom_fps = args[0]
    appdom_radius = args[1]
    pred_mols = args[2]
    pred_names = args[3]
    if len(args) == 5:
        pred_acts = args[4]

    accept_mols = []
    accept_names = []
    reject_mols = []
    reject_names = []

    if 'pred_acts' in locals():
        accept_acts = []
        reject_acts = []

    for i in range(len(pred_mols)):
        test_fp = FingerprintMols.FingerprintMol(pred_mols[i])
        distances = []
        for training_fp in appdom_fps:
            distances.append(1.0 - (DataStructs.FingerprintSimilarity(training_fp, test_fp)))

        distances = np.array(distances)
        if np.min(distances) <= appdom_radius:
            accept_mols.append(pred_mols[i])
            accept_names.append(pred_names[i])
            if 'pred_acts' in locals():
                accept_acts.append(pred_acts[i])
        else:
            reject_mols.append(pred_mols[i])
            reject_names.append(pred_names[i])
            if 'pred_acts' in locals() and len(args) > 4:
                reject_acts.append(pred_acts[i])

            print("Compound %s is out of the AD for this model", pred_names[i])

    if step.startswith('p'):
        if len(reject_names) == 0:
            print("No molecules rejected for prediction by AD")
        return_dict = {}

        return_dict['test_mols'] = accept_mols
        return_dict['test_names'] = accept_names
        return_dict['rej_mols'] = reject_mols
        return_dict['rej_names'] = reject_names

        if 'pred_acts' in locals() and len(args) > 4:
            return_dict['test_acts'] = accept_acts
            return_dict['rej_acts'] = reject_acts

        return return_dict
    elif step.startswith('b'):
        if len(reject_names) == 0:
            print('All molecules in test set within applicability domain.')
        return accept_mols, accept_acts, accept_names, reject_mols, reject_acts, reject_names


def calc_topo_descs(mols, indexes=None):
    descs = np.zeros((len(mols), len(dlist)))
    for i in range(len(mols)):
        for j in range(len(dlist)):
            descs[i, j] = dlist[j](mols[i])
    if indexes is not None:
        # Select descriptors
        del_indexes = []
        for i in range(len(dlist)):
            if i not in indexes:
                del_indexes.append(i)

        del_indexes.reverse()

        for i in del_indexes:
            descs = np.delete(descs, [i], axis=1)

    return descs


def prune_topo_descs(mode, input_descs, acts_train, out_model_dir):
    # output_ext = "%s_%s_%d_%d" % (mode, method, int(rand_split), int(rand_state))
    local_dnames = copy.copy(dnames)  # create a local copy of dnames
    del_indices = []  # array of indices of deleted entries from dnames

    hold = np.copy(input_descs)

    # Remove zero-variance descriptors
    for i in range((len(local_dnames) - 1), -1, -1):
        if min(input_descs[:, i]) == max(input_descs[:, i]):
            hold = np.delete(hold, [i], axis=1)
            del_indices.append(i)
            del local_dnames[i]
    descs = hold

    # Remove highly correlated descriptors
    correl = np.zeros((len(local_dnames), len(local_dnames)))
    hcpairs = []
    hcdescs = {}
    descs_to_del = []

    for i in range(len(local_dnames)):
        for j in range(len(local_dnames)):
            correl[i, j] = pearsonr(descs[:, i], descs[:, j])[0]

    for i in range((len(local_dnames) - 1)):
        for j in range(i + 1, len(local_dnames)):
            if correl[i, j] > 0.99:
                hcpairs.append((i, j))

    for pair in hcpairs:
        if pair[0] not in hcdescs.keys():
            hcdescs[pair[0]] = 1
        else:
            hcdescs[pair[0]] = hcdescs[pair[0]] + 1
        if pair[1] not in hcdescs.keys():
            hcdescs[pair[1]] = 1
        else:
            hcdescs[pair[1]] = hcdescs[pair[1]] + 1

    sorted_hcdescs = sorted(hcdescs.items(), key=operator.itemgetter(0))
    sorted_hcdescs.reverse()

    while len(sorted_hcdescs) > 0:
        foo = sorted_hcdescs[0][0]
        descs_to_del.append(foo)

        for i in range((len(hcpairs) - 1), -1, -1):
            if foo in hcpairs[i]:
                del hcpairs[i]

        hcdescs = {}
        for pair in hcpairs:
            if pair[0] not in hcdescs.keys():
                hcdescs[pair[0]] = 1
            else:
                hcdescs[pair[0]] = hcdescs[pair[0]] + 1
            if pair[1] not in hcdescs.keys():
                hcdescs[pair[1]] = 1
            else:
                hcdescs[pair[1]] = hcdescs[pair[1]] + 1

        sorted_hcdescs = sorted(hcdescs.items(), key=operator.itemgetter(1))
        sorted_hcdescs.reverse()

    descs_to_del.sort()
    descs_to_del.reverse()

    hold = np.copy(descs)

    del_indices = []
    for i in descs_to_del:
        hold = np.delete(hold, [i], axis=1)
        del_indices.append(i)
        del local_dnames[i]
    descs = hold
    if mode.startswith('reg'):
        # Select predictive variables for further modeling
        acts_train_binarized = copy.copy(acts_train)
        descs_to_del = []

        for foo in range(len(acts_train_binarized)):
            if float(acts_train_binarized[foo]) < 5:
                acts_train_binarized[foo] = 0
            else:
                acts_train_binarized[foo] = 1

        for desc_num in range(len(local_dnames)):
            aroc = roc_auc_score(acts_train_binarized, descs[:, desc_num])
            if aroc < 0.55 and aroc > 0.45:
                descs_to_del.append(desc_num)

        descs_to_del.reverse()
        hold = np.copy(descs)

        del_indices = []
        for i in descs_to_del:
            hold = np.delete(hold, [i], axis=1)
            del_indices.append(i)
            del local_dnames[i]
    cleaned_descriptors = hold

    # Get complete list of indices of descriptors in use
    indices = []

    for desc in local_dnames:
        i = dnames.index(desc)
        indices.append(i)
    indices.sort()

    # Write list of descriptor indices and significant bits to files
    with open(out_model_dir + "/indices_%s.dat" % output_ext, 'wb') as f:
        pickle.dump(indices, f)

    return cleaned_descriptors, indices, local_dnames


def calc_phore_descs(mols, significant_bits=None, testing=False):
    fp_holding = []
    accumulated_bits_set = {}

    for mol in mols:
        fp = Generate.Gen2DFingerprint(mol, Gobbi_Pharm2D.factory)
        fp_holding.append(fp)
        if significant_bits is not None:
            bits_set = list(fp.GetOnBits())
            for fp_bit in bits_set:
                if fp_bit not in accumulated_bits_set.keys():
                    accumulated_bits_set[fp_bit] = 1
                else:
                    accumulated_bits_set[fp_bit] = accumulated_bits_set[fp_bit] + 1

    if significant_bits is not None:
        phore_descs = np.zeros((len(mols), len(significant_bits)))

        for mol_num in range(len(mols)):
            for bit_num in range(len(significant_bits)):
                if significant_bits[bit_num] in fp_holding[mol_num].GetOnBits():
                    phore_descs[mol_num, bit_num] = 1
        if testing:
            return "significant_bits: %d" % len(significant_bits), "fp_descriptors: %s" % str(phore_descs.shape)
        print("significant_bits:", len(significant_bits))
        print("fp_descriptors:", phore_descs.shape)
        return phore_descs
    else:
        return fp_holding


def prune_phore_descs(input_descs, out_model_dir):
    # output_ext = "%s_%s_%d_%d" % (mode, method, int(rand_split), int(rand_state))
    # Extract pharmacophore bits that occur often enough to be useful
    accumulated_bits_set = {}

    for fp in input_descs:
        bits_set = list(fp.GetOnBits())
        for fp_bit in bits_set:
            if fp_bit not in accumulated_bits_set.keys():
                accumulated_bits_set[fp_bit] = 1
            else:
                accumulated_bits_set[fp_bit] = accumulated_bits_set[fp_bit] + 1

    significant_bits = []
    all_set_bits = list(accumulated_bits_set.keys())
    all_set_bits.sort()
    fp_names = []

    for bit in all_set_bits:
        if accumulated_bits_set[bit] >= 100:
            significant_bits.append(bit)
            bit_name = 'Ph2D_' + str(bit)
            fp_names.append(bit_name)

    nmols = np.shape(input_descs)[0]

    fp_descriptors = np.zeros((nmols, len(significant_bits)))

    for mol_num in range(nmols):
        for bit_num in range(len(significant_bits)):
            if significant_bits[bit_num] in input_descs[mol_num].GetOnBits():
                fp_descriptors[mol_num, bit_num] = 1

    print("significant_bits:", len(significant_bits))
    print("fp_descriptors:", fp_descriptors.shape)

    with open(out_model_dir + "/sigbits_%s.dat" % output_ext, 'wb') as f:
        pickle.dump(significant_bits, f)

    return fp_descriptors, significant_bits, fp_names


def build_model(mode, ml, rand_state, training_descs, training_acts, out_model_dir):
    # output_ext = "%s_%s_%d_%d" % (mode, method, int(rand_split), int(rand_state))
    if ml == 'xgb':
        # Parameter sets as specified in Noskov and Wacker
        # For parameters, random_state is a random seed, default is 0
        # If random_state is not changed, then the results will be the same
        parameters_xgb = dict(colsample_bytree=[0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                              subsample=[0.2, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1],
                              max_depth=[2, 3, 4, 5, 6, 7, 8],
                              learning_rate=[0.001, 0.01, 0.02, 0.08, 0.1])
        # Train XGBoost model and estimate best parameters using GridSearchCV with 10-fold cross validation
        if mode.startswith('reg'):
            mode = 'reg'
            reg = GridSearchCV(XGBRegressor(random_state=rand_state), parameters_xgb, cv=10,
                               scoring='neg_mean_squared_error', n_jobs=-1)
        elif mode.startswith('class'):
            mode = 'class'
            reg = GridSearchCV(XGBClassifier(random_state=rand_state), parameters_xgb, cv=10, scoring='accuracy',
                               n_jobs=-1)

        # Uncomment the following line for y-randomization tests only
        # random.shuffle(acts_train)

        reg.fit(training_descs, training_acts)
        best_params = reg.best_params_
        candidate_model = reg.best_estimator_
        score = reg.best_score_
        if mode.startswith('class'):
            # oob = model.oob_score_
            print("Best model params: ", best_params)
            print("Final Accuracy Score: ", score)
            # print("Out-Of-Bag error for this model is %5.3f" % oob)

        # Save model generated by GridSearchCV with XGBoost regressor to a file
        with open(str(out_model_dir) + "/model_%s.dat" % output_ext, "wb") as f:
            pickle.dump(candidate_model, f)
        if mode.startswith('reg'):
            print("Best parameters score:", str(best_params), abs(score))
    elif ml == 'rf':
        # Train RandomForest model and estimate best parameters using GridSearchCV with 10-fold cross validation
        if mode.startswith('reg'):
            mode = 'reg'
            parameters = {'n_estimators': [30, 100, 300, 1000, 3000, 10000, 30000, 100000]}
            mod = GridSearchCV(RandomForestRegressor(random_state=rand_state,oob_score=True),
                               parameters,
                               cv=10,
                               scoring='neg_mean_squared_error',
                               n_jobs=-1)
        elif mode.startswith('class'):
            mode = 'class'
            parameters = {'n_estimators': [30, 100, 300, 1000, 3000, 10000, 30000, 100000],
                          'max_features': ['auto', 'sqrt', 'log2']
                          }
            mod = GridSearchCV(RandomForestClassifier(random_state=rand_state,oob_score=True),
                               parameters,
                               cv=10,
                               scoring='accuracy',
                               n_jobs=-1)
        mod.fit(training_descs, training_acts)

        best_params = mod.best_params_
        candidate_model = mod.best_estimator_
        score = mod.best_score_
        oob = candidate_model.oob_score_

        print("Best model params: ", best_params)
        if mode.startswith('class'):
            print("Final Accuracy Score: ", score)
        else:
            print("Final MSE Score: ", score)
        print("Out-Of-Bag error for this model is %5.3f" % oob)

        # Save model generated by GridSearchCV with XGBoost regressor to a file
        with open(out_model_dir + "/model_%s.dat" % output_ext, "wb") as f:
            pickle.dump(candidate_model, f)

    return candidate_model, score, best_params


def predict_model(candidate_model, test_descs, test_acts, train_acts, split, out_model_dir, mode,
                  ml_type, rand_state, Verbose=False):
    if mode.startswith('reg'):
        # Make prediction on test set using model
        y_true, y_pred = test_acts, candidate_model.predict(test_descs)

        f = open(out_model_dir + "/stat_" + str(split), "w")
        f.write('random_split %s R2 %s RMSE %s MSE %s\n' % (
            split, R2(y_true, y_pred), RMSE(y_true, y_pred), MSE(y_true, y_pred)))
        f.close()

        if Verbose:
            # Evaluate model using measures of fit
            print("R2: %2.3f" % R2(y_true, y_pred))
            print("RMSE: %2.3f" % RMSE(y_true, y_pred))
            print("MSE: %2.3f" % MSE(y_true, y_pred))

        # Save model generated by GridSearchCV with XGBoost regressor to a file
        # pickle.dump(candidate_model,
        #             open(out_model_dir + "/model_reg_%s_%d_%d.dat" % (ml_type, split, rand_state), "wb"))

        # Plot data: y_pred vs. y_true
        plt.xlim(1, 10)
        plt.ylim(1, 10)
        plt.scatter(y_true, y_pred)
        plt.title("Test Data Results")
        plt.xlabel("True Output")
        plt.ylabel("Predicted Output")
        plt.savefig(out_model_dir + "/graph_%d.png" % split)
        plt.clf()

        return R2(y_true, y_pred), RMSE(y_true, y_pred), MSE(y_true, y_pred)

    elif mode.startswith('class'):
        test_results = candidate_model.predict(test_descs)

        acc = accuracy_score(test_acts, test_results)
        sens = recall_score(test_acts, test_results, pos_label=1)
        spec = recall_score(test_acts, test_results, pos_label=0)

        print('Accuracy: %5.3f' % acc)
        print('Sensitivity: %5.3f ' % sens)
        print('Specificity: %5.3f ' % spec)

        confmat = confusion_matrix(test_acts, test_results, labels=[1, 0])

        print('')
        print(confmat)
        print('')

        train_pos = 0
        train_neg = 0
        test_pos = 0
        test_neg = 0

        for i in train_acts:
            if i:
                train_pos += 1
            else:
                train_neg += 1

        for i in test_acts:
            if i:
                test_pos += 1
            else:
                test_neg += 1

        print('Training Positives: ', train_pos)
        print('Training Negatives: ', train_neg)
        print('Testing Positives:', test_pos)
        print('Testing Negatives: ', test_neg)

        sample = open(out_model_dir + '/prediction_class_%s_%d_%d.dat' % (ml_type, split, rand_state), 'w')
        print('Accuracy: %5.3f' % acc, file=sample)
        print('Sensitivity: %5.3f ' % sens, file=sample)
        print('Specificity: %5.3f ' % spec, file=sample)
        print('', file=sample)
        print(confmat, file=sample)
        print('', file=sample)
        print('Training Positives: ', train_pos, file=sample)
        print('Training Negatives: ', train_neg, file=sample)
        print('Testing Positives:', test_pos, file=sample)
        print('Testing Negatives: ', test_neg, file=sample)
        sample.close()

        return test_results, acc


def MSE(y_true, y_pred):
    mse = mean_squared_error(y_true, y_pred)
    return mse


def R2(y_true, y_pred):
    r2 = r2_score(y_true, y_pred)
    return r2


def RMSE(y_true, y_pred):
    rmse = sqrt(MSE(y_true, y_pred))
    return rmse


def read_mols(mode, method, basename, datadir='Default', modeldir='Default'):
    currworkdir = os.getcwd()
    if datadir == 'Default':
        datadir = os.path.join(currworkdir, 'data')
    else:
        if not os.path.isdir(datadir):
            print("error: ", datadir, " is not a directory. exiting.")
            exit(2)

    if modeldir == 'Default':
        modeldir = os.path.join(currworkdir, 'models')
    else:
        if not os.path.isdir(modeldir):
            print("error: ", modeldir, " is not a directory. exiting.")
            exit(2)
        else:
            print('setting modeldir to ', modeldir, '.')
            print('Have you set the random splits to be correct for the model?')

    mol_data_filename = basename + '.smi'
    act_data_filename = basename + '.act'
    moldatafile = os.path.join(datadir, mol_data_filename)
    actdatafile = os.path.join(datadir, act_data_filename)

    # output_ext = "%s_%s_%d_%d" % (mode, method, int(rand_split), int(rand_state))
    model_filename = "model_%s.dat" % output_ext
    index_filename = "indices_%s.dat" % output_ext
    appdom_fp_filename = "training-FPs_%s.dat" % output_ext
    appdom_rad_filename = "AD-radius_%s.dat" % output_ext

    if mode.startswith('class'):
        if os.path.isfile(actdatafile):
            actfh = open(actdatafile)

            activities = []  # array of tuples: (activity, molecule name)

            for actline in actfh:
                line = actline.split()
                act = float(line[1])
                actname = line[0]
                activities.append((act, actname))

            actfh.close()

    elif mode.startswith('reg') and method == 'xgb':

        bits_filename = "sigbits_%s.dat" % output_ext
        bits_file = os.path.join(modeldir, bits_filename)
        with open(bits_file, 'rb') as f:
            significant_bits = pickle.load(f)

    model_file = os.path.join(modeldir, model_filename)
    loaded_model = pickle.load(open(model_file, "rb"))

    index_file = os.path.join(modeldir, index_filename)
    with open(index_file, 'rb') as f:
        indexes = pickle.load(f)

    appdom_fp_file = os.path.join(modeldir, appdom_fp_filename)
    with open(appdom_fp_file, 'rb') as f:
        appdom_fps = pickle.load(f)

    appdom_rad_file = os.path.join(modeldir, appdom_rad_filename)
    with open(appdom_rad_file, 'rb') as f:
        appdom_radius = pickle.load(f)

    # Read in molecules from test set
    molfh = open(moldatafile)

    molecules = []  # array of tuples: (molecule, molecule name)

    for molline in molfh:
        line = molline.split()
        mol = Chem.MolFromSmiles(line[0])
        molname = line[1]
        molecules.append((mol, molname))

    molfh.close()

    mols_train = []
    molnames_train = []

    if 'activities' in locals():
        acts_train = []
        actnames_train = []

    for i in range(len(molecules)):
        mols_train.append(molecules[i][0])
        molnames_train.append(molecules[i][1])
        if mode.startswith('class') and 'activities' in locals():
            acts_train.append(activities[i][0])
            actnames_train.append(activities[i][1])

    # Standardize structures
    s = Standardizer()
    standard_mols_train = []
    for mol in mols_train:
        standard_mols_train.append(s.standardize(mol))

    return_dict = {}

    return_dict['molnames'] = molnames_train
    return_dict['molecules'] = standard_mols_train
    return_dict['model'] = loaded_model
    return_dict['inds'] = indexes
    if mode.startswith('reg') and method == 'xgb':
        return_dict['sigbits'] = significant_bits
    elif mode.startswith('class') and 'activities' in locals():
        return_dict['activities'] = acts_train
    return_dict['ad_fps'] = appdom_fps
    return_dict['ad_radius'] = appdom_radius

    return return_dict


def make_preds(*args, mode):

    if len(args) < 3 and len(args) > 4:
        print("Error: incorrect number of arguments passed to check_appdom()")

    mol_names = args[0]
    descs = args[1]
    saved_model = args[2]
    if len(args) == 4:
        if mode.startswith('reg'):
            split = args[3]
        elif mode.startswith('class'):
            y_true = args[3]

    # Make predictions for test data using previous model
    y_pred = saved_model.predict(descs)


    prediction_results = {}
    prediction_results['predictions'] = y_pred

    if 'y_true' in locals() and mode.startswith('class'):
        print("Molecule\tActual Act.\tPredicted Act.")
        for out_line in range(len(mol_names)):
            print(mol_names[out_line], "\t", y_true[out_line], "\t\t", y_pred[out_line])
        print("")
        print("Accuracy Score:",accuracy_score(y_true,y_pred))
        print("")
        confmat = confusion_matrix(y_true,y_pred, labels=[1,0])
        print(confmat)
        prediction_results['accuracy'] = accuracy_score(y_true, y_pred)
    else:
        print("Molecule\tPredicted Act.")
        for out_line in range(len(mol_names)):
            print(mol_names[out_line], "\t", y_pred[out_line])

    return prediction_results


def summarize_preds(*args):
    if len(args) != 2:
        print("Error: incorrect number of arguments passed to summarize_preds()")

    names = args[0]
    preds_list = args[1]

    nsplits = len(preds_list)
    npreds = len(preds_list[0])
    pivot_list = []
    for i in range(npreds):
        pivot_list.append([])

    for pred in range(npreds):
        for trial in range(nsplits):
            pivot_list[pred].append(preds_list[trial][pred])

    predmean_list = []
    predstd_list = []

    for i in range(npreds):
        predmean_list.append(np.mean(pivot_list[i]))
        predstd_list.append(np.std(pivot_list[i]))

    print("Compound\tPredicted\tStdDev")
    for i in range(len(names)):
        print("%s \t %2.3f \t %2.3f" % (names[i], predmean_list[i], predstd_list[i]))

    return names, predmean_list, predstd_list
