"""
Written by Andy Guan
This file contains functions to perform unit tests on our buildmodel module
"""
## define enviroment

import sys
from pathlib import Path
import unittest
import pickle
import os
import numpy as np

home = str(Path.home())
base_dir = home+'/repositories/QSAR_DAT-hERG'
core_dir = base_dir+'/core'
conf_dir = core_dir+'/conf'
unittest_data_dir = core_dir+'/unittest_data'
sys.path.insert(0, conf_dir)
sys.path.insert(0, core_dir)

from buildmodel import *

""" Define core files used to run unittest
1) in_file - contains the data used to build model during unit test
2) reference - directory with reference output files to ensure fidelity between test run and previous runs
3) negcon - shorter version of in_file, used as a negative condition for testing read_data function
"""

in_file = unittest_data_dir+"/data4buildmodels/pubdata_40"
reference = unittest_data_dir+"/reference"
negcon = unittest_data_dir+"/data4buildmodels/pubdata_20_ctrl"


class TestBuildModel(unittest.TestCase):
    def setUp(self):
        """ Function used to setup for unit testing """
        self.mode = 'reg'
        self.method = 'xgb'
        self.output_dir = "test_output"
        self.rand_split = gen_random_splits(control_seed=2020, num_splits=1)
        self.rand_states = gen_random_splits(control_seed=501, num_splits=1)

        # create a temporary output directory for testing
        if not os.path.exists(self.output_dir):
            Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    # remove all files and directories created during testing
    def tearDown(self):
        """ Function used to remove intermediate elements created from unit testing """
        if os.path.isdir(self.output_dir):
            os.system("rm -rf %s" % self.output_dir)

    # setup helper function for buildmodel
    def startUp(self):
        """ Setup helper function

        Performs the early steps in buildmodel
        Used to setup for testing of some functions
        Prevents redundant code lines in testing functions
        """
        mols, acts, deletes, changes = read_data4buildmodel(in_file, self.mode)
        mols, acts = curate_mols(mols, acts, deletes, changes)
        train_mols, train_names, train_acts = all_data(mols, acts)
        # output_ext is required to be defined to be able to read back results during testing
        output_ext = get_output_ext('a', 'b', 0, 1, 2)

        return train_mols, train_names, train_acts, output_ext

    def startUp1(self):
        """ Setup helper function for late-stage buildmodel functions """
        train_mols, train_names, train_acts, output_ext = self.startUp()
        train_topo_descs = calc_topo_descs(train_mols)
        train_topo_descs, topo_index, topo_names = prune_topo_descs(self.mode, train_topo_descs,
                                                                    train_acts, self.output_dir)
        train_phore_descs = calc_phore_descs(train_mols)
        train_phore_descs, phore_sigbits, phore_names = prune_phore_descs(train_phore_descs, self.output_dir)
        train_descs = np.concatenate((train_topo_descs, train_phore_descs), axis=1)

        return train_mols, train_acts, train_names, train_descs, output_ext, topo_index, phore_sigbits

    def startUp2(self):
        """ Setup helper function for prediction functions """
        input_data = read_mols(self.mode, self.method, "pred", datadir="core/unittest_data/data4buildmodels",
                               modeldir=reference)
        molnames = input_data['molnames']
        mols = input_data['molecules']
        model = input_data['model']
        inds = input_data['inds']
        sigbits = input_data['sigbits']
        ad_fps = input_data['ad_fps']
        ad_radius = input_data['ad_radius']
        appdom_results = check_appdom(ad_fps, ad_radius, mols, molnames, step="pred")
        mols = appdom_results['test_mols']
        molnames = appdom_results['test_names']
        molecules_rej = appdom_results['rej_mols']
        molnames_rej = appdom_results['rej_names']
        descriptors = calc_topo_descs(mols, inds)
        result_list = []
        phore_descriptors = calc_phore_descs(mols, sigbits)
        descriptors = np.concatenate((descriptors, phore_descriptors), axis=1)

        return molnames, descriptors, model, result_list

    # Check using negative control
    def test_read_data(self):
        """ Test read_data function """
        mols1, acts1, deletes1, changes1 = read_data4buildmodel(in_file, self.mode)
        mols2, acts2, deletes2, changes2 = read_data4buildmodel(negcon, self.mode)
        mols3, acts3, deletes3, changes3 = read_data4buildmodel(in_file, self.mode)
        mols4, acts4, deletes4, changes4 = read_data4buildmodel(negcon, self.mode)

        self.assertNotEqual(mols1, mols2)
        self.assertNotEqual(acts1, acts2)

        self.assertEqual(len(mols1), len(mols3))
        self.assertEqual(acts1, acts3)

    def test_curate_mols(self):
        """ Test curate_mols function

        Input arrays must have same length (from smile and activity data)
        """
        # define test/dummy arrays to test functionality of curate_mols
        a = [(1, 2), (3, 4)]
        b = [(5, 6), (7, 8)]
        c = [(1, 2)]

        with self.assertRaises(SystemExit):
            # test to make sure that function checks for equal length
            curate_mols(a, [], [], [])
            curate_mols(a, b, [], [])

        # we remove nothing, so the resulting data should just be untouched
        self.assertEqual(curate_mols(a, a, {}, {}), (a, a))

        # check that delete list works as intended
        # we will delete 4, so there should only be (1,2) in the output
        self.assertEqual(curate_mols(a, a, {4}, {}), (c, c))

        mols, acts, deletes, changes = read_data4buildmodel(in_file, 'reg')
        mols, acts = curate_mols(mols, acts, deletes, changes)

        # TO-DO: leave more notes

    def test_split_data(self):
        """ Test split_data function """
        mols, acts, deletes, changes = read_data4buildmodel(in_file, self.mode)
        mols, acts = curate_mols(mols, acts, deletes, changes)
        train_mols, train_names, train_acts, test_mols, test_names, test_acts = split_data(mols, acts, 0.10, 0)
        self.assertGreater(len(mols), len(train_mols))
        self.assertGreater(len(mols), len(test_mols))
        self.assertEqual(len(train_mols) + len(test_mols), len(mols))

    def test_all_data(self):
        """ Test all_data function """
        a = [[1, 2], [3, 4]]
        b = [1, 3]
        c = [2, 4]
        self.assertEquals(all_data(a, a), (b, c, b))

    def test_get_output_ext(self):
        """ Test get_output_ext function """
        self.assertEquals(get_output_ext('a', 'b', 0, 1, 2), "a_b_0.00_1_2")

    def test_get_output_dir(self):
        """ Test get_output_dir function """
        self.assertEquals(get_output_dir('a', 'b', 0), "a_b_0.00")

    def test_calc_appdom(self):
        """ Test calc_appdom function

        Uses a reference results file to ensure fidelity of results
        """
        train_mols, train_names, train_acts, output_ext = self.startUp()
        self.assertEqual(len(calc_appdom(train_mols, self.output_dir)), 2)

        ad_fps, ad_rad = calc_appdom(train_mols, self.output_dir)
        f = open(self.output_dir+('/AD-radius_%s.dat' % output_ext), 'rb')
        g = open(reference+'/AD-radius_ref.dat', 'rb')
        self.assertEquals(pickle.load(f), pickle.load(g))
        f.close()
        g.close()

        f = open(self.output_dir + ('/training-FPs_%s.dat' % output_ext), 'rb')
        g = open(reference + '/training-FPs_ref.dat', 'rb')
        self.assertEquals(pickle.load(f), pickle.load(g))
        f.close()
        g.close()

    def test_check_appdom(self):
        with self.assertRaises(SystemExit):
            check_appdom(1, step='b')
            check_appdom(1, 2, 3, 4, 5, 6, step='b')

    def test_calc_topo_descs(self):
        """ Test calc_topo_descs function """
        train_mols, train_names, train_acts, output_ext = self.startUp()
        train_topo_descs = calc_topo_descs(train_mols)
        sel = [1, 2, 3]
        # check that entries are properly removed when specified
        # first entry should have higher count than second
        self.assertGreater(train_topo_descs.shape[1], calc_topo_descs(train_mols, sel).shape[1])
        # check that the topo_descs have 200 entries (same as descriptor number)
        self.assertEqual(train_topo_descs.shape[1], 200)
        # check that the right number of entries are maintained when remove is specified
        self.assertEquals(calc_topo_descs(train_mols, sel).shape[1], len(sel))

    def test_prune_topo_descs(self):
        """ Test prune_topo_descs function

        Uses a reference results file to ensure fidelity of results
        """
        train_mols, train_names, train_acts, output_ext = self.startUp()
        train_topo_descs = calc_topo_descs(train_mols)
        self.assertEqual(len(prune_topo_descs(self.mode, train_topo_descs, train_acts, self.output_dir)), 3)
        a, b, c = prune_topo_descs(self.mode, train_topo_descs, train_acts, self.output_dir)
        f = open(self.output_dir+("/indices_%s.dat" % output_ext), 'rb')
        g = open(reference+"/indices_ref.dat", 'rb')
        self.assertEquals(pickle.load(f), pickle.load(g))
        f.close()
        g.close()

    def test_calc_phore_descs(self):
        """ Test calc_phore_descs function

        Uses a reference results file to ensure fidelity of results
        """
        train_mols, train_names, train_acts, output_ext = self.startUp()
        train_phore_descs = calc_phore_descs(train_mols)
        train_phore_descs, phore_sigbits, phore_names = prune_phore_descs(train_phore_descs, self.output_dir)
        a, b = calc_phore_descs(train_mols, phore_sigbits, testing=True)

        f = open(reference + "/regression_xgb_0.00.log", 'r')
        lines = f.readlines()
        self.assertEquals(lines[1].replace("\n", ""), a)
        self.assertEquals(lines[2].replace("\n", ""), b)
        f.close()

    def test_prune_phore_descs(self):
        """ Test prune_phore_descs function

        Uses a reference results file to ensure fidelity of results
        """
        train_mols, train_names, train_acts, output_ext = self.startUp()
        train_phore_descs = calc_phore_descs(train_mols)
        self.assertEqual(len(prune_phore_descs(train_phore_descs, self.output_dir)), 3)
        train_phore_descs, phore_sigbits, phore_names = prune_phore_descs(train_phore_descs, self.output_dir)
        f = open(self.output_dir+('/sigbits_%s.dat' % output_ext), 'rb')
        g = open(reference+'/sigbits_ref.dat', 'rb')
        self.assertEquals(pickle.load(f), pickle.load(g))
        f.close()
        g.close()

    def test_build_model(self):
        """ Test build_model function

        Uses a reference results file to ensure fidelity of results
        """
        train_mols, train_acts, train_names, train_descs, output_ext, topo_index, phore_sigbits = self.startUp1()
        model, model_score, best_params = build_model(self.mode, self.method, self.rand_states[0],
                                                      train_descs, train_acts, self.output_dir)

        # compare the string from model output to make sure parameters are the same
        f = open(self.output_dir + ('/model_%s.dat' % output_ext), 'rb')
        # f = open(reference + '/model_ref.dat', 'rb')
        # g = open(self.output_dir + ('/model_%s.dat' % output_ext), 'rb')
        g = open(reference + '/model_ref.dat', 'rb')
        self.assertEqual(str(pickle.load(f)), str(pickle.load(g)))
        f.close()
        g.close()

    def test_read_mols(self):
        """ Test read_mols function

        Uses a reference results file to ensure fidelity of results
        """
        a = np.load(reference+"/readmols4pred.npy", allow_pickle=True)
        input_data = read_mols(self.mode, self.method, "pred", datadir="core/unittest_data/data4buildmodels",
                               modeldir=reference)
        ref = dict(enumerate(a.flatten()))[0]
        self.assertEqual(ref["molnames"], input_data["molnames"])
        self.assertEqual(ref["inds"], input_data["inds"])
        self.assertEqual(ref["sigbits"], input_data["sigbits"])
        self.assertEqual(str(ref["model"]), str(input_data["model"]))

    def test_make_preds(self):
        """ Test make_preds function

        Uses a reference results file to ensure fidelity of results
        """
        a = np.load(reference+"/makepreds.npy", allow_pickle=True)
        ref = dict(enumerate(a.flatten()))[0]
        molnames, descriptors, model, result_list = self.startUp2()
        pred_results = make_preds(molnames, descriptors, model, self.rand_split[0], mode=self.mode)

        self.assertTrue((pred_results["predictions"] == ref["predictions"]).all())

    def test_summarize_preds(self):
        """ Test summarize_preds function

        Uses a reference results file to ensure fidelity of results
        """
        molnames, descriptors, model, result_list = self.startUp2()
        pred_results = make_preds(molnames, descriptors, model, self.rand_split[0], mode=self.mode)
        result_list.append(pred_results['predictions'])
        compound, pred_mean, pred_error = summarize_preds(molnames, result_list)
        data = {"compound": compound, "mean": pred_mean, "stdev": pred_error}

        a = np.load(reference+"/summarizepreds.npy", allow_pickle=True)
        ref = dict(enumerate(a.flatten()))[0]

        self.assertTrue(data["compound"] == ref["compound"])
        self.assertTrue(data["mean"] == ref["mean"])
        self.assertTrue(data["stdev"] == ref["stdev"])

    def test_predict_model(self):
        """ Test predict_model function

        Uses a reference results file to ensure fidelity of results
        """
        a = np.load(reference + "/predictmodel.npy", allow_pickle=True)
        ref = dict(enumerate(a.flatten()))[0]

        train_mols, train_acts, train_names, train_descs, output_ext, topo_index, phore_sigbits = self.startUp1()
        model, model_score, best_params = build_model(self.mode, self.method, self.rand_states[0],
                                                      train_descs, train_acts, self.output_dir)

        test_mols, test_acts, test_deletes, test_changes = read_data4buildmodel(negcon, self.mode)
        test_mols, test_acts = curate_mols(test_mols, test_acts, test_deletes, test_changes)
        test_mols, test_names, test_acts = all_data(test_mols, test_acts)

        test_topo_descs = calc_topo_descs(test_mols, topo_index)
        test_phore_descriptors = calc_phore_descs(test_mols, phore_sigbits)
        test_descs = np.concatenate((test_topo_descs, test_phore_descriptors), axis=1)

        test_r2, test_rmse, test_mse = predict_model(model, test_descs, test_acts, train_acts, self.rand_split[0],
                                                     self.output_dir, self.mode, self.method, self.rand_states[0])

        self.assertEqual(ref["r2"], test_r2)
        self.assertEqual(ref["rmse"], test_rmse)
        self.assertEqual(ref["mse"], test_mse)




