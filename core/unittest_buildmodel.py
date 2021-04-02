## define enviroment
import sys, os
from pathlib import Path
home = str(Path.home())
base_dir = home+'/repositories/herg/hERGvDAT/'
core_dir = base_dir+'/core'
conf_dir = core_dir+'/conf'
unittest_data_dir = core_dir+'/unittest_data'
sys.path.insert(0, conf_dir)
sys.path.insert(0, core_dir)


import unittest
import os
import io
import shutil
from pathlib import Path
from buildmodel import *
from descriptor_setup import dnames, dlist
import pickle


in_file = unittest_data_dir+"/data4buildmodels/pubdata_40"
reference = unittest_data_dir+"/reference"



class TestBuildModel(unittest.TestCase):
    def setUp(self):
        self.mode = 'reg'
        self.output_dir = "test_output"
        self.rand_split = gen_random_splits(control_seed=2020, num_splits=1)
        self.rand_states = gen_random_splits(control_seed=501, num_splits=1)
        if not os.path.exists(self.output_dir):
            Path(self.output_dir).mkdir(parents=True, exist_ok=True)

    def tearDown(self):
        if os.path.isdir(self.output_dir):
            shutil.rmtree(self.output_dir)

    # setup helper function
    def startUp(self):
        mols, acts, deletes, changes = read_data4buildmodel(in_file, self.mode)
        mols, acts = curate_mols(mols, acts, deletes, changes)
        train_mols, train_names, train_acts = all_data(mols, acts)
        output_ext = get_output_ext('a', 'b', 0, 1, 2)

        return train_mols, train_names, train_acts, output_ext

    # Still thinking about how to do this one
    def test_read_data(self):
        a = read_data4buildmodel(in_file, self.mode)

    def test_curate_mols(self):
        a = [(1, 2), (3, 4)]
        b = [(5, 6), (7, 8)]
        c = [(1, 2)]
        with self.assertRaises(SystemExit):
            curate_mols(a, [], [], [])
            curate_mols(a, b, [], [])

        self.assertEqual(curate_mols(a, a, {}, {}), (a, a))
        self.assertEqual(curate_mols(a, a, {4}, {}), (c, c))

        mols, acts, deletes, changes = read_data4buildmodel(in_file, 'reg')
        mols, acts = curate_mols(mols, acts, deletes, changes)

    def test_split_data(self):
        mols, acts, deletes, changes = read_data4buildmodel(in_file, self.mode)
        mols, acts = curate_mols(mols, acts, deletes, changes)
        train_mols, train_names, train_acts, test_mols, test_names, test_acts = split_data(mols, acts, 0.10, 0)
        self.assertGreater(len(mols), len(train_mols))
        self.assertGreater(len(mols), len(test_mols))
        self.assertEqual(len(train_mols) + len(test_mols), len(mols))

    def test_all_data(self):
        a = [[1, 2], [3, 4]]
        b = [1, 3]
        c = [2, 4]
        self.assertEquals(all_data(a, a), (b, c, b))

    def test_get_output_ext(self):
        self.assertEquals(get_output_ext('a', 'b', 0, 1, 2), "a_b_0.00_1_2")

    def test_get_output_dir(self):
        self.assertEquals(get_output_dir('a', 'b', 0), "a_b_0.00")

    def test_calc_appdom(self):
        train_mols, train_names, train_acts, output_ext = self.startUp()
        self.assertEqual(len(calc_appdom(train_mols, self.output_dir)), 2)

        ad_fps, ad_rad = calc_appdom(train_mols, self.output_dir)
        f = open(self.output_dir+('/AD-radius_%s.dat' % output_ext), 'rb')
        g = open(reference+'/AD-radius_ref.dat', 'rb')
        self.assertEquals(pickle.load(f), pickle.load(g))
        f.close()
        g.close()

    def test_check_appdom(self):
        with self.assertRaises(SystemExit):
            check_appdom(1, step='b')
            check_appdom(1, 2, 3, 4, 5, 6, step='b')

    def test_calc_topo_descs(self):
        train_mols, train_names, train_acts, output_ext = self.startUp()
        train_topo_descs = calc_topo_descs(train_mols)
        sel = [1, 2, 3]
        self.assertGreater(train_topo_descs.shape[1], calc_topo_descs(train_mols, sel).shape[1])
        self.assertEqual(train_topo_descs.shape[1], 200)
        self.assertEquals(calc_topo_descs(train_mols, sel).shape[1], len(sel))

    def test_prune_topo_descs(self):
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
        pass

    def test_prune_phore_descs(self):
        train_mols, train_names, train_acts, output_ext = self.startUp()
        train_phore_descs = calc_phore_descs(train_mols)
        self.assertEqual(len(prune_phore_descs(train_phore_descs, self.output_dir)), 3)
        train_phore_descs, phore_sigbits, phore_names = prune_phore_descs(train_phore_descs, self.output_dir)
        f = open(self.output_dir+('/sigbits_%s.dat' % output_ext), 'rb')
        g = open(reference+'/sigbits_ref.dat', 'rb')
        self.assertEquals(pickle.load(f), pickle.load(g))
        f.close()
        g.close()
