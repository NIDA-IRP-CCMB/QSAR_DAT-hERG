import unittest
import os
import io
import shutil
from pathlib import Path
from buildmodel import *
# from django.test.util import tag
from descriptor_setup import dnames, dlist

in_file = "./unittest/buildmodel2/dataset/pubdata_40"
reference = "./unittest/reference"


class TestBuildModel(unittest.TestCase):
    def setUp(self):
        # method = getattr(self, self._testMethodName)
        # tags = getattr(method, 'tags', {})
        # if 'skip_setup' in tags:
        #     return
        self.mode = 'reg'
        self.output_dir = "test_output"
        self.rand_split = gen_random_splits(control_seed=2020, num_splits=1)
        self.rand_states = gen_random_splits(control_seed=501, num_splits=1)
        if not os.path.exists(self.output_dir):
            Path(self.output_dir).mkdir(parents=True, exist_ok=True)

        # mols, acts, deletes, changes = read_data4buildmodel(in_file, self.mode)
        # mols, acts = curate_mols(mols, acts, deletes, changes)
        # self.train_mols, self.train_names, self.train_acts = all_data(mols, acts)
        # self.output_ext = get_output_ext('a', 'b', 0, 1, 2)

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
        pass

    # @tag('skip_setup')
    def test_curate_mols(self):
        a = [(1, 2), (3, 4)]
        b = [(5, 6), (7, 8)]
        c = [(1, 2)]
        with self.assertRaises(SystemExit):
            curate_mols(a, [], [], [])
            curate_mols(a, b, [], [])

        self.assertEquals(curate_mols(a, a, {}, {}), (a, a))
        self.assertEquals(curate_mols(a, a, {4}, {}), (c, c))

        mols, acts, deletes, changes = read_data4buildmodel(in_file, 'reg')
        mols, acts = curate_mols(mols, acts, deletes, changes)
        print(mols)

    # @tag('skip_setup')
    def test_split_data(self):
        pass

    # @tag('skip_setup')
    def test_all_data(self):
        a = [[1, 2], [3, 4]]
        b = [1, 3]
        c = [2, 4]
        self.assertEquals(all_data(a, a), (b, c, b))

    # @tag('skip_setup')
    def test_get_output_ext(self):
        self.assertEquals(get_output_ext('a', 'b', 0, 1, 2), "a_b_0.00_1_2")

    # @tag('skip_setup')
    def test_get_output_dir(self):
        self.assertEquals(get_output_dir('a', 'b', 0), "a_b_0.00")

    def test_calc_appdom(self):
        # mols, acts, deletes, changes = read_data4buildmodel(in_file, self.mode)
        # mols, acts = curate_mols(mols, acts, deletes, changes)
        # train_mols, train_names, train_acts = all_data(mols, acts)
        train_mols, train_names, train_acts, output_ext = self.startUp()
        print(output_ext)

        self.assertEquals(len(calc_appdom(train_mols, self.output_dir)), 2)

        ad_fps, ad_rad = calc_appdom(train_mols, self.output_dir)

        # self.assertEquals(list(io.open(self.output_dir+('/AD-radius_%s.dat' % output_ext))),
        #                   list(io.open(reference+'./AD-radius_ref.dat')))

    def test_check_appdom(self):
        with self.assertRaises(SystemExit):
            check_appdom(1, step='b')
            check_appdom(1, 2, 3, 4, 5, 6, step='b')

    def test_calc_topo_descs(self):
        train_topo_descs = calc_topo_descs(train_mols)
        pass

    # def test_prune_topo_descs(self):
    #     train_mols, train_names, train_acts, output_ext = self.startUp()
    #     train_topo_descs = calc_topo_descs(train_mols)
    #     self.assertEquals(len(prune_topo_descs(self.mode, train_topo_descs, train_acts, self.output_dir)), 3)
    #     a, b, c = prune_topo_descs(self.mode, train_topo_descs, train_acts, self.output_dir)
    #
    #     self.assertEquals(list(io.open(self.output_dir+("/indices_%s.dat" % output_ext))),
    #                       list(io.open(reference+"/indices_ref.dat")))

    def test_calc_phore_descs(self):
        pass

    # def test_prune_phore_descs(self):
    #     train_mols, train_names, train_acts, output_ext = self.startUp()
    #     train_phore_descs = calc_phore_descs(train_mols)
    #     self.assertEquals(len(prune_phore_descs(train_phore_descs, self.output_dir)), 3)
    #     train_phore_descs, phore_sigbits, phore_names = prune_phore_descs(train_phore_descs, self.output_dir)
    #     self.assertEquals(list(io.open(self.output_dir+('/sigbits_%s.dat' % output_ext))),
    #                       list(io.open(reference+'/sigbits_ref.dat')))


