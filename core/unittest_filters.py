import unittest
from filters import *

# Initialize global variables
buffer = None
chembl_tsv_file = 'unittest/data/chembl25_DAT.tsv'
chembl_tsv_file2 = 'unittest/data/chembl25_hERG.tsv'
ref = [13273, 8392, 8369, 3418, 2659, 2368, 1616, 1561, 1559, 1559, 1540, 1540, 1540, 1189, 887, 798, 89]
ref2 = [20695, 20056, 19330, 10454, 8957, 6685, 2021, 1968, 1542, 1542, 1542, 1519, 1519, 1519, 1350, 737, 263, 474]


class TestFilters(unittest.TestCase):

    def test_01_read_data(self):
        global buffer
        buffer = read_data(chembl_tsv_file)
        self.assertEqual(len(buffer), ref[0])

    def test_02_confidence_filter(self):
        global buffer
        buffer = filter_confidence(buffer)
        self.assertEqual(len(buffer), ref[1])

    def test_03_assay_filter(self):
        global buffer
        buffer = filter_assay_type(buffer)
        self.assertEqual(len(buffer), ref[2])

    def test_04_affinity_filter(self):
        global buffer
        buffer = filter_affinity(buffer, keepIC50=False, keepKi=True)
        self.assertEqual(len(buffer), ref[3])

    def test_05_units_filter(self):
        global buffer
        buffer = filter_units(buffer)
        self.assertEqual(len(buffer), ref[4])

    def test_06_exact_filter(self):
        global buffer
        buffer = filter_exact(buffer)
        self.assertEqual(len(buffer), ref[5])

    # TODO: Find a way to specify DAT vs hERG test and inhibitor vs uptake test
    def test_07_assay_definition_filter(self):
        global buffer
        buffer, unused = filter_assaydefinition(buffer, 'DAT', 'inhibitor')
        self.assertEqual(len(buffer), ref[6])

    # TODO: Include test for filter secondary test set

    def test_08_small_set_filter(self):
        global buffer
        buffer = filter_small_sets(buffer, threshold=4)
        self.assertEqual(len(buffer), ref[7])

    def test_09_salts_filter(self):
        global buffer
        buffer = filter_salts(buffer)
        self.assertEqual(len(buffer), ref[8])

    def test_10_elements_filter(self):
        global buffer
        buffer = filter_elements(buffer)
        self.assertEqual(len(buffer), ref[9])

    def test_11_size_filter(self):
        global buffer
        buffer = filter_size(buffer)
        self.assertEqual(len(buffer), ref[10])

    def test_12_pchembl_filter(self):
        global buffer
        buffer = filter_pchembl_values(buffer, replace=True)
        self.assertEqual(len(buffer), ref[11])

    def test_13_weirdos_filter(self):
        global buffer
        buffer = filter_weirdos(buffer)
        self.assertEqual(len(buffer), ref[12])

    def test_14_deduplicate_filter(self):
        global buffer
        buffer = deduplicate_mols(buffer)
        self.assertEqual(len(buffer), ref[13])

    # TODO: Maybe also check the write functions but I'm not sure what a good way is yet


class TestFilters2(unittest.TestCase):

    def test_01_read_data(self):
        global buffer
        buffer = read_data(chembl_tsv_file2)
        self.assertEqual(len(buffer), ref2[0])

    def test_02_confidence_filter(self):
        global buffer
        buffer = filter_confidence(buffer)
        self.assertEqual(len(buffer), ref2[1])

    def test_03_assay_filter(self):
        global buffer
        buffer = filter_assay_type(buffer)
        self.assertEqual(len(buffer), ref2[2])

    def test_04_affinity_filter(self):
        global buffer
        buffer = filter_affinity(buffer, keepIC50=True, keepKi=False)
        self.assertEqual(len(buffer), ref2[3])

    def test_05_units_filter(self):
        global buffer
        buffer = filter_units(buffer)
        self.assertEqual(len(buffer), ref2[4])

    def test_06_exact_filter(self):
        global buffer
        buffer = filter_exact(buffer)
        self.assertEqual(len(buffer), ref2[5])

    def test_07_assay_definition_filter(self):
        global buffer
        buffer, unused = filter_assaydefinition(buffer, 'hERG', 'clamp')
        # buffer, unused = filter_assaydefinition(buffer, 'hERG', 'binding')
        self.assertEqual(len(buffer), ref2[6])

    def test_08_secondary_set_filter(self):
        global buffer
        buffer = filter_secondary_test_set(buffer)
        self.assertEqual(len(buffer), ref2[7])

    def test_09_small_set_filter(self):
        global buffer
        buffer = filter_small_sets(buffer, threshold=4)
        self.assertEqual(len(buffer), ref2[8])

    def test_10_salts_filter(self):
        global buffer
        buffer = filter_salts(buffer)
        self.assertEqual(len(buffer), ref2[9])

    def test_11_elements_filter(self):
        global buffer
        buffer = filter_elements(buffer)
        self.assertEqual(len(buffer), ref2[10])

    def test_12_size_filter(self):
        global buffer
        buffer = filter_size(buffer)
        self.assertEqual(len(buffer), ref2[11])

    def test_13_pchembl_filter(self):
        global buffer
        buffer = filter_pchembl_values(buffer, replace=True)
        self.assertEqual(len(buffer), ref2[12])

    def test_14_weirdos_filter(self):
        global buffer
        buffer = filter_weirdos(buffer)
        self.assertEqual(len(buffer), ref2[13])

    def test_15_deduplicate_filter(self):
        global buffer
        buffer = deduplicate_mols(buffer)
        self.assertEqual(len(buffer), ref2[14])


if __name__ == '__main__':
    unittest.defaultTestLoader.sortTestMethodsUsing = None
    unittest.main(testLoader=unittest.defaultTestLoader)
