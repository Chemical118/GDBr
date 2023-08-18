from unittest import TestCase
from gdbr.correct import correct_main, get_correctrd_location_by_idx

import shutil

class PreprocessTest(TestCase):
    def test_main_example(self):
        correct_main('test/example/ref.fa', ['test/example/qry.fa'], ['test/example/variants.vcf'], 'human', min_sv_size=50)

        # clean test
        shutil.rmtree('sv')
        shutil.rmtree('data')

    def test_corrected_index(self):
        self.assertEqual((0, 9, 10, 0, 4, 5), get_correctrd_location_by_idx("ATGCTATGCT", "GCACA", 0, 9, 0, 4))
        self.assertEqual((2, 9, 8, 2 ,4, 3), get_correctrd_location_by_idx("GCGCTATGCT", "GCACA", 0, 9, 0, 4))
        self.assertEqual((0, 7, 8, 0 ,2, 3), get_correctrd_location_by_idx("ATGCTATGCA", "GCACA", 0, 9, 0, 4))
        self.assertEqual((3, 7, 5, 3 ,2, 0), get_correctrd_location_by_idx("GCATTATGCA", "GCACA", 0, 9, 0, 4))