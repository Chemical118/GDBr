from unittest import TestCase
from gdbr.annotate import annotate_main

import os
import csv
import shutil


class HomologyTest(TestCase):
    def test_main_example_run(self):
        annotate_main('test/example/ref.fa', ['test/example/qry.fa'], ['test/example/example_sv.csv'], workdir='data_', pbar=False)

        # clean test
        shutil.rmtree('data_')
        shutil.rmtree('dsbr')


