from unittest import TestCase
from gdbr.correct import correct_main, get_correctrd_location_by_idx

import shutil


test_ans = [[0, 'FND_IDX:(0, 1)'],
        [1, 'DEL', 'chr12:1-2500000', 95831, 96147, 95496, 95495],
        [2, 'SUB', 'chr12:1-2500000', 96424, 96783, 96129, 96131],
        [3, 'DEL', 'chr12:1-2500000', 123974, 124029, 123332, 123331],
        [4, 'INS', 'chr12:1-2500000', 144004, 144003, 143298, 143454],
        [5, 'FND_IDX:(1, 0)'],
        [6, 'FND_IDX:(1, 0)'],
        [7, 'FND_IDX:(1, 0)'],
        [8, 'FND_IDX:(0, 1)'],
        [9, 'FND_IDX:(0, 1)'],
        [10, 'SUB', 'chr12:1-2500000', 173356, 173365, 170488, 170554],
        [11, 'DEL', 'chr12:1-2500000', 212941, 213035, 209942, 209941],
        [12, 'SUB', 'chr12:1-2500000', 236113, 236114, 233018, 233217],
        [13, 'INS', 'chr12:1-2500000', 417102, 417101, 414280, 414719],
        [14, 'DEL', 'chr12:1-2500000', 448795, 448864, 446392, 446391],
        [15, 'DEL', 'chr12:1-2500000', 503357, 503428, 500850, 500849],
        [16, 'INS', 'chr12:1-2500000', 516232, 516231, 513477, 513568],
        [17, 'DEL', 'chr12:1-2500000', 525358, 525677, 522694, 522693],
        [18, 'DEL', 'chr12:1-2500000', 540040, 540081, 537048, 537047],
        [19, 'INS', 'chr12:1-2500000', 541866, 541865, 538832, 538965],
        [20, 'DEL', 'chr12:1-2500000', 573393, 573478, 570510, 570509],
        [21, 'DEL', 'chr12:1-2500000', 578352, 578435, 575374, 575373],
        [22, 'INS', 'chr12:1-2500000', 586645, 586644, 583583, 583768],
        [23, 'INS', 'chr12:1-2500000', 779627, 779626, 776785, 776923],
        [24, 'DEL', 'chr12:1-2500000', 915378, 916183, 912655, 912654],
        [25, 'INS', 'chr12:1-2500000', 1047858, 1047857, 1044362, 1044449],
        [26, 'INS', 'chr12:1-2500000', 1064042, 1064041, 1060632, 1067596],
        [27, 'DEL', 'chr12:1-2500000', 1075405, 1075545, 1078968, 1078967],
        [28, 'DEL', 'chr12:1-2500000', 1115565, 1116496, 1115123, 1115122],
        [29, 'FND_IDX:(0, 1)'],
        [30, 'FND_IDX:(0, 1)'],
        [31, 'DEL', 'chr12:1-2500000', 1148492, 1148545, 1147180, 1147179],
        [32, 'SUB', 'chr12:1-2500000', 1222653, 1222669, 1221229, 1221317],
        [33, 'DEL', 'chr12:1-2500000', 1287697, 1287872, 1286301, 1286300],
        [34, 'DEL', 'chr12:1-2500000', 1323134, 1323753, 1321587, 1321586],
        [35, 'INS', 'chr12:1-2500000', 1530544, 1530543, 1528341, 1528602],
        [36, 'SUB', 'chr12:1-2500000', 1537658, 1537659, 1535714, 1536138],
        [37, 'DEL', 'chr12:1-2500000', 1569026, 1569456, 1567515, 1567514],
        [38, 'INS', 'chr12:1-2500000', 1657402, 1657401, 1655438, 1655492],
        [39, 'FND_IDX:(0, 1)'],
        [40, 'FND_IDX:(0, 1)'],
        [41, 'INS', 'chr12:1-2500000', 1751468, 1751467, 1750185, 1750499],
        [42, 'INS', 'chr12:1-2500000', 1771135, 1771134, 1770174, 1770306],
        [43, 'DEL', 'chr12:1-2500000', 1904818, 1905090, 1903985, 1903984],
        [44, 'DEL', 'chr12:1-2500000', 1917808, 1917995, 1916700, 1916699],
        [45, 'INS', 'chr12:1-2500000', 2101863, 2101862, 2100532, 2100903],
        [46, 'INS', 'chr12:1-2500000', 2102205, 2102204, 2100874, 2101245],
        [47, 'SUB', 'chr12:1-2500000', 2103112, 2103113, 2101779, 2102157],
        [48, 'INS', 'chr12:1-2500000', 2210541, 2210540, 2209560, 2210884],
        [49, 'FND_IDX:(147, 148)'],
        [50, 'FND_IDX:(148, 1)'],
        [51, 'INS', 'chr12:1-2500000', 2353810, 2353809, 2354444, 2354576],
        [52, 'INS', 'chr12:1-2500000', 2408476, 2408475, 2409260, 2409774]]

class PreprocessTest(TestCase):
    def test_main_example(self):
        sv_data = correct_main('test/example/ref.fa', ['test/example/qry.fa'], ['test/example/variants.vcf'], file=False, pbar=False)
        for ans, func_ans in zip(test_ans, sv_data[0]):
            self.assertEqual(ans, func_ans)

        # clean test
        shutil.rmtree('data')

    def test_corrected_index(self):
        self.assertEqual((0, 9, 10, 0, 4, 5), get_correctrd_location_by_idx("ATGCTATGCT", "GCACA", 0, 9, 0, 4))
        self.assertEqual((2, 9, 8, 2 ,4, 3), get_correctrd_location_by_idx("GCGCTATGCT", "GCACA", 0, 9, 0, 4))
        self.assertEqual((0, 7, 8, 0 ,2, 3), get_correctrd_location_by_idx("ATGCTATGCA", "GCACA", 0, 9, 0, 4))
        self.assertEqual((3, 7, 5, 3 ,2, 0), get_correctrd_location_by_idx("GCATTATGCA", "GCACA", 0, 9, 0, 4))