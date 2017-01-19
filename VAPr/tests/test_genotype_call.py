import unittest
from variantannotation import genotype_calling

class TestGenCall(unittest.TestCase):

    def setUp(self):
        self.gen_call = [
                          ['GT:GQ:PL', '412434', '1/1:6:69,6,0'],
                          ['GT:AD:DP:GQ:PL', 'e21r3dc/', '0/1:10,3:13:87:87,0,379'],
                          ['GT:AD:DP:GQ:PL', 'r3dedwe', '0/1:97,97:194:99:3469,0,3458'],
                          ['GT:AD:DP:GQ:PL', 'frwvervef', '1/1:0,267:267:99:11155,803,0'],
                          ['3rfwewev', 'GT:AD:DP:GQ:PL', '1/2:23,48,156:227:99:3983,2789,3888,475,0,411'],
                          ['0000', 'GT:AD:DP:GQ:PL', '1/2:29,158,85:272:99:6336,1149,2024,4122,0,6055']
                        ]

        self.gen_call_split = [
                                {'GQ': 6,
                                 'GT': '1/1',
                                 'PL': [69.0, 6.0, 0.0]},

                                {'AD': [10.0, 3.0],
                                 'DP': 13,
                                 'GQ': 87,
                                 'GT': '0/1',
                                 'PL': [87.0, 0.0, 379.0]},

                                {'AD': [97.0, 97.0],
                                 'DP': 194,
                                 'GQ': 99,
                                 'GT': '0/1',
                                 'PL': [3469.0, 0.0, 3458.0]},

                                {'AD': [0.0, 267.0],
                                 'DP': 267,
                                 'GQ': 99,
                                 'GT': '1/1',
                                 'PL': [11155.0, 803.0, 0.0]},

                                {'AD': [23.0, 48.0, 156.0],
                                 'DP': 227,
                                 'GQ': 99,
                                 'GT': '1/2',
                                 'PL': [3983.0, 2789.0, 3888.0, 475.0, 0.0, 411.0]},

                                {'AD': [29.0, 158.0, 85.0],
                                 'DP': 272,
                                 'GQ': 99,
                                 'GT': '1/2',
                                 'PL': [6336.0, 1149.0, 2024.0, 4122.0, 0.0, 6055.0]}
                                ]


    def test_get_gen_call(self):

        dict_info = []
        for i in self.gen_call:
            dict_info.append(genotype_calling.return_dict(i))

        self.assertEqual(dict_info, self.gen_call_split)



if __name__ == '__main__':
    unittest.main()
