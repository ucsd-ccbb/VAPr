import unittest
import sys
import os
sys.path.append('/Users/carlomazzaferro/Documents/Code/variantannotation-master/variantannotation')
#quick and dirty way of importing functions


from variantannotation import csv_to_df

file_name = os.path.dirname(os.path.realpath('__file__')) + '/test_file.csv'
sample_list = csv_to_df.open_and_parse(file_name)
sample_df = csv_to_df.parse_to_df(sample_list)


class OpenParseTest(unittest.TestCase):

    def setUp(self):
        self.data = file_name

    def test_csv_read_data_headers(self):
        list2 = csv_to_df.open_and_parse(self.data)
        self.assertEqual(
            ['Chr', 'Start', 'End', 'Ref', 'Alt', 'Func.knownGene', 'Gene.knownGene', 'GeneDetail.knownGene',
             'ExonicFunc.knownGene', 'AAChange.knownGene', 'tfbsConsSites', 'cytoBand', 'targetScanS',
             'genomicSuperDups', 'gwasCatalog', 'esp6500siv2_all', '1000g2015aug_all', 'PopFreqMax', '1000G2012APR_ALL',
             '1000G2012APR_AFR', '1000G2012APR_AMR', '1000G2012APR_ASN', '1000G2012APR_EUR', 'ESP6500si_ALL',
             'ESP6500si_AA', 'ESP6500si_EA', 'CG46', 'clinvar_20140929', 'cosmic70', 'nci60', 'Otherinfo'],

            list2[0]
            )

    def test_parse_to_df(self):
        self.assertEqual(list(sample_df.columns), list(csv_to_df.parse_to_df(sample_list).columns))

    def test_parse_csv_chunks(self):
        list3 = csv_to_df.open_and_parse_chunks(self.data[0:1000], 1000, 0)
        list2 = csv_to_df.open_and_parse(self.data)

        #Assert that the first 1000 elements in list are the same for both processes
        self.assertEqual(list3[0:1000],list2[0:1000])


"""
    def test_csv_read_data_points(self):
        self.assertEqual(read_data(self.data)[1][7], '87')

    def test_get_min_score_difference(self):
        self.assertEqual(get_min_score_difference(self.parsed_data), 1)

    def test_get_team(self):
        index_value = get_min_score_difference(self.parsed_data)
        self.assertEqual(get_team(index_value), 'Liverpool')
"""

if __name__ == '__main__':
    unittest.main()