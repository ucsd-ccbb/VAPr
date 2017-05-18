import unittest
import sys
import os
import csv
import myvariant
#quick and dirty way of importing functions
from VAPr import parser_models

sys.path.append('/Users/carlomazzaferro/Documents/Code/variant-annotation/variantannotation')
vcf_file = os.path.dirname(os.path.realpath('__file__')) + '/Normal_targeted_seq.vcf'

txt_file = os.path.dirname(os.path.realpath('__file__')) + '/annotated.hg19_multianno.txt'


class OpenParseTest(unittest.TestCase):

    def setUp(self):
        self.txt_file = txt_file
        self.vcf_file = vcf_file
        self.step = 0
        self.txt_len = sum(1 for _ in open(self.txt_file))
        self.chunksize = 950
        self.hgvs = parser_models.HgvsParser(self.vcf_file)
        self.csv_parsing = parser_models.TxtParser(self.txt_file)
#locals()

    def test_csv_read_data_headers(self):
        dat = self.csv_parsing.open_and_parse_chunks(self.step)
        l = ['chr',
             'start',
             'end',
             'ref',
             'alt',
             'func_knowngene',
             'gene_knowngene',
             'genedetail_knowngene',
             'exonicfunc_knowngene',
             'tfbsconssites',
             'cytoband',
             'genomicsuperdups',
             '1000g2015aug_all',
             'esp6500siv2_all',
             'cosmic70',
             'nci60',
             'otherinfo',
             'genotype']

        for dictionary in dat:
            for key in list(dictionary.keys()):
                self.assertTrue(key in l)

    """
    def test_hgvs_parsing(self):

        step = 0
        chunksize = 600
        list_hgvs_ids = []

        while self.txt_len > step * chunksize:
            list_hgvs_ids.extend(self.hgvs.get_variants_from_vcf(step))
            step += 1

        self.assertEqual(len(list_hgvs_ids), self.txt_len - 1) #to account for the header



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