import unittest
import sys
import os
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
        self.chunksize = 1000
        self.hgvs = parser_models.HgvsParser(self.vcf_file)
        self.csv_parsing = parser_models.TxtParser(self.txt_file)

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
             '1000g20XX',
             'esp6500siv2_all',
             'cosmic70',
             'nci60',
             'otherinfo',
             'genotype']

        for dictionary in dat:
            for key in list(dictionary.keys()):
                self.assertTrue(key in l)

    def test_hgvs_parsing(self):
        return 'a'
        ls = []
        while self.hgvs.num_lines > self.step * self.chunksize:
            ls.extend(self.hgvs.get_variants_from_vcf(self.step))
            self.step += 1

        self.assertEqual(len(ls), self.txt_len)

    def test_my_var_query(self):
        #open_file = parser_models.VariantParsing()
        pass

    def test_parse_csv_chunks(self):
        list3 = 3 #csv_to_df.open_and_parse_chunks(self.data, 1000, 0)
        list2 = 4 #csv_to_df.open_and_parse(self.data)
        self.assertEqual(len(list3[0:1000]), len(list2[0:1000]))

    def push_myvar_to_db_test(self):

        chunk_list = []
        while self.hgvs.num_lines > self.step * self.chunksize:

            list_hgvs_ids = self.hgvs.get_variants_from_vcf(self.step)
            myvariants_variants = self.get_dict_myvariant(list_hgvs_ids)

            chunk_list.append(myvariants_variants)
            self.step += 1

        return chunk_list

    def push_to_db(self):

        chunk_list = []
        while self.csv_parsing.num_lines > self.step*self.chunksize:

            list_hgvs_ids = self.hgvs.get_variants_from_vcf(self.step)
            myvariants_variants = self.get_dict_myvariant(list_hgvs_ids)
            csv_variants = self.csv_parsing.open_and_parse_chunks(self.step)

            merged_list = []
            for i, _ in enumerate(myvariants_variants):
                merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

            chunk_list.append(merged_list)

    @staticmethod
    def merge_dict_lists(*dict_args):

        result = {}
        for dictionary in dict_args:
            result.update(dictionary)
        return result

    def get_dict_myvariant(self, variant_list):


        mv = myvariant.MyVariantInfo()
        # This will retrieve a list of dictionaries
        variant_data = mv.getvariants(variant_list, as_dataframe=False)
        variant_data = self.remove_id_key(variant_data)
        return variant_data

    @staticmethod
    def remove_id_key(variant_data):

        for dic in variant_data:
            dic['hgvs_id'] = dic.pop("_id", None)
            dic['hgvs_id'] = dic.pop("query", None)

        return variant_data


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