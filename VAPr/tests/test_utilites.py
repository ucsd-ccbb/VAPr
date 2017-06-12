import unittest
import sys
import os
sys.path.append("/Users/carlomazzaferro/Documents/Code/variantannotation-master/variantannotation")
from VAPr import models as parser_models


file_name = os.path.dirname(os.path.realpath('__file__')) + '/test_file.csv'

"""
sample_list = csv_to_df.open_and_parse(file_name)
sample_df = csv_to_df.parse_to_df(sample_list)

larger_list = csv_to_df.open_and_parse(file_name)
larger_df = csv_to_df.parse_to_df(larger_list)
"""



class TestUtilities(unittest.TestCase):

    def setUp(self):
        self.data = file_name
        self.cytoBand = ['1p36.33', '16p11.1', '16q11.2', '16q21', 'Xp22.32', 'Xp22.2', 'Xp22.11',
                         'Xq12', 'Xq13.1']

        self.cytoBand_split = [
                                ['1', 'p', '36', '.', '33'],
                                ['16', 'p', '11', '.', '1'],
                                ['16', 'q', '11', '.', '2'],
                                ['16', 'q', '21'],
                                ['X', 'p', '22', '.', '32'],
                                ['X', 'p', '22', '.', '2'],
                                ['X', 'p', '22', '.', '11'],
                                ['X', 'q', '12'],
                                ['X', 'q', '13', '.', '1']
                              ]

        self.cytoBand_dict = [
                                {'Band': 'p', 'Name': '1p36.33', 'Chromosome': '1', 'Region': 36, 'Sub_Band': 33},
                                {'Band': 'p', 'Name': '16p11.1', 'Chromosome': '16', 'Region': 11, 'Sub_Band': 1},
                                {'Band': 'q', 'Name': '16q11.2', 'Chromosome': '16', 'Region': 11, 'Sub_Band': 2},
                                {'Band': 'q', 'Name': '16q21', 'Chromosome': '16', 'Region': 21},
                                {'Band': 'p', 'Name': 'Xp22.32', 'Chromosome': 'X', 'Region': 22, 'Sub_Band': 32},
                                {'Band': 'p', 'Name': 'Xp22.2', 'Chromosome': 'X', 'Region': 22, 'Sub_Band': 2},
                                {'Band': 'p', 'Name': 'Xp22.11', 'Chromosome': 'X', 'Region': 22, 'Sub_Band': 11},
                                {'Band': 'q', 'Name': 'Xq12', 'Chromosome': 'X', 'Region': 12},
                                {'Band': 'q', 'Name': 'Xq13.1', 'Chromosome': 'X', 'Region': 13, 'Sub_Band': 1}
                             ]


        self.ids_to_expand = ['chrMT:g.16225TinsC,TCC', 'chr16:g.16225GTTdelinsCGTT,T', 'chr2:g.16225CdelinsCCT,CT',
                              'chr14:g.16225T>C', 'chrMT:g.16225TdelinsC,T']

        self.expanded_ids = ['chrMT:g.16225TinsC', 'chrMT:g.16225TinsTCC', 'chr16:g.16225GTTdelinsCGTT',
                             'chr16:g.16225GTTdelinsT', 'chr2:g.16225CdelinsCCT', 'chr2:g.16225CdelinsCT',
                             'chr14:g.16225T>C', 'chrMT:g.16225TdelinsC', 'chrMT:g.16225TdelinsT']

        self.input_names = ['/data/ccbb_internal/interns/Carlo/test_vcf/Tumor_RNAseq_rare_variants_VCF.vcf',
                            '/data/ccbb_internal/interns/Carlo/test_vcf/Tumor_targeted_seq.vqsr.vcf',
                            '/data/ccbb_internal/interns/Carlo/test_vcf/normal_targeted_seq.vcf',
                            '/data/ccbb_internal/interns/Carlo/test_vcf/normal_blood_WGS.vqsr.vcf',
                            's/data/ccbb_internal/interns/Carlo/test_vcf/somatic_mutect_old.vcf']


        self.output_path = '/data/ccbb_internal/interns/Carlo/test_vcf_out/'

        self.expected_outfile_names = ['/data/ccbb_internal/interns/Carlo/test_vcf_out/Tumor_RNAseq_rare_variants_VCF',
                                       '/data/ccbb_internal/interns/Carlo/test_vcf_out/Tumor_targeted_seq',
                                       '/data/ccbb_internal/interns/Carlo/test_vcf_out/normal_targeted_seq',
                                       '/data/ccbb_internal/interns/Carlo/test_vcf_out/normal_blood_WGS',
                                       '/data/ccbb_internal/interns/Carlo/test_vcf_out/somatic_mutect_old']


    def test_cytoBand_to_dict(self):

        cytoband_list = []
        for i in self.cytoBand:
            cyto = parser_models.CytoBand(i)
            cytoband_list.append(cyto.processed)

        for idx,_ in enumerate(self.cytoBand_dict):
            print(cytoband_list[idx])
            self.assertEqual(cytoband_list[idx].keys(), self.cytoBand_dict[idx].keys())


    #def test_expand_list_ids(self):
    #    self.assertEqual(utilities.expand_list(self.ids_to_expand), self.expanded_ids)


#HGVS_id creation not tested since it comes straight from myvariant.info's implementation: tested before.
    def test_annovar_subprocess_naming(self):

        output_csv_path = []
        for i in range(0, len(self.input_names)):
            if '/' in self.input_names[i]:
                output_name = self.input_names[i].split('/')[-1].split('.')[0]
            else:
                output_name = self.input_names[i].split('.')[0]
            output_csv_path.append(self.output_path + output_name)

        self.assertEqual(output_csv_path, self.expected_outfile_names)





    #def test_split_string(self):
    #    list1 = larger_df["ExonicFunc.knownGene"].dropna()


       # self.assertEqual()ls




if __name__ == '__main__':
    unittest.main()
