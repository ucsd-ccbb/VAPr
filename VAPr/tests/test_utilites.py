import unittest
import sys
import os
sys.path.append("/Users/carlomazzaferro/Documents/Code/variantannotation-master/variantannotation")
from variantannotation import utilities
from variantannotation import csv_to_df


file_name = os.path.dirname(os.path.realpath('__file__')) + '/test_file.csv'
sample_list = csv_to_df.open_and_parse(file_name)
sample_df = csv_to_df.parse_to_df(sample_list)

larger_list = csv_to_df.open_and_parse(file_name)
larger_df = csv_to_df.parse_to_df(larger_list)


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
                                {'Arm': 'p', 'Band': 6, 'Chromosome': '1', 'Region': 3, 'Sub_Band': 33},
                                {'Arm': 'p', 'Band': 1, 'Chromosome': '16', 'Region': 1, 'Sub_Band': 1},
                                {'Arm': 'q', 'Band': 1, 'Chromosome': '16', 'Region': 1, 'Sub_Band': 2},
                                {'Arm': 'q', 'Band': 1, 'Chromosome': '16', 'Region': 2},
                                {'Arm': 'p', 'Band': 2, 'Chromosome': 'X', 'Region': 2, 'Sub_Band': 32},
                                {'Arm': 'p', 'Band': 2, 'Chromosome': 'X', 'Region': 2, 'Sub_Band': 2},
                                {'Arm': 'p', 'Band': 2, 'Chromosome': 'X', 'Region': 2, 'Sub_Band': 11},
                                {'Arm': 'q', 'Band': 2, 'Chromosome': 'X', 'Region': 1},
                                {'Arm': 'q', 'Band': 3, 'Chromosome': 'X', 'Region': 1, 'Sub_Band': 1}
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


    def test_split_cytoBand(self):

        cytolist = []
        for i in range(0, len(self.cytoBand)):
            cytolist.append(utilities.split_cytoband(self.cytoBand[i]))

        self.assertEqual(cytolist, self.cytoBand_split)

    def test_cytoBand_to_dict(self):

        cytodict = []
        for i in self.cytoBand:
            print cytodict
            cytodict.append(utilities.lists_to_dict(utilities.split_cytoband(i)))

        self.assertEqual(cytodict, self.cytoBand_dict)


    def test_expand_list_ids(self):
        self.assertEqual(utilities.expand_list(self.ids_to_expand), self.expanded_ids)


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
