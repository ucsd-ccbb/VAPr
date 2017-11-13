# standard libraries
import unittest
import os
from VAPr.vcf_mappings_maker import SingleVcfFileMappingMaker

__author__ = 'Mark'


class TestFunctions(unittest.TestCase):
    # TODO: no tests for _store_mapping_for_single_vcf

    def setUp(self):
        self.base_dir = os.getcwd()
        self.input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir/G1000')
        self.output_dir = os.path.join(self.base_dir, 'test_files/test_output_dir/G1000')
        self.vcf_path = os.path.join(self.base_dir, 'test_files/test_output_dir/G1000/merged.vcf')

        self.vcf_mapping_dict = {'sample_names': ['HG00097', 'HG00096'],
                                 'num_samples_in_csv': 2,
                                 'raw_vcf_file_full_path': self.vcf_path,
                                 'csv_file_basename': 'merged_annotated',
                                 'vcf_file_basename': 'merged.vcf',
                                 'csv_file_full_path': self.output_dir,
                                 'extra_data': None,
                                 'vcf_sample_dir': self.input_dir}

        self.map = SingleVcfFileMappingMaker(self.vcf_path, self.input_dir, self.output_dir)

    def test__add_extra_data(self):
        self.assertEqual(self.map.vcf_mapping_dict['extra_data'], self.vcf_mapping_dict['extra_data'])

    def test__fill_sample_names(self):
        self.assertListEqual(self.map.vcf_mapping_dict['sample_names'], self.vcf_mapping_dict['sample_names'])

    def test__fill_csv_file_basename(self):
        self.assertEqual(self.map.vcf_mapping_dict['csv_file_basename'], self.vcf_mapping_dict['csv_file_basename'])
