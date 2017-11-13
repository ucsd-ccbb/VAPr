# standard libraries
import unittest
import os
import logging

# standard libraries

# project-specific libraries
from VAPr.vcf_merge import MergeVcfs
from VAPr.vcf_mappings_maker import SingleVcfFileMappingMaker

logger = logging.getLogger()
logger.setLevel(logging.INFO)

__author__ = 'Adam Mark'


# TODO: Figure out why this test case seems to basically reimplement AnnovarWrapper
class TestMergeVcfs(unittest.TestCase):

    def setUp(self):

        self.base_dir = os.getcwd()
        self.files_input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir/G1000')
        self.output_dir = os.path.join(self.base_dir, 'test_files/test_output_dir/G1000')
        self.merged_vcf = "merged"
        self.annovar = os.path.join(self.base_dir, 'test_files/annovar_dir')
        self.vcf_file_extension = ".vcf.gz"

        self.vcf_files = [os.path.join(self.files_input_dir, 'HG00096.vcf'),
                          os.path.join( self.files_input_dir, 'HG00097.vcf')]
        self.design_file = os.path.join(self.base_dir, 'test_files/design_file_G1000_files.csv')
        with open(self.design_file, 'w') as file:
            file.write('Sample_Names\n')
            for vcf in self.vcf_files:
                file.write(vcf)
                file.write('\n')

        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}
        self.genome_build_version = 'hg19'

        self.vcf_mapping_dict = {'sample_names': ['HG00096', 'HG00097'], 'num_samples_in_csv': 2,
                                 'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                                        'test_files/test_output_dir/G1000/merged.vcf'),
                                         'csv_file_basename': 'merged_annotated', 'vcf_file_basename': 'merged.vcf',
                                         'csv_file_full_path': os.path.join(self.base_dir,
                                                                            'test_files/test_output_dir/G1000/'),
                                         'extra_data': None,
                                         'vcf_sample_dir': os.path.join(self.base_dir,
                                                                        'test_files/test_input_dir/G1000/')}
        self.merge_command = 'bcftools merge {0}'.format(" ".join(self.vcf_files))
        self.bgzip_command = 'bgzip -c {0}'.format(self.vcf_files[0])
        self.index_command = 'tabix -p vcf {0}'.format(self.vcf_files[0] + '.gz')

    def test_merge_vcfs(self):
        merged = MergeVcfs(input_dir=self.files_input_dir,
                           output_dir=self.output_dir,
                           design_file=None,
                           analysis_name=self.merged_vcf,
                           vcf_file_extension=".vcf.gz").merge_vcfs()
        self.assertDictEqual(merged, self.vcf_mapping_dict)


    def test_merge_vcfs_design_file(self):
        merged = MergeVcfs(input_dir=self.files_input_dir,
                           output_dir=self.output_dir,
                           design_file=self.design_file,
                           analysis_name=self.merged_vcf,
                           vcf_file_extension=".vcf.gz").merge_vcfs()
        self.assertDictEqual(merged, self.vcf_mapping_dict)

    def test__build_merge_vcf_command_str(self):
        merge_command = MergeVcfs(input_dir=self.files_input_dir,
                           output_dir=self.output_dir,
                           design_file=self.design_file,
                           analysis_name=self.merged_vcf,
                           vcf_file_extension=".vcf.gz")._build_merge_vcf_command_str(self.vcf_files)
        self.assertEqual(self.merge_command, merge_command)

    def test__build_bgzip_vcf_command_str(self):
        bgzip_command = MergeVcfs(input_dir=self.files_input_dir,
                           output_dir=self.output_dir,
                           design_file=self.design_file,
                           analysis_name=self.merged_vcf,
                           vcf_file_extension=".vcf.gz")._build_bgzip_vcf_command_str(self.vcf_files[0])
        self.assertEqual(self.bgzip_command, bgzip_command)

    def test__build_index_vcf_command_str(self):
        index_command = MergeVcfs(input_dir=self.files_input_dir,
                           output_dir=self.output_dir,
                           design_file=self.design_file,
                           analysis_name=self.merged_vcf,
                           vcf_file_extension=".vcf.gz")._build_index_vcf_command_str(self.vcf_files[0] + '.gz')
        self.assertEqual(self.index_command, index_command)