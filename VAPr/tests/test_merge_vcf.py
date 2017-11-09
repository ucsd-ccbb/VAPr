# standard libraries
import unittest
import os
import logging

# standard libraries

# project-specific libraries
from VAPr.vcf_merge import MergeVcfs

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
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}
        self.genome_build_version = 'hg19'
        self.design_file = None
        self.vcf_file_extension = ".vcf.gz"
        #self.vcf_file_paths = [os.path.join(self.base_dir, 'test_files/test_input_dir/G1000/HG00096.vcf.gz', os.path.join(self.base_dir, 'test_files/test_input_dir/G1000/HG00097.vcf.gz']

    def test_merge_vcfs(self):
        merged = MergeVcfs(self.files_input_dir,
                          self.output_dir,
                          self.merged_vcf,
                          self.design_file,
                          self.vcf_file_extension).merge_vcfs()
        print(merged)
        #self.assertEqual(merged[0]['vcf_file_basename'], self.merged_vcf + ".vcf")

