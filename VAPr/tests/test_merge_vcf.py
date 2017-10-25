# standard libraries
import unittest
import os
import logging

# standard libraries

# project-specific libraries
from VAPr.annotation_project import AnnotationProject
from VAPr.annovar import AnnovarWrapper, AnnovarJobHandler
from VAPr.vcf_merge import MergeVcfs
from VAPr.vcf_mappings_maker import SingleVcfFileMappingMaker

logger = logging.getLogger()
logger.setLevel(logging.INFO)

__author__ = 'Adam Mark'


# TODO: Figure out why this test case seems to basically reimplement AnnovarWrapper
class TestMergeVcfs(unittest.TestCase):

    def setUp(self):

        self.base_dir = os.getcwd()
        self.files_input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir')
        self.output_dir = os.path.join(self.base_dir, 'test_files/test_out_csv_path/merged')
        self.merged_vcf = "merged.vcf"
        self.annovar = os.path.join(self.base_dir, 'test_files/annovar_dir')
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}
        self.genome_build_version = 'hg19'

        self.vcf_mapping_dicts_files = [{'sample_names': ['X45'], 'num_samples_in_csv': 1,
                                         'raw_vcf_file_full_path': os.path.join(self.base_dir, 'test_files/test_input_dir/X45/X45.raw.21.vcf'),
                                         'csv_file_basename': 'X45.raw.21_annotated', 'vcf_file_basename': 'X45.raw.21.vcf',
                                         'csv_file_full_path': os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_files/'),
                                         'extra_data': None,
                                         'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir')},
                                        {'sample_names': ['X7'], 'num_samples_in_csv': 1,
                                         'raw_vcf_file_full_path': os.path.join(self.base_dir, 'test_files/test_input_dir/X7/X7.raw.22.vcf'),
                                         'csv_file_basename': 'X7.raw.22_annotated', 'vcf_file_basename': 'X7.raw.22.vcf',
                                         'csv_file_full_path': os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_files/'),
                                         'extra_data': None,
                                         'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir')},
                                        ]

    # def test_merge_vcfs(self):
    #     print(AnnotationProject(input_dir=self.files_input_dir,
    #                             output_dir=self.output_dir,
    #                             merged_vcf_name=self.merged_vcf,
    #                             annovar_path=self.annovar,
    #                             mongo_db_and_collection_names_dict=self.project_data,
    #                             design_file=None, build_ver=self.genome_build_version, split_vcf=False).merged_vcf_mapping_dict)

    def test_merge_vcfs(self):
        print(MergeVcfs(input_dir=self.files_input_dir,
                                output_dir=self.output_dir,
                                list_of_vcf_mapping_dicts=self.vcf_mapping_dicts_files,
                                merged_vcf_name=self.merged_vcf).merge_vcfs())


    # def test_single_vcf_mapping_dict(self):
    #     print(SingleVcfFileMappingMaker(single_input_file_path=os.path.join(self.output_dir, self.merged_vcf),
    #                               input_dir=self.files_input_dir,
    #                               out_dir=self.output_dir,
    #                               sample_id='infer',
    #                               sample_id_type='files',
    #                               extra_data=None))