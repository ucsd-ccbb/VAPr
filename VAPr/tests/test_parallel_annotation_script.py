# standard libraries
import unittest
import os
import logging

# project-specific libraries
from VAPr.annotation_project import AnnotationProject
from VAPr import definitions
from VAPr.parsers import HgvsParser, AnnovarTxtParser
from VAPr.parsers import VariantParsing
from VAPr.parsers import parse_by_step

logger = logging.getLogger()
logger.setLevel(logging.INFO)


__author__ = 'Mazzaferro'


class TestParallelAnnotationFunctions(unittest.TestCase):

    def setUp(self):

        self.base_dir = os.getcwd()
        self.files_input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir/G1000')
        self.output_dir = os.path.join(self.base_dir, 'test_files/test_output_dir/G1000')
        self.merged_vcf = "merged"
        self.annovar = os.path.join(self.base_dir, 'test_files/annovar_dir')
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}
        self.genome_build_version = 'hg19'
        self.fake_genome = "hg18"
        self.vcf_file_extension = ".vcf.gz"
        self.vcf_mapping_dict = {'sample_names': ['HG00096', 'HG00097'], 'num_samples_in_csv': 2,
                                         'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                                                'test_files/test_output_dir/G1000/merged.vcf'),
                                         'csv_file_basename': 'merged_annotated', 'vcf_file_basename': 'merged.vcf',
                                         'csv_file_full_path': os.path.join(self.base_dir,
                                                                            'test_files/test_output_dir/G1000'),
                                         'extra_data': None,
                                         'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir')}

        self.sample_csv_vcf_tuple = (['merged'],
                                     os.path.join(self.base_dir,
                                                  'test_files/test_out_csv_path/des_file_dirs/G1000/merged_annotated.hg19_multianno.vcf'),
                                     os.path.join(self.base_dir,
                                                  'test_files/test_out_csv_path/des_file_dirs/G1000/merged_annotated.hg19_multianno.txt'),
                                     'VariantDatabase', 'collect', None)

    def test_ensure_good_input(self):
        project = AnnotationProject(input_dir=self.files_input_dir,
                                    output_dir=self.output_dir,
                                    analysis_name=self.merged_vcf,
                                    vcf_file_extension=self.vcf_file_extension,
                                    annovar_path=self.annovar,
                                    build_ver=self.genome_build_version,
                                    mongo_db_and_collection_names_dict = self.project_data)
        #self.assertEqual(project.vcf_mapping_dict['raw_vcf_file_full_path'],
        #                 self.vcf_mapping_dict['raw_vcf_file_full_path'])

    def test_parallel_annotator_mapper(self):
        hgvs = HgvsParser(self.sample_csv_vcf_tuple[1])
        csv_parsing = AnnovarTxtParser(self.sample_csv_vcf_tuple[2], samples=hgvs.samples,
                                       extra_data=self.sample_csv_vcf_tuple[5])
        annotator_wrapper = VariantParsing(self.files_input_dir, self.output_dir, self.annovar,
                                                self.project_data, [self.vcf_mapping_dicts[0]],
                                                build_ver=self.genome_build_version)
        num_lines = csv_parsing.num_lines
        n_steps = int(num_lines / annotator_wrapper.chunksize) + 1
        mapped = annotator_wrapper.parallel_annotator_mapper(self.sample_csv_vcf_tuple, n_steps,
                                                                  extra_data=self.sample_csv_vcf_tuple[5],
                                                                  mongod_cmd=self.mapper_output[7])
        for i, _m in enumerate(mapped):
            self.assertEqual(self.mapper_output[0], _m[0])
            self.assertEqual(self.mapper_output[1], _m[1])
            self.assertEqual(self.mapper_output[2], _m[2])
            self.assertEqual(self.mapper_output[3], _m[3])
            self.assertEqual(self.mapper_output[4], _m[4])
            self.assertEqual(self.mapper_output[5], _m[5])
            self.assertEqual(i, _m[6])
            self.assertEqual(self.mapper_output[7], _m[7])




