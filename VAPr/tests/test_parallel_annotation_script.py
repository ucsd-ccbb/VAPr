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

        self.base_dir = '/Volumes/Carlo_HD1/CCBB/VAPr_files/'
        self.input_dir = os.path.join(self.base_dir, 'vcf_files/multi_sample')
        self.design_file_dirs = os.path.join(self.base_dir, 'design_file_three_samples.csv')
        self.out_path = os.path.join(self.base_dir, 'csv_multisample')
        self.annovar = os.path.join(self.base_dir, '../annovar')
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}

        self.project = AnnotationProject(self.input_dir,
                                         self.out_path,
                                         self.annovar,
                                         self.project_data,
                                         design_file=self.design_file_dirs,
                                         build_ver='hg19',
                                         mongod_cmd='/usr/local/bin/mongod --dbpath /Volumes/Carlo_HD1/data/db/ '
                                                    '--storageEngine wiredTiger')

        self.annotator_wrapper = VariantParsing(self.input_dir,
                                                self.out_path,
                                                self.annovar,
                                                self.project_data,
                                                self.project.list_of_vcf_mapping_dicts,
                                                design_file=self.design_file_dirs,
                                                build_ver='hg19',
                                                mongod_cmd='mongod --dbpath /Volumes/Carlo_HD1/data/db/ '
                                                           '--storageEngine wiredTiger')

        self.x_7_raw_X = {'sample_names': ['X7'],
                          'num_samples_in_csv': 1,
                          'raw_vcf_file_full_path': '/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_files/multi_sample/X7/'
                                                    'X7.raw.X.vcf', 'csv_file_basename': 'X7.raw.X_annotated',
                          'vcf_file_basename': 'X7.raw.X.vcf',
                          'csv_file_full_path': '/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/X7',
                          'extra_data': {'libType': 'singleend', 'Tissue': 'lymphoblast', 'Patient': 'JNJ005',
                                         'Treatment': 'VPA', 'Condition': 'BD_lithium_responder'},
                          'vcf_sample_dir': '/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_files/multi_sample/X7'}

        self.sample_csv_vcf_tuple = (['X7'],
                                     '/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/X7/X7.raw.X_annotated.hg19_multianno.vcf',
                                     '/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/X7/X7.raw.X_annotated.hg19_multianno.txt',
                                     'VariantDatabase', 'collect',
                                     {'libType': 'singleend', 'Tissue': 'lymphoblast', 'Patient': 'JNJ005',
                                      'Treatment': 'VPA', 'Condition': 'BD_lithium_responder'})

        self.mapper_output = (['X7'],
                              '/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/X7/'
                              'X7.raw.X_annotated.hg19_multianno.vcf',
                              '/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/X7/'
                              'X7.raw.X_annotated.hg19_multianno.txt',
                              {'libType': 'singleend', 'Tissue': 'lymphoblast', 'Patient': 'JNJ005', 'Treatment': 'VPA',
                              'Condition': 'BD_lithium_responder'},
                              'VariantDatabase',
                              'collect',
                              0,
                              '/usr/local/bin/mongod --dbpath /Volumes/Carlo_HD1/data/db/ --storageEngine wiredTiger')

        self.minimapper = (['RAND'],
                           '/Users/carlomazzaferro/Desktop/sample_N15_vs_T15sample_mutect.13_annotated.hg19_multianno.vcf',
                           '/Users/carlomazzaferro/Desktop/sample_N15_vs_T15sample_mutect.13_annotated.hg19_multianno.txt',
                           'VariantDatabase',
                           'collect',
                           {'libType': 'singleend', 'Tissue': 'lymphoblast', 'Patient': 'JNJ005', 'Treatment': 'VPA',
                            'Condition': 'BD_lithium_responder'})

    def test_ensure_good_input(self):
        self.assertEqual(self.project.list_of_vcf_mapping_dicts[0], self.x_7_raw_X)

    def test_get_sample_csv_vcf_tuple(self):
        self.assertEqual(self.sample_csv_vcf_tuple, self.annotator_wrapper._get_sample_csv_vcf_tuple()[0])

    def test_parallel_annotator_mapper(self):
        hgvs = HgvsParser(self.sample_csv_vcf_tuple[1])
        csv_parsing = AnnovarTxtParser(self.sample_csv_vcf_tuple[2], samples=hgvs.samples,
                                       extra_data=self.sample_csv_vcf_tuple[5])
        num_lines = csv_parsing.num_lines
        n_steps = int(num_lines / self.annotator_wrapper.chunksize) + 1
        mapped = self.annotator_wrapper.parallel_annotator_mapper(self.sample_csv_vcf_tuple, n_steps,
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

    def test_parse_by_step(self):
        pass
        """
        hgvs = HgvsParser(self.sample_csv_vcf_tuple[1])
        csv_parsing = AnnovarTxtParser(self.sample_csv_vcf_tuple[2], samples=hgvs.samples,
                                extra_data=self.sample_csv_vcf_tuple[5])
        num_lines = csv_parsing.num_lines
        n_steps = int(num_lines / self.annotator_wrapper.chunksize) + 1
        mapped = self.annotator_wrapper.parallel_annotator_mapper(self.sample_csv_vcf_tuple, n_steps,
                                                                  extra_data=self.sample_csv_vcf_tuple[5])

        for i in range(mapped[-1][10]):
            parse_by_step(mapped[i])
        """
    def test_parse_by_step_2(self):
        hgvs = HgvsParser(self.minimapper[1])
        csv_parsing = AnnovarTxtParser(self.minimapper[2], samples=hgvs.samples,
                                       extra_data=self.minimapper[5])
        num_lines = csv_parsing.num_lines
        n_steps = int(num_lines / self.annotator_wrapper.chunksize) + 1
        mapped = self.annotator_wrapper.parallel_annotator_mapper(self.minimapper, n_steps,
                                                                  extra_data=self.minimapper[5],
                                                                  mongod_cmd=self.mapper_output[7])

        for i in range(mapped[-1][-2]):
            parse_by_step(mapped[i])