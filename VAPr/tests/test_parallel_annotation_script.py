# standard libraries
import unittest
import os
import logging
import zipfile

# project-specific libraries
from VAPr.annotation_project import AnnotationProject
from VAPr import definitions
from VAPr.parsers import HgvsParser, TxtParser
from VAPr.parsers import VariantParsing
from VAPr.parsers import parse_by_step
from pprint import pprint


logger = logging.getLogger()
logger.setLevel(logging.INFO)


__author__ = 'Mazzaferro'


"""
Here are 3 sample vcf files, randomly picked from a collection of vcf files representing a single sample (X7), and
representing an area of a chromosome of a vcf file. Despite this not being an usual use-case, the test is implemented
to ensure that the annotation is going on correctly on multiple vcf files found in a specific directory which represents
a sample. This directory is also specified in the design file as the sample of the vcf files in question.
"""

# FILE: X7.raw.X.vcf (Chromosome X)
test_file_1_vcf_annotated = os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'
                                                      'X7.raw.X_annotated.hg19_multianno.vcf')
test_file_1_txt = os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'
                                            'X7.raw.X_annotated.hg19_multianno.txt')


# FILE: X7.raw.11.vcf (Chromosome 11)
test_file_2_vcf_annotated = os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'
                                                      'X7.raw.11_annotated.hg19_multianno.vcf')
test_file_2_txt = os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'
                                            'X7.raw.11_annotated.hg19_multianno.txt')


# FILE: X7.raw.7.vcf (Chromosome 7)
test_file_3_vcf_annotated = os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'
                                                      'X7.raw.7_annotated.hg19_multianno.vcf')
test_file_3_txt = os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'
                                            'X7.raw.7_annotated.hg19_multianno.txt')


class TestParallelAnnotationFunctions(unittest.TestCase):

    def setUp(self):
        self.unzip_all()   # files need to be posted on GitHub compressed since files > 100 MB are not allowed
        self.base_dir = os.getcwd()
        self.input_dir = os.path.join(self.base_dir, 'test_files/real_vcf_files')
        self.design_file_dirs = os.path.join(self.base_dir, 'test_files/design_file_single_sample.csv')
        self.out_path = os.path.join(self.base_dir, 'test_files/test_out_csv_path/real_files')
        self.annovar = os.path.join(self.base_dir, 'test_files/annovar_dir')
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

        self.x_7_raw_X = {'sample_names': ['X7'], 'num_samples_in_csv': 1,
                          'raw_vcf_file_full_path': os.path.join(os.getcwd(), 'test_files/real_vcf_files/X7/X7.raw.X.vcf'),
                          'csv_file_basename': 'X7.raw.X_annotated',
                          'vcf_file_basename': 'X7.raw.X.vcf',
                          'csv_file_full_path':  os.path.join(self.base_dir, 'test_files/test_out_csv_path/real_files/X7'),
                          'extra_data': {'libType': 'singleend', 'Tissue': 'lymphoblast', 'Patient': 'JNJ005',
                                         'Treatment': 'VPA', 'Condition': 'BD_lithium_responder'},
                          'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/real_vcf_files/X7')}

        self.sample_csv_vcf_tuple = (['X7'],
                                     test_file_1_vcf_annotated,
                                     test_file_1_txt,
                                     'VariantDatabase', 'collect',
                                     {'libType': 'singleend', 'Tissue': 'lymphoblast', 'Patient': 'JNJ005',
                                      'Treatment': 'VPA', 'Condition': 'BD_lithium_responder'})

        self.mapper_output = (['X7'],
                              test_file_1_vcf_annotated,
                              test_file_1_txt,
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

    def tearDown(self):
        os.remove(test_file_1_vcf_annotated)
        os.remove(test_file_1_txt)
        os.remove(test_file_2_vcf_annotated)
        os.remove(test_file_2_txt)
        os.remove(test_file_3_vcf_annotated)
        os.remove(test_file_3_txt)

    @staticmethod
    def unzip_all():
        zip_csv_1 = zipfile.ZipFile(test_file_1_txt + '.zip', 'r')
        zip_csv_1.extractall(os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'))
        zip_vcf_1_anno = zipfile.ZipFile(test_file_1_vcf_annotated + '.zip', 'r')
        zip_vcf_1_anno.extractall(os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'))

        zip_csv_2 = zipfile.ZipFile(test_file_2_txt + '.zip', 'r')
        zip_csv_2.extractall(os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'))
        zip_vcf_2_anno = zipfile.ZipFile(test_file_2_vcf_annotated + '.zip', 'r')
        zip_vcf_2_anno.extractall(os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'))

        zip_csv_3 = zipfile.ZipFile(test_file_3_txt + '.zip', 'r')
        zip_csv_3.extractall(os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'))
        zip_vcf_3_anno = zipfile.ZipFile(test_file_3_vcf_annotated + '.zip', 'r')
        zip_vcf_3_anno.extractall(os.path.join(os.getcwd(), 'test_files/test_out_csv_path/real_files/X7/'))

    def test_ensure_good_input(self):

        self.assertEqual(self.project.list_of_vcf_mapping_dicts[-1], self.x_7_raw_X)

    def test_get_sample_csv_vcf_tuple(self):
        pprint(self.sample_csv_vcf_tuple)
        pprint(self.annotator_wrapper._get_sample_csv_vcf_tuple()[2])
        self.assertEqual(self.sample_csv_vcf_tuple, self.annotator_wrapper._get_sample_csv_vcf_tuple()[2])

    def test_parallel_annotator_mapper(self):
        hgvs = HgvsParser(self.sample_csv_vcf_tuple[1])
        csv_parsing = TxtParser(self.sample_csv_vcf_tuple[2], samples=hgvs.samples,
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
        csv_parsing = TxtParser(self.sample_csv_vcf_tuple[2], samples=hgvs.samples,
                                extra_data=self.sample_csv_vcf_tuple[5])
        num_lines = csv_parsing.num_lines
        n_steps = int(num_lines / self.annotator_wrapper.chunksize) + 1
        mapped = self.annotator_wrapper.parallel_annotator_mapper(self.sample_csv_vcf_tuple, n_steps,
                                                                  extra_data=self.sample_csv_vcf_tuple[5])

        for i in range(mapped[-1][10]):
            parse_by_step(mapped[i])
        """
