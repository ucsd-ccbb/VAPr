# standard libraries
import unittest
import os
import shlex
import logging

# standard libraries

# project-specific libraries
import subprocess
from VAPr.annotation_project import AnnotationProject
from VAPr import definitions
from VAPr.annovar import listen, AnnovarWrapper, AnnovarJobHandler

logger = logging.getLogger()
logger.setLevel(logging.INFO)

__author__ = 'Mazzaferro'


# TODO: Figure out why this test case seems to basically reimplement AnnovarWrapper
class TestAnnovar(unittest.TestCase):

    def setUp(self):

        self.base_dir = os.getcwd()
        self.files_input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir')
        self.samples_input_dir = os.path.join(self.base_dir, 'test_files/test_input_sample_dir')
        self.design_file_files = os.path.join(self.base_dir, 'test_files/design_file_by_file_name.csv')
        self.design_file_dirs = os.path.join(self.base_dir, 'test_files/design_file_by_dir_name.csv')
        self.output_csv_path_files = os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_files')
        self.output_csv_path_dirs = os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_dirs')
        self.annovar = os.path.join(self.base_dir, 'test_files/annovar_dir')
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}
        self.genome_build_version = 'hg19'
        self.relevant_dbs = ['genomicSuperDups']
        self.database_download_command_list = ['perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb genomicSuperDups /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb cytoBand /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb tfbsConsSites /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar knownGene /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar popfreq_all_20150413 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar cosmic70 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb targetScanS /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar clinvar_20161128 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar 1000g2015aug /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar esp6500siv2_all /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
                                               'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar nci60 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/']

        # use raw dictionary rather than using AnnotationProject to create input for project 1 and 2
        self.project_1 = AnnotationProject(self.files_input_dir,
                                           self.output_csv_path_files,
                                           self.annovar,
                                           self.project_data,
                                           build_ver='hg19')

        self.project_2 = AnnotationProject(self.files_input_dir,
                                           self.output_csv_path_dirs,
                                           self.annovar,
                                           self.project_data,
                                           design_file=self.design_file_dirs,
                                           build_ver='hg19')

        self.X45_22 = {'sample_names': ['X45'],
                      'num_samples_in_csv': 1,
                      'raw_vcf_file_full_path': os.path.join(self.base_dir, 'test_files/test_input_dir/X45/X45.raw.22.vcf'),
                      'csv_file_basename': 'X45.raw.22_annotated',
                      'vcf_file_basename': 'X45.raw.22.vcf',
                      'csv_file_full_path': os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_files/'),
                      'extra_data': None,
                      'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir/')}

    # def annovar_runner_stub(self, batch_jobs,  multisample=False):
    #
    #     handler = AnnovarJobHandler(batch_jobs, multisample, self.project_2.list_of_vcf_mapping_dicts)
    #     n_files_created = 0
    #     for index, job in enumerate(handler.chunkenize):
    #         logging.info('Job %i/%i sent for processing' % (index + 1, len(self.project_2.list_of_vcf_mapping_dicts) / batch_jobs + 1))
    #         n_files_created += len(job)
    #
    #         for idx, _map in enumerate(job):
    #             annotation_dir = _map['csv_file_full_path']
    #             if os.path.isdir(annotation_dir):
    #                 logging.info('Directory already exists for %s. '
    #                              'Writing output files there for file %s.' % (annotation_dir,
    #                                                                           _map['raw_vcf_file_full_path']))
    #             else:
    #                 os.makedirs(annotation_dir)
    #
    #             vcf_path = _map['raw_vcf_file_full_path']
    #             csv_path = os.path.join(_map['csv_file_full_path'], _map['csv_file_basename'])
    #             my_cmd = ['echo'] + ['test file for sample %s, job # %i' % (_map['sample_names'][0], idx)]
    #             with open(csv_path + '.txt', "w") as outfile:
    #                 subprocess.call(my_cmd, stdout=outfile)
    #
    #         logging.info('Annovar jobs submitted for %i files: %s' % (len(job),
    #                      ', '.join([os.path.basename(i['raw_vcf_file_full_path']) for i in job])))
    #
    #         listen(self.output_csv_path_dirs, len(job), n_files_created)
    #         logging.info('Finished running Annovar on this batch')
    #
    #     logging.info('Finished running Annovar on all files')
    #
    # def annovar_runner_stub_no_des_file(self, batch_jobs,  multisample=False):
    #
    #     handler = AnnovarJobHandler(batch_jobs, multisample, self.project_1.list_of_vcf_mapping_dicts)
    #     n_files_created = 0
    #     for index, job in enumerate(handler.chunkenize):
    #         logging.info('Job %i/%i sent for processing' % (index + 1, len(self.project_1.list_of_vcf_mapping_dicts) / batch_jobs + 1))
    #         n_files_created += len(job)
    #
    #         for idx, _map in enumerate(job):
    #             annotation_dir = _map['csv_file_full_path']
    #             if os.path.isdir(annotation_dir):
    #                 logging.info('Directory already exists for %s. '
    #                              'Writing output files there for file %s.' % (annotation_dir,
    #                                                                           _map['raw_vcf_file_full_path']))
    #             else:
    #                 os.makedirs(annotation_dir)
    #
    #             csv_path = os.path.join(_map['csv_file_full_path'], _map['csv_file_basename'])
    #
    #             my_cmd = ['echo'] + ['test file for sample %s, job # %i' % (_map['sample_names'], idx)]
    #             with open(csv_path + '.txt', "w") as outfile:
    #                 subprocess.call(my_cmd, stdout=outfile)
    #
    #         logging.info('Annovar jobs submitted for %i files: %s' % (len(job),
    #                      ', '.join([os.path.basename(i['raw_vcf_file_full_path']) for i in job])))
    #
    #         listen(self.output_csv_path_files, len(job), n_files_created)
    #         logging.info('Finished running Annovar on this batch')
    #
    #     logging.info('Finished running Annovar on all files')
    #
    # def test___get_annovar_dbs_to_use_for_build_version(self):
    #     self.assertEqual(AnnovarWrapper.hg_19_databases.iterkeys().next(), self.relevant_dbs[0])
    #     self.assertRaises(NameError, lambda: AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None, self.X45_22, None, 'FakeGenome99')[:1])
    #
    # def test_ensure_input_validity(self):
    #     self.assertEqual(AnnotationProject(self.files_input_dir, self.output_csv_path_files, self.annovar, self.project_data, build_ver='hg19').list_of_vcf_mapping_dicts[0], self.X45_22)

    # def test_run_annovar(self):
    #     self.annovar_runner_stub(5)
    #     self.annovar_runner_stub_no_des_file(5)
    #
    # def test_annovar_job_handler_init(self):
    #     ann = AnnovarJobHandler(10, False, self.project_1.list_of_vcf_mapping_dicts)
    #     self.assertEqual(ann._next()[0]['raw_vcf_file_full_path'], self.X45_22['raw_vcf_file_full_path'])
    #
    # def test_annovar_job_handler_init_dirs(self):
    #     ann = AnnovarJobHandler(10, False, self.project_2.list_of_vcf_mapping_dicts)
    #     self.assertEqual(ann._next()[0]['raw_vcf_file_full_path'], self.X45_22['raw_vcf_file_full_path'])
    #
    #
    # def test__build_table_annovar_command_str(self):
    #     wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None, self.X45_22, None, self.genome_build_version)
    #     expected_cmd_str = 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/table_annovar.pl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_input_dir/X45/X45.raw.22.vcf /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/ -genome_build_version hg19 -out /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_out_csv_path/des_file_files -remove -protocol genomicSuperDups,cytoBand,tfbsConsSites,knownGene,popfreq_all_20150413,cosmic70,targetScanS,clinvar_20161128,1000g2015aug_all,esp6500siv2_all,nci60 -operation r,r,r,g,f,f,r,f,f,f,f -nastring . -otherinfo -vcfinput'
    #     output_cmd_str = wrapper._build_table_annovar_command_str(self.X45_22['raw_vcf_file_full_path'], self.output_csv_path_files)
    #     self.assertEqual(expected_cmd_str, output_cmd_str)

    # def test__build_table_annovar_command_str_multisample(self):
    #     wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None, self.X45_22, None, self.genome_build_version)
    #     expected_cmd_str = 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/table_annovar.pl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_input_dir/X45/X45.raw.22.vcf /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/ -genome_build_version hg19 -out /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_out_csv_path/des_file_files -remove -protocol genomicSuperDups,cytoBand,tfbsConsSites,knownGene,popfreq_all_20150413,cosmic70,targetScanS,clinvar_20161128,1000g2015aug_all,esp6500siv2_all,nci60 -operation r,r,r,g,f,f,r,f,f,f,f -nastring . -otherinfo -vcfinput -format vcf4 -allsample -withfreq'
    #     output_cmd_str = wrapper._build_table_annovar_command_str(self.X45_22['raw_vcf_file_full_path'], self.output_csv_path_files, vcf_is_multisample=True)
    #     self.assertEqual(expected_cmd_str, output_cmd_str)
    #
    # # test _build_annovar_database_download_command_str(True, self.relevant_dbs)

    def test__build_annovar_database_download_command_str(self):
        wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None, self.X45_22,
                                 None, self.genome_build_version)
        output_cmd_str = wrapper._build_annovar_database_download_command_str(self.relevant_dbs)
        self.assertEqual([self.database_download_command_list[0]], output_cmd_str)

    def test__build_annovar_database_download_command_str_none(self):
        wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None, self.X45_22,
                                 None, self.genome_build_version)
        output_cmd_str = wrapper._build_annovar_database_download_command_str(None)
        self.assertListEqual(self.database_download_command_list, output_cmd_str)

    # def test_annovar_db_download(self):
    #     self.fail("Test not implemented")