# standard libraries
import logging
import os
import unittest

# project-specific libraries
from VAPr.annovar_runner import AnnovarWrapper

logger = logging.getLogger()
logger.setLevel(logging.INFO)


class TestableAnnovarWrapper(AnnovarWrapper):
    _test_file_num = 0

    def _build_table_annovar_command_str(self, vcf_path, csv_path, vcf_is_multisample=False):
        self._test_file_num += 1
        temp_filename = "testable_annovar_wrapper_file_{0}.txt".format(self._test_file_num)
        temp_filepath = os.path.join(csv_path, temp_filename)
        cmd_string = "echo hello > {0}".format(temp_filepath)
        return cmd_string


# TODO: Figure out why this test case seems to basically reimplement AnnovarWrapper
class TestAnnovar(unittest.TestCase):
    def setUp(self):
        self.base_dir = os.getcwd()
        self.files_input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir')
        self.samples_input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir')
        self.single_input_dir = os.path.join(self.base_dir, 'test_files/test_input_dir/C90')
        self.design_file_files = os.path.join(self.base_dir, 'test_files/design_file_by_file_name.csv')
        self.design_file_dirs = os.path.join(self.base_dir, 'test_files/design_file_by_dir_name.csv')
        self.output_csv_path_files = os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_files')
        self.output_csv_path_dirs = os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_dirs')
        self.analysis_name = "annovar_test"
        self.annovar = os.path.join(self.base_dir, 'test_files/annovar_dir')
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}
        self.genome_build_version = 'hg19'
        self.relevant_dbs = ['knownGene']
        self.database_download_command_list = [
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb genomicSuperDups /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb cytoBand /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb tfbsConsSites /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar knownGene /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/'
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar popfreq_all_20150413 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar cosmic70 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb targetScanS /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar clinvar_20161128 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar 1000g2015aug /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/',
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar esp6500siv2_all /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/'
            # 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/annotate_variation.pl -build hg19 -downdb -webfrom annovar nci60 /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/'
        ]

        self.expected_table_cmd_str = 'perl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/table_annovar.pl /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_input_dir/X45/X45.raw.22.vcf /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir/humandb/ -genome_build_version hg19 -out /Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_out_csv_path/des_file_files -remove -protocol knownGene -operation g -nastring . -otherinfo -vcfinput'

        self.vcf_mapping_dicts_files = [{'sample_names': ['X45'], 'num_samples_in_csv': 1,
                                         'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                                                'test_files/test_input_dir/X45/X45.raw.22.vcf'),
                                         'csv_file_basename': 'X45.raw.22_annotated',
                                         'vcf_file_basename': 'X45.raw.22.vcf',
                                         'csv_file_full_path': os.path.join(self.base_dir,
                                                                            'test_files/test_out_csv_path/des_file_files/'),
                                         'extra_data': None,
                                         'vcf_sample_dir': os.path.join(self.base_dir,
                                                                        'test_files/test_input_dir/X45')},
                                        {'sample_names': ['X7'], 'num_samples_in_csv': 1,
                                         'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                                                'test_files/test_input_dir/X7/X7.raw.22.vcf'),
                                         'csv_file_basename': 'X7.raw.Y_annotated',
                                         'vcf_file_basename': 'X7.raw.22.vcf',
                                         'csv_file_full_path': os.path.join(self.base_dir,
                                                                            'test_files/test_out_csv_path/des_file_files/'),
                                         'extra_data': None,
                                         'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir/X7')}

                                        ]

        self.vcf_mapping_dicts_dirs = [{'sample_names': ['X45'], 'num_samples_in_csv': 1,
                                        'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                                               'test_files/test_input_dir/X45/X45.raw.22.vcf'),
                                        'csv_file_basename': 'X45.raw.22_annotated',
                                        'vcf_file_basename': 'X45.raw.22.vcf',
                                        'csv_file_full_path': os.path.join(self.base_dir,
                                                                           'test_files/test_out_csv_path/des_file_dirs/X45'),
                                        'extra_data': {'libType': 'singleend', 'Tissue': 'lymphoblast',
                                                       'Patient': 'JNJ005', 'Treatment': 'Li',
                                                       'Condition': 'BD_lithium_responder'},
                                        'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir/X45')},
                                       {'sample_names': ['X7'], 'num_samples_in_csv': 1,
                                        'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                                               'test_files/test_input_dir/X7/X7.raw.X.vcf'),
                                        'csv_file_basename': 'X7.raw.Y_annotated', 'vcf_file_basename': 'X7.raw.X.vcf',
                                        'csv_file_full_path': os.path.join(self.base_dir,
                                                                           'test_files/test_out_csv_path/des_file_dirs/X7'),
                                        'extra_data': {'libType': 'singleend', 'Tissue': 'lymphoblast',
                                                       'Patient': 'JNJ005', 'Treatment': 'VPA',
                                                       'Condition': 'BD_lithium_responder'},
                                        'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir/X7')}

                                       ]

        self.X45_22 = {'sample_names': ['X45'],
                       'num_samples_in_csv': 1,
                       'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                              'test_files/test_input_dir/X45/X45.raw.22.vcf'),
                       'csv_file_basename': 'X45.raw.22_annotated',
                       'vcf_file_basename': 'X45.raw.22.vcf',
                       'csv_file_full_path': os.path.join(self.base_dir,
                                                          'test_files/test_out_csv_path/des_file_files/'),
                       'extra_data': None,
                       'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir/X45')}

        self.C90_22 = {'sample_names': ['C90'],
                       'num_samples_in_csv': 1,
                       'raw_vcf_file_full_path': os.path.join(self.base_dir,
                                                              'test_files/test_input_dir/C90/C90.raw.22.vcf'),
                       'csv_file_basename': 'C90.raw.22_annotated',
                       'vcf_file_basename': 'C90.raw.22.vcf',
                       'csv_file_full_path': os.path.join(self.base_dir,
                                                          'test_files/test_out_csv_path/des_file_files/'),
                       'extra_data': None,
                       'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir/C90')}

    def test__get_annovar_dbs_to_use_for_build_version(self):
        # use known-good genome version in init
        test_wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None,
                                      self.vcf_mapping_dicts_files, None, 'hg19')
        self.assertDictEqual(AnnovarWrapper.hg_19_databases, test_wrapper.annovar_dbs_to_use)

    def test__get_annovar_dbs_to_use_for_build_version_error(self):
        test_wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None,
                                      self.vcf_mapping_dicts_files, None, 'hg19')
        with self.assertRaises(NameError):
            # Note: would never access _genome_build_version directly from outside like this except for testing:
            # is wrapped in a property specifically to PREVENT this kind of corruption :)
            test_wrapper._genome_build_version = "Fakevalue"
            test_wrapper._get_annovar_dbs_to_use_for_build_version()

    # Ok, I don't know what this test code was trying to do, but it definitely shouldn't be testing anything in
    # AnnovarWrapper by calling it *through AnnotationProject*.  Needs rewrite.
    # def test_ensure_input_validity(self):
    #     self.assertEqual(
    #         AnnotationProject(self.files_input_dir, self.output_csv_path_files, self.analysis_name, self.annovar,
    #                           self.project_data, build_ver='hg19').list_of_vcf_mapping_dicts[0]['sample_names'][0],
    #         self.vcf_mapping_dicts_files[1]['sample_names'][0])
    #
    # def test_ensure_input_validity(self):
    #     self.assertEqual(
    #         AnnotationProject(self.single_input_dir, self.output_csv_path_files, self.analysis_name, self.annovar,
    #                           self.project_data, build_ver='hg19').list_of_vcf_mapping_dicts[0]['sample_names'][0],
    #         self.C90_22['sample_names'][0])

    def test_run_annovar(self):
        self.fail("test not implemented")
        # Ok, I don't know what this test code was trying to do, but it definitely shouldn't be testing run_annovar
        # in AnnovarWrapper by calling it *through AnnotationProject*.  Needs rewrite.
        # #     try:
        # #         wrapper = TestableAnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None, self.vcf_mapping_dicts_files, None, self.genome_build_version)
        # #         self.assertEqual(self.expected_table_cmd_str, wrapper._run_annovar_annotation(2))
        # #     finally:
        # #         # clean up test files here
        # Project = AnnotationProject(self.files_input_dir, self.output_csv_path_dirs, self.analysis_name, self.annovar,
        #                             self.project_data, build_ver="hg19")
        # Project._run_annovar_annotation(multisample=True)

    def test__build_table_annovar_command_str(self):
        wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None,
                                 self.vcf_mapping_dicts_files, None, self.genome_build_version)
        output_cmd_str = wrapper._build_table_annovar_command_str(self.X45_22['raw_vcf_file_full_path'],
                                                                  self.output_csv_path_files)
        self.assertEqual(self.expected_table_cmd_str, output_cmd_str)

    def test__build_table_annovar_command_str_multisample(self):
        wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None, self.X45_22,
                                 None, self.genome_build_version)
        expected_cmd_str = self.expected_table_cmd_str + ' -format vcf4 -allsample -withfreq'
        output_cmd_str = wrapper._build_table_annovar_command_str(self.X45_22['raw_vcf_file_full_path'],
                                                                  self.output_csv_path_files, vcf_is_multisample=True)
        self.assertEqual(expected_cmd_str, output_cmd_str)

    def test__build_annovar_database_download_command_str(self):
        wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None,
                                 self.vcf_mapping_dicts_files,
                                 None, self.genome_build_version, self.relevant_dbs)
        output_cmd_str = wrapper._build_annovar_database_download_command_str()
        self.assertEqual([self.database_download_command_list[0]], output_cmd_str)

    def test__build_annovar_database_download_command_str_none(self):
        wrapper = AnnovarWrapper(self.files_input_dir, self.output_csv_path_files, self.annovar, None,
                                 self.vcf_mapping_dicts_files,
                                 None, self.genome_build_version)
        output_cmd_str = wrapper._build_annovar_database_download_command_str()
        self.assertListEqual(self.database_download_command_list, output_cmd_str)
