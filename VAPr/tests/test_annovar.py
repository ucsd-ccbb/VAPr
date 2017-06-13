# standard libraries
import unittest
import os
import shlex
import logging

# project-specific libraries
import subprocess
from VAPr.base import AnnotationProject
from VAPr import definitions
from VAPr.annovar import listen, AnnovarJobHandler

logger = logging.getLogger()
logger.setLevel(logging.INFO)


__author__ = 'Mazzaferro'


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
                             'project_name': 'collect'}
        self.buildver = 'hg19'
        self.hg_18_databases = definitions.hg_18_databases
        self.hg_19_databases = definitions.hg_19_databases
        self.hg_38_databases = definitions.hg_38_databases
        self.databases = self.get_databases()
        self.project_1 = AnnotationProject(self.files_input_dir,
                                           self.output_csv_path_files,
                                           self.annovar,
                                           self.project_data,
                                           build_ver='hg19')

        self.project_2 = AnnotationProject(self.samples_input_dir,
                                           self.output_csv_path_dirs,
                                           self.annovar,
                                           self.project_data,
                                           design_file=self.design_file_dirs,
                                           build_ver='hg19')

        self.mini1 = {'sample_names': [],
                      'num_samples_in_csv': 0,
                      'raw_vcf_file_full_path': os.path.join(self.base_dir, 'test_files/test_input_dir/mini1.vcf'),
                      'csv_file_basename': 'mini1_annotated',
                      'vcf_file_basename': 'mini1.vcf',
                      'csv_file_full_path': os.path.join(self.base_dir, 'test_files/test_out_csv_path/des_file_files/'),
                      'extra_data': None,
                      'vcf_sample_dir': os.path.join(self.base_dir, 'test_files/test_input_dir/')}

    def annovar_runner_stub(self, batch_jobs,  multisample=False):

        handler = AnnovarJobHandler(batch_jobs, multisample, self.project_2.mapping)
        n_files_created = 0
        for index, job in enumerate(handler.chunkenize):
            logging.info('Job %i/%i sent for processing' % (index + 1, len(self.project_2.mapping)/batch_jobs + 1))
            n_files_created += len(job)

            for idx, _map in enumerate(job):
                annotation_dir = _map['csv_file_full_path']
                if os.path.isdir(annotation_dir):
                    logging.info('Directory already exists for %s. '
                                 'Writing output files there for file %s.' % (annotation_dir,
                                                                              _map['raw_vcf_file_full_path']))
                else:
                    os.makedirs(annotation_dir)

                vcf_path = _map['raw_vcf_file_full_path']
                csv_path = os.path.join(_map['csv_file_full_path'], _map['csv_file_basename'])
                cmd_string = self.build_annovar_command_str(vcf_path, csv_path, multisample=multisample)
                args = shlex.split(cmd_string)

                my_cmd = ['echo'] + ['test file for sample %s, job # %i' % (_map['sample_names'][0], idx)]
                with open(csv_path + '.txt', "w") as outfile:
                    subprocess.call(my_cmd, stdout=outfile)

            logging.info('Annovar jobs submitted for %i files: %s' % (len(job),
                         ', '.join([os.path.basename(i['raw_vcf_file_full_path']) for i in job])))

            listen(self.output_csv_path_dirs, len(job), n_files_created)
            logging.info('Finished running Annovar on this batch')

        logging.info('Finished running Annovar on all files')

    def annovar_runner_stub_no_des_file(self, batch_jobs,  multisample=False):

        handler = AnnovarJobHandler(batch_jobs, multisample, self.project_1.mapping)
        n_files_created = 0
        for index, job in enumerate(handler.chunkenize):
            logging.info('Job %i/%i sent for processing' % (index + 1, len(self.project_1.mapping)/batch_jobs + 1))
            n_files_created += len(job)

            for idx, _map in enumerate(job):
                annotation_dir = _map['csv_file_full_path']
                if os.path.isdir(annotation_dir):
                    logging.info('Directory already exists for %s. '
                                 'Writing output files there for file %s.' % (annotation_dir,
                                                                              _map['raw_vcf_file_full_path']))
                else:
                    os.makedirs(annotation_dir)

                vcf_path = _map['raw_vcf_file_full_path']
                csv_path = os.path.join(_map['csv_file_full_path'], _map['csv_file_basename'])
                cmd_string = self.build_annovar_command_str(vcf_path, csv_path, multisample=multisample)
                args = shlex.split(cmd_string)

                my_cmd = ['echo'] + ['test file for sample %s, job # %i' % (_map['sample_names'], idx)]
                with open(csv_path + '.txt', "w") as outfile:
                    subprocess.call(my_cmd, stdout=outfile)

            logging.info('Annovar jobs submitted for %i files: %s' % (len(job),
                         ', '.join([os.path.basename(i['raw_vcf_file_full_path']) for i in job])))

            listen(self.output_csv_path_files, len(job), n_files_created)
            logging.info('Finished running Annovar on this batch')

        logging.info('Finished running Annovar on all files')

    def build_annovar_command_str(self, _vcf, _csv, multisample=False):
        """ Concatenate command string arguments for Annovar jobs """

        # TODO: check for newer version of databases

        dbs = ",".join(list(self.databases.keys()))
        dbs_args = ",".join(list(self.databases.values()))

        if '1000g2015aug' in dbs:
            dbs = dbs.replace('1000g2015aug', '1000g2015aug_all')
        command = " ".join(['perl', os.path.join(self.annovar, 'table_annovar.pl'), _vcf,
                            os.path.join(self.annovar, 'humandb/'), '-buildver', self.buildver, '-out',
                            _csv, '-remove -protocol', dbs,  '-operation',
                            dbs_args, '-nastring .', '-otherinfo -vcfinput'])
        if multisample:
            command += ' -format vcf4 -allsample -withfreq'

        return command

    def get_databases(self):

        if self.buildver == 'hg18':
            databases = self.hg_18_databases
        elif self.buildver == 'hg19':
            databases = self.hg_19_databases
        else:
            databases = self.hg_38_databases

        return databases

    def test_ensure_input_validity(self):
        self.assertEqual(self.project_1.mapping[0], self.mini1)

    def test_run_annovar(self):
        self.annovar_runner_stub(5)
        self.annovar_runner_stub_no_des_file(5)

    def test_annovar_job_handler(self):
        ann = AnnovarJobHandler(10, False, self.project_1.mapping)
        self.assertEqual(ann._next()[0], self.mini1)

    def test_annovar_job_handler_dirs(self):
        ann = AnnovarJobHandler(10, False, self.project_2.mapping)
        print(len(ann._next()))
        print(len(ann._next()))

    def test_annovar_db_download(self):
        pass  # fill
