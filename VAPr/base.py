from __future__ import division, print_function
import os
import sys
import pandas
import logging
from VAPr.annovar import AnnovarWrapper
from VAPr.parsers import VariantParsing
from VAPr.ingester import Ingester
from pymongo import MongoClient

logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'


class ProjectData(object):
    """ Base class for the project instantiation classes that will run the necessary functions. Mostly a container. """

    supported_build_vers = ['hg19', 'hg18', 'hg38']

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 split_vcf=False,
                 design_file=None,
                 build_ver=None,
                 mongod_cmd=None):

        self.input_dir = input_dir
        self.output_csv_path = output_dir
        self.design_file = design_file
        self.annovar = annovar_path
        self.project_data = project_data
        self.buildver = self.check_ver(build_ver)
        self.times_called = 0
        self.mapping = self.get_mapping()
        self.split = split_vcf
        self.mongod = mongod_cmd

    def get_mapping(self):
        """
        Get mapping of vcf file

        Each vcf file has a dictionary of values associated with it of the following form:

             mapping = {'raw_vcf_file_full_path': full path of vcf file,
                        'vcf_file_basename': vcf file base name,
                        'csv_file_basename': csv file base name,
                        'sample_names': samples in vcf file,
                        'num_samples_in_csv': number of samples,
                        'csv_file_full_path': output of sample file where csv files will be stored,
                        'vcf_sample_dir': directory of sample where vcf file lives,
                        'extra_data': None or some dictionary of data contained in design file
                        }

        Each vcf file will be moved to a directory named after the samples it contains

        """
        if self.design_file:
            design_df = pandas.read_csv(self.design_file)
            sample_to_vcf_map = self.check_file_existance_return_mapping(design_df)

        else:
            sample_to_vcf_map = self.check_file_existance_return_mapping_no_design_file()

        return sample_to_vcf_map

    @staticmethod
    def listdir_fullpath(d):
        """ Helper function to list full path of files in directory """
        return [os.path.join(d, f) for f in os.listdir(d)]

    def check_file_existance_return_mapping_no_design_file(self):
        """ Ingest all files in specified directory and return mapping """

        organizer = Ingester(self.input_dir, self.output_csv_path)
        organizer.walk()
        return organizer.mapping_list

    def check_file_existance_return_mapping(self, design_df):
        """ Ingest design file and the directories/files referenced in it and return mapping """

        organizer = Ingester(self.input_dir, self.output_csv_path)
        organizer.digest_design_file(design_df)
        return organizer.mapping_list

    def check_ver(self, build_ver):
        """ Make sure genome version is acceptable """

        if not build_ver:
            self.buildver = 'hg19'  # Default genome build version

        if build_ver not in self.supported_build_vers:
            raise ValueError('Build version must not recognized. Supported builds are'
                             ' %s, %s, %s' % (self.supported_build_vers[0],
                                              self.supported_build_vers[1],
                                              self.supported_build_vers[2]))
        else:
            self.buildver = build_ver

        return self.buildver


class AnnotationProject(ProjectData):

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None,
                 mongod_cmd=None):

        """ Class that implements the API and the major annotation/saving methods  """

        super(AnnotationProject, self).__init__(input_dir,
                                                output_dir,
                                                annovar_path,
                                                project_data,
                                                design_file=design_file,
                                                build_ver=build_ver,
                                                mongod_cmd=mongod_cmd)

        self.annovar_wrapper = AnnovarWrapper(self.input_dir,
                                              self.output_csv_path,
                                              self.annovar,
                                              self.project_data,
                                              self.mapping,
                                              design_file=self.design_file,
                                              build_ver=self.buildver)

        self.annotator_wrapper = VariantParsing(self.input_dir,
                                                self.output_csv_path,
                                                self.annovar,
                                                self.project_data,
                                                self.mapping,
                                                design_file=self.design_file,
                                                build_ver=self.buildver,
                                                mongod_cmd=self.mongod)

    def download_dbs(self, all_dbs=True, dbs=None):
        """ Wrapper around Annovar database downloading function """
        self.annovar_wrapper.download_dbs(all_dbs=all_dbs, dbs=dbs)

    def run_annovar(self, batch_jobs=10, multisample=False):
        """ Wrapper around multiprocess Annovar annotation  """
        self.annovar_wrapper.run_annovar(batch_jobs=batch_jobs, multisample=multisample)

    def annotate_and_save(self, buffer_vars=False, verbose=2):
        """ Wrapper around annotation runner. Deprecated in favour of parallel parsing """
        self.annotator_wrapper.annotate_and_saving(buffer_vars=buffer_vars, verbose=verbose)

    def parallel_annotation_and_saving(self, n_processes=4, verbose=1):
        """ Wrapper around parallel annotation multiprocess runner  """
        self.annotator_wrapper.parallel_annotation(n_processes=n_processes, verbose=verbose)

    def quick_annotate_and_save(self, n_processes=8):
        """ Wrapper around parallel annotation multiprocess runner using MyVariant solely """
        self.annotator_wrapper.quick_annotate_and_save(n_processes=n_processes)

    def write_output_files_by_sample(self):
        """ Wrapper around function that implemts the writing of csv files for each sample in collection """
        self.annotator_wrapper.generate_output_files_by_sample()

