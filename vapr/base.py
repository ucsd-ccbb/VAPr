from __future__ import division, print_function
import os
import sys
import pandas
import logging
from vapr.annovar import AnnovarWrapper
from vapr.parsers import VariantParsing
from vapr.ingester import Ingester
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class ProjectData(object):

    supported_build_vers = ['hg19', 'hg18', 'hg38']

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 split_vcf=False,
                 design_file=None,
                 build_ver=None):

        """
        Base class for the project instantiation classes that will run the necessary functions. Mostly a container.
        """

        self.input_dir = input_dir
        self.output_csv_path = output_dir
        self.design_file = design_file
        self.annovar = annovar_path
        self.project_data = project_data
        self.buildver = self.check_ver(build_ver)
        self.times_called = 0
        self.mapping = self.get_mapping()
        self.split = split_vcf

    def get_mapping(self):

        if self.design_file:
            design_df = pandas.read_csv(self.design_file)
            sample_to_vcf_map = self.check_file_existance_return_mapping(design_df)

        else:
            sample_to_vcf_map = self.check_file_existance_return_mapping_no_design_file()

        return sample_to_vcf_map

    @staticmethod
    def listdir_fullpath(d):
        return [os.path.join(d, f) for f in os.listdir(d)]

    def check_file_existance_return_mapping_no_design_file(self):
        organizer = Ingester(self.input_dir, self.output_csv_path)
        organizer.walk()
        return organizer.mapping_list

    def check_file_existance_return_mapping(self, design_df):

        organizer = Ingester(self.input_dir, self.output_csv_path)
        organizer.digest_design_file(design_df)
        return organizer.mapping_list

    def check_ver(self, build_ver):

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

    """
        def split_vcf(self):
        if self.split:
            new_input_dir = os.path.join(self.input_dir, 'split_files')
            # self.input_dir = new_input_dir
            os.mkdir(os.path.join(new_input_dir))
            for sample in self.mapping.keys():
                os.mkdir(os.path.join(new_input_dir, sample))
                for file_pair in self.mapping[sample]['vcf_csv']:
                    vcf_reader = vcf.Reader(filename=os.path.join(self.input_dir, sample, file_pair[0]))

            vcf_writer = vcf.Writer(open('/dev/null', 'w'), vcf_reader)
            for record in vcf_reader:
                vcf_writer.write_record(record)
    """


class AnnotationProject(ProjectData):

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None):

        """ Class that implements the API and the major annotation/saving methods  """

        super(AnnotationProject, self).__init__(input_dir,
                                                output_dir,
                                                annovar_path,
                                                project_data,
                                                design_file=design_file,
                                                build_ver=build_ver)

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
                                                build_ver=self.buildver)

    def download_dbs(self, all_dbs=True, dbs=None):
        """ Wrapper around Annovar database downloading function """
        self.annovar_wrapper.download_dbs(all_dbs=all_dbs, dbs=dbs)

    def run_annovar(self, multisample=False):
        """ Wrapper around multiprocess Annovar annotation  """
        self.annovar_wrapper.run_annovar(multisample=multisample)

    def annotate_and_save(self, buffer_vars=False, verbose=2):
        """ Wrapper around annotation runner  """
        self.annotator_wrapper.annotate_and_saving(buffer_vars=buffer_vars, verbose=verbose)

    def parallel_annotation_and_saving(self, n_processes=4, verbose=1):
        """ Wrapper around parallel annotation multiprocess runner  """
        self.annotator_wrapper.parallel_annotation(n_processes=n_processes, verbose=verbose)


if __name__ == '__main__':

    IN_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_files/multi_sample/"
    OUT_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample_2/"
    ing = Ingester(IN_PATH, OUT_PATH)
    ing.walk()
    print(ing.mapping_list)


    """

    # Directory of input files to be annotated
    IN_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_files/multi_sample/"

    # Output file directory
    OUT_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/"
    # Location of your annovar dowload. The folder should contain the following files/directories:
    ANNOVAR_PATH = '/Volumes/Carlo_HD1/CCBB/annovar/'

    # Design File (optional)
    # design_file = '/Volumes/Carlo_HD1/CCBB/VAPr_files/guorong_single_sample.csv'

    # Databse and Collection names (optional)
    proj_data = {'db_name': 'VariantDBMultiBenchmark',
                 'project_name': 'collect'}

    Project = AnnotationProject(IN_PATH,
                                OUT_PATH,
                                ANNOVAR_PATH,
                                proj_data,
                                # design_file=design_file,
                                build_ver='hg19')
    """