from __future__ import division, print_function
import sys
import pandas
import logging
from VAPr.vcf_merge import MergeVcfs
from VAPr.annovar import AnnovarWrapper
from VAPr.parsers import VariantParsing
import VAPr.vcf_mappings_maker

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'

# TODO: is this really necessary?
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


# TODO: In long term, would prefer to pass db_name and collection_name as individual arguments to AnnotationProject
# init method.  However, for now, to avoid breaking existing interface, I am adding this function to ensure
# project_data can easily be created correctly.
def make_mongo_db_and_collection_names_dict(mongo_db_name, mongo_collection_name):
    return {'db_name': mongo_db_name, 'collection_name': mongo_collection_name}


class AnnotationProject:
    # TODO: Consider making this an enum?
    supported_build_vers = ['hg19', 'hg38']

    def __init__(self, input_dir, output_dir, analysis_name, annovar_path, mongo_db_and_collection_names_dict,
                 design_file=None, build_ver=None, mongod_cmd=None, split_vcf=False):
        # type: (object, object, object, object, object, object, object, object, object) -> object
        # TODO: Why is the above here?
        """ Class that implements the API and the major annotation/saving methods  
        :rtype: object
        """

        self.input_dir = input_dir
        self.output_csv_path = output_dir
        self.analysis_name = analysis_name
        self.design_file = design_file
        self.annovar_path = annovar_path
        self.mongo_db_and_collection_names_dict = mongo_db_and_collection_names_dict
        self.genome_build_version = self._validate_genome_version(build_ver)
        self.times_called = 0
        self.split = split_vcf
        self.mongod = mongod_cmd

        vcf_file_path_list = VAPr.vcf_mappings_maker.get_vcf_file_paths_list(input_dir, design_file)

        # return vcf mapping dict
        self.vcf_mapping_dict = MergeVcfs(self.input_dir,
                                          self.output_csv_path,
                                          vcf_file_path_list,
                                          self.analysis_name).merge_vcfs()

        # TODO: These two calls takes in exactly the same parameters; Consider storing them to an object and passing
        # object around?
        self.annovar_wrapper = AnnovarWrapper(self.input_dir, self.output_csv_path, self.annovar_path,
                                              self.mongo_db_and_collection_names_dict, self.vcf_mapping_dict,
                                              design_file=self.design_file,
                                              genome_build_version=self.genome_build_version)

        self.annotator_wrapper = VariantParsing(self.input_dir, self.output_csv_path, self.annovar_path,
                                                self.mongo_db_and_collection_names_dict, self.vcf_mapping_dict,
                                                design_file=self.design_file,
                                                build_ver=self.genome_build_version,
                                                mongod_cmd=self.mongod)

    def download_dbs(self):
        """ Wrapper around Annovar database downloading function """
        self.annovar_wrapper.download_dbs()

    def run_annovar(self, multisample=False):
        """ Wrapper around multiprocess Annovar annotation  """
        self.annovar_wrapper.run_annovar(vcf_is_multisample=multisample)

    def parallel_annotation(self, n_processes=4, verbose=1):
        """ Wrapper around parallel annotation multiprocess runner  """
        self.annotator_wrapper.parallel_annotation(num_processes=n_processes, verbose=verbose)

    def quick_annotate(self, n_processes=8):
        """ Wrapper around parallel annotation multiprocess runner using MyVariant solely """
        self.annotator_wrapper.quick_annotate_and_save(n_processes=n_processes)

    def write_output_files_by_sample(self):
        """ Wrapper around function that implemts the writing of csv files for each sample in collection """
        self.annotator_wrapper.generate_output_files_by_sample()

    # TODO: Seems like this functionality (and the definition of supported build versions) really belongs in
    # AnnovarWrapper?
    def _validate_genome_version(self, build_ver):
        """ Make sure genome version is acceptable """

        if not build_ver:
            self.genome_build_version = 'hg19'  # Default genome build version

        if build_ver not in self.supported_build_vers:
            str_of_acceptable_versions = ", ".join(self.supported_build_vers)
            raise ValueError('Build version must not recognized. Supported builds are {0}'.format(
                str_of_acceptable_versions))

        else:
            self.genome_build_version = build_ver

        return self.genome_build_version
