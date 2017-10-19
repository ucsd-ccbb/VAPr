from __future__ import division, print_function
import sys
import pandas
import logging
from VAPr.annovar import AnnovarWrapper
from VAPr.parsers import VariantParsing
from VAPr.vcf_mappings_maker import VcfMappingsMaker

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
    supported_build_vers = ['hg19', 'hg18', 'hg38']

    def __init__(self, input_dir, output_dir, annovar_path, mongo_db_and_collection_names_dict,
                 design_file=None, build_ver=None, mongod_cmd=None, split_vcf=False):
        """ Class that implements the API and the major annotation/saving methods  """

        self.input_dir = input_dir
        self.output_csv_path = output_dir
        self.design_file = design_file
        self.annovar_path = annovar_path
        self.mongo_db_and_collection_names_dict = mongo_db_and_collection_names_dict
        self.genome_build_version = self._validate_genome_version(build_ver)
        self.times_called = 0
        self.list_of_vcf_mapping_dicts = self._get_list_of_vcf_mapping_dicts()
        self.split = split_vcf
        self.mongod = mongod_cmd

        # TODO: These two calls takes in exactly the same parameters; Consider storing them to an object and passing
        # object around?
        self.annovar_wrapper = AnnovarWrapper(self.input_dir, self.output_csv_path, self.annovar_path,
                                              self.mongo_db_and_collection_names_dict, self.list_of_vcf_mapping_dicts,
                                              design_file=self.design_file,
                                              genome_build_version=self.genome_build_version)

        self.annotator_wrapper = VariantParsing(self.input_dir, self.output_csv_path, self.annovar_path,
                                                self.mongo_db_and_collection_names_dict, self.list_of_vcf_mapping_dicts,
                                                design_file=self.design_file,
                                                build_ver=self.genome_build_version,
                                                mongod_cmd=self.mongod)

    def download_dbs(self, dbs=None):
        """ Wrapper around Annovar database downloading function """
        self.annovar_wrapper.download_dbs(annovar_dbs_to_get=dbs)

    def run_annovar(self, batch_jobs=10, multisample=False):
        """ Wrapper around multiprocess Annovar annotation  """
        self.annovar_wrapper.run_annovar(num_batch_jobs=batch_jobs, vcf_is_multisample=multisample)

    # TODO: Apparent bug
    # I cannot find any definition of a method called annotate_and_saving anywhere in the project, so
    # I can't see how this method (the key method of the entire codebase as described in the documentation)
    # could possibly run.  Perhaps the author's intention was that it would be superceded by a choice of either
    # parallel_annotation_and_saving or quick_annotate_and_save?
    def annotate_and_save(self, buffer_vars=False, verbose=2):
        """ Wrapper around annotation runner. Deprecated in favour of parallel parsing """
        self.annotator_wrapper.annotate_and_saving(buffer_vars=buffer_vars, verbose=verbose)

    def parallel_annotation_and_saving(self, n_processes=4, verbose=1):
        """ Wrapper around parallel annotation multiprocess runner  """
        self.annotator_wrapper.parallel_annotation(num_processes=n_processes, verbose=verbose)

    def quick_annotate_and_save(self, n_processes=8):
        """ Wrapper around parallel annotation multiprocess runner using MyVariant solely """
        self.annotator_wrapper.quick_annotate_and_save(n_processes=n_processes)

    def write_output_files_by_sample(self):
        """ Wrapper around function that implemts the writing of csv files for each sample in collection """
        self.annotator_wrapper.generate_output_files_by_sample()

    def _get_list_of_vcf_mapping_dicts(self):
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
            list_of_vcf_mapping_dicts = self._get_vcf_mappings_from_design_file(design_df)
        else:
            list_of_vcf_mapping_dicts = self._get_vcf_mappings_from_directory()

        return list_of_vcf_mapping_dicts

    # TODO: It appears this method is never used
    # @staticmethod
    # def listdir_fullpath(d):
    #     """ Helper function to list full path of files in directory """
    #     return [os.path.join(d, f) for f in os.listdir(d)]

    def _get_vcf_mappings_from_directory(self):
        """ Ingest all files in specified directory and return mapping """

        # TODO: This is not the most suitable pattern for this task--given that python allows top-level functions,
        # having an object that must be instantiated and then dereferenced for its method(s) and propertie(s) is
        # unnecessary here.  Simply have a method in the ingester module called get_file_mapping_from_directory that
        # takes in the input dir and output csv path, and returns the mapping. (The VcfMappingsMaker object isn't persisted
        # anyway, so its saved state must be irrelevant.)  Added advantage is that this method (and one below) go away
        # altogether, replaced by a call (as to ingester.get_file_mapping_from_directory) in the above
        # _get_mapping_of_vcf_file method, which contains some actual logic.
        organizer = VcfMappingsMaker(self.input_dir, self.output_csv_path)
        organizer.get_mappings_from_directory()
        return organizer.list_of_vcf_mapping_dicts

    def _get_vcf_mappings_from_design_file(self, design_df):
        """ Ingest design file and the directories/files referenced in it and return mapping """

        # TODO: See comment above in _get_vcf_mappings_from_directory
        organizer = VcfMappingsMaker(self.input_dir, self.output_csv_path)
        organizer.get_mappings_from_design_file(design_df)
        return organizer.list_of_vcf_mapping_dicts

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
