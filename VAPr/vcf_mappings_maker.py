import os
import logging
import vcf
import sys
import pandas

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'

# TODO: is this really necessary?
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class SingleVcfFileMappingMaker:
    """ Populate mapping dictionary for single VCF file """

    VCF_FILE_PATH_KEY = 'raw_vcf_file_full_path'
    VCF_FILE_BASENAME_KEY = 'vcf_file_basename'
    ANNOVAR_OUTPUT_BASENAME_KEY = 'csv_file_basename'
    SAMPLE_NAMES_LIST_KEY = 'sample_names'
    NUM_SAMPLES_IN_CSV_KEY = 'num_samples_in_csv'
    ANNOVAR_OUTPUT_DIR = 'csv_file_full_path'
    VCF_SAMPLE_DIR_KEY = 'vcf_sample_dir'

    def __init__(self, single_input_file_path, input_dir, out_dir, sample_id='infer', sample_id_type='files',
                 extra_data=None):

        self.single_vcf_file_path = single_input_file_path
        self.out_dir = out_dir
        self.input_dir = input_dir
        # TODO: does this open call need to be closed later?
        self.reader = vcf.Reader(open(single_input_file_path, 'r'))
        self.sample_id = sample_id
        self.sample_id_type = sample_id_type
        self.extra_data = extra_data

        # TODO: Why store a bunch of this info in properties of the object instance but ALSO store it in a dictionary?
        # Increases chance for bugs--if data changed in one place but not another, data in object will be inconsistent
        # and output will depend on whether data was accessed through property or dictionary.
        self.vcf_mapping_dict = {self.VCF_FILE_PATH_KEY: os.path.abspath(self.single_vcf_file_path),
                                 self.VCF_FILE_BASENAME_KEY: os.path.basename(self.single_vcf_file_path),
                                 self.ANNOVAR_OUTPUT_BASENAME_KEY: self._fill_csv_file_basename(),
                                 self.SAMPLE_NAMES_LIST_KEY: self._fill_sample_names(),
                                 self.NUM_SAMPLES_IN_CSV_KEY: len(self._fill_sample_names()),
                                 self.ANNOVAR_OUTPUT_DIR: os.path.join(self.out_dir, self._sample_dir_name()),
                                 self.VCF_SAMPLE_DIR_KEY: os.path.join(self.input_dir, self._sample_dir_name())}

        self._add_extra_data()

    def _add_extra_data(self):
        if self.extra_data:
            self.vcf_mapping_dict['extra_data'] = self.extra_data
        else:
            self.vcf_mapping_dict['extra_data'] = None

    # TODO: Refactor string values of sample into symbolic constants

    def _fill_sample_names(self):
        if self.sample_id == 'infer':
            return self.reader.samples
        else:
            return [self.sample_id]

    # TODO: Refactor into property!
    def _sample_dir_name(self):
        if self.sample_id_type == 'dirs':
            return self.sample_id
        else:
            return ''

    # TODO: Refactor string in file name into symbolic constant
    def _fill_csv_file_basename(self):
        return os.path.splitext(os.path.basename(self.single_vcf_file_path))[0] + '_annotated'
