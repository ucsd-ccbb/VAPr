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

    def __init__(self, single_input_file_path, input_dir, out_dir, extra_data=None):

        self.single_vcf_file_path = single_input_file_path
        self.output_dir = out_dir
        self.input_dir = input_dir
        # TODO: does this open call need to be closed later?
        self.reader = vcf.Reader(open(single_input_file_path, 'r'))
        self.extra_data = extra_data

        # TODO: Why store a bunch of this info in properties of the object instance but ALSO store it in a dictionary?
        # Increases chance for bugs--if data changed in one place but not another, data in object will be inconsistent
        # and output will depend on whether data was accessed through property or dictionary.
        self.vcf_mapping_dict = {'raw_vcf_file_full_path': os.path.abspath(self.single_vcf_file_path),
                                 'vcf_file_basename': os.path.basename(self.single_vcf_file_path),
                                 'csv_file_basename': self._fill_csv_file_basename(),
                                 'sample_names': self._fill_sample_names(),
                                 'num_samples_in_csv': len(self._fill_sample_names()),
                                 'csv_file_full_path': self.output_dir,
                                 'vcf_sample_dir': self.input_dir}

        self._add_extra_data()

    def _add_extra_data(self):
        if self.extra_data:
            self.vcf_mapping_dict['extra_data'] = self.extra_data
        else:
            self.vcf_mapping_dict['extra_data'] = None

    # TODO: Refactor string values of sample into symbolic constants

    def _fill_sample_names(self):
        return self.reader.samples

    # TODO: Refactor string in file name into symbolic constant
    def _fill_csv_file_basename(self):
        return os.path.splitext(os.path.basename(self.single_vcf_file_path))[0] + '_annotated'
