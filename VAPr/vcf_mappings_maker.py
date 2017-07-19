import os
import logging
import vcf
import sys

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'

# TODO: is this really necessary?
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class VcfMappingsMaker:
    def __init__(self, input_dir, output_dir):
        self.base_dir = input_dir
        self.out_dir = output_dir

        self.list_of_vcf_mapping_dicts = []
        self.single_vcf_mapping_maker = SingleVcfFileMappingMaker

    # TODO: It appears this method is only ever used by check_contents, which is itself never used
    # @staticmethod
    # def listdir_fullpath(d):
    #     return [os.path.join(d, f) for f in os.listdir(d)]

    # TODO: It does not appear this code is ever used!
    # def check_contents(self):
    #     """ Basic checker for input files """
    #     types = {'files': [], 'dirs': [], 'other': 0}
    #     contents = self.listdir_fullpath(self.base_dir)
    #
    #     # TODO: Given that contents are going into labeled lists, why do they need to be a label+path tuple?
    #     for path in contents:
    #         if os.path.isfile(path):
    #             types['files'].append(('file', path))
    #         elif os.path.isdir(path):
    #             types['dirs'].append(('dir', path))
    #         else:
    #             types['other'] += 1
    #
    #     if len(set(types)) > 1:
    #         raise IOError("Detected both files and directories at path '{0}'".format(self.base_dir))

    def get_mappings_from_directory(self):
        """ Creates list of mapping dicts for VCF files found by walking the files in the base directory"""
        walker = os.walk(self.base_dir)
        for folder, _, files in walker:
            for curr_file in files:
                full_path_single_file = os.path.join(os.path.abspath(folder), curr_file)
                self._store_mapping_for_single_vcf(full_path_single_file)

    def get_mappings_from_design_file(self, design_df):
        """ Creates list of mapping dicts for VCF files referenced in dataframe containing contents of design file """

        VCF_EXTENSION = ".vcf"

        # TODO: Find out required format for a design file; apparently it has to have a column named Sample_Names ...
        # what else?
        # TODO: Figure out what this line is doing
        design_file_mapping = design_df.set_index('Sample_Names').T.to_dict()

        # TODO: I personally think we shouldn't support a design file that mixes the two, but
        # the code as it is allows that without error or comment
        for sample_identifier in design_file_mapping.keys():
            if sample_identifier.endswith(VCF_EXTENSION):  # Design file contains file names
                sample_dir = self.base_dir
                vcf_file = [i for i in os.listdir(sample_dir) if i.startswith(sample_identifier)]
                if len(vcf_file) > 1:
                    raise NameError("More than one vcf file found that starts with sample identifier '{0}' from the "
                                    "design file".format(sample_identifier))
                else:
                    logging.info('Found %i unique vcf files for sample %s' % (len(set(vcf_file)), sample_identifier))
                    full_path_single_file = os.path.join(os.path.abspath(sample_dir), vcf_file[0])
                    self._store_mapping_for_single_vcf(full_path_single_file, sample_id=sample_identifier,
                                                       sample_id_type='files',
                                                       extra_data=design_file_mapping[sample_identifier])

            else:  # Design file contains directory (sample) names
                sample_dir = os.path.join(self.base_dir, sample_identifier)
                if not os.path.exists(sample_dir) or os.listdir(sample_dir) == []:
                    raise NameError('Could not find directory named %s as provided in design file' % sample_identifier)

                vcf_files = [i for i in os.listdir(sample_dir) if i.endswith(VCF_EXTENSION)]
                for vcf_file in vcf_files:
                    full_path_single_file = os.path.join(os.path.abspath(sample_dir), vcf_file)
                    self._store_mapping_for_single_vcf(full_path_single_file, sample_id=sample_identifier,
                                                       sample_id_type='dirs',
                                                       extra_data=design_file_mapping[sample_identifier])

    def _store_mapping_for_single_vcf(self, single_input_file_path, sample_id='infer', sample_id_type='files',
                                      extra_data=None):
        """ Digests input data from input directory """

        vcf_mapping_maker = self.single_vcf_mapping_maker(single_input_file_path, self.base_dir, self.out_dir,
                                                          sample_id=sample_id, sample_id_type=sample_id_type,
                                                          extra_data=extra_data)

        # reach into mapping maker above to get mapping for the single input vcf file, add that mapping to the list
        # of mappings for all vcf files
        self.list_of_vcf_mapping_dicts.append(vcf_mapping_maker.vcf_mapping_dict)

        # Eliminates possible duplicates--by creating a dictionary in which each value is an original dictionary
        # from self.mapping_list and each key is the value of 'raw_vcf_file_full_path' in that dictionary,
        # *then* pulling out only the values (which is to say, original dictionaries) from that new dictionary.
        # This has the effect of ensuring that there are no duplicates in the value of 'raw_vcf_file_full_path',
        # although if two dictionaries have the same 'raw_vcf_file_full_path' value but different other values,
        # this approach will pick one arbitrarily (and probably not stably across runs).
        # TODO: Figure out whether this approach really accomplishes what we want here
        # TODO: Pull out string keys into symbolic constants
        self.list_of_vcf_mapping_dicts = list(
            {v['raw_vcf_file_full_path']: v for v in self.list_of_vcf_mapping_dicts}.values()
        )


# TODO: Again, not sure that class-based approach is the best way to go here, as we don't really need this object to
# stick around and hold state over time--we just want a one-shot, get-it-and-go approach.  Consider refactoring as
# top-level function in its own module, supported by private functions as necessary.
class SingleVcfFileMappingMaker:
    """ Populate mapping dictionary for single VCF file """

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
        self.vcf_mapping_dict = {'raw_vcf_file_full_path': os.path.abspath(self.single_vcf_file_path),
                        'vcf_file_basename': os.path.basename(self.single_vcf_file_path),
                        'csv_file_basename': self._fill_csv_file_basename(),
                        'sample_names': self._fill_sample_names(),
                        'num_samples_in_csv': len(self._fill_sample_names()),
                        'csv_file_full_path': os.path.join(self.out_dir, self._sample_dir_name()),
                        'vcf_sample_dir': os.path.join(self.input_dir, self._sample_dir_name())}

        self._add_extra_data()
        self._create_csv_output_dir()

    def _add_extra_data(self):
        if self.extra_data:
            self.vcf_mapping_dict['extra_data'] = self.extra_data
        else:
            self.vcf_mapping_dict['extra_data'] = None

    # TODO: It appears this method is never used
    # def add_key(self, key, value):
    #     self.mapping[key] = value

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

    # TODO: Refactor string key into symbolic constant
    def _create_csv_output_dir(self):
        """ Creates directory named after the sample in a vcf file """
        try:
            os.mkdir(self.vcf_mapping_dict['csv_file_full_path'])
        except OSError:
            # catch error, don't propagate
            logging.info('Csv output directory %s for sample already exists; using existing directory' %
                         self.vcf_mapping_dict['csv_file_full_path'])

    # TODO: It appears this method is never used
    # def move_file(self):
    #     """ Moves files to newly created directory """
    #     try:
    #         os.rename(self.mapping['raw_vcf_file_full_path'], os.path.join(self.mapping['vcf_sample_dir'],
    #                                                                        self.mapping['vcf_file_basename']))
    #
    #         self.mapping['raw_vcf_file_full_path'] = os.path.join(self.mapping['vcf_sample_dir'],
    #                                                               self.mapping['vcf_file_basename'])
    #     except OSError:
    #         logging.info('Files are organized already')
    #
    #     print(os.path.isfile(self.mapping['raw_vcf_file_full_path']), self.mapping['raw_vcf_file_full_path'])