import os
import logging
import vcf
import sys
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'


class Ingester:

    def __init__(self, input_dir, output_dir):

        self.base_dir = input_dir
        self.out_dir = output_dir
        self.mapping_list = []
        self.atomic_ingester = SingleInstance

    @staticmethod
    def listdir_fullpath(d):
        return [os.path.join(d, f) for f in os.listdir(d)]

    def check_contents(self):
        """ Basic checker for input files """
        types = {'files': [], 'dirs': [], 'other': 0}
        contents = self.listdir_fullpath(self.base_dir)

        for path in contents:
            if os.path.isfile(path):
                types['files'].append(('file', path))
            elif os.path.isdir(path):
                types['dirs'].append(('dir', path))
            else:
                types['other'] += 1

        if len(set(types)) > 1:
            raise IOError('Files and dir')

    def walk(self):
        """ Walks the files, subdirs and dirs in input directory"""
        walker = os.walk(self.base_dir)
        for folder, subfolders, files in walker:
            for _file in files:
                full_path_single_file = os.path.join(os.path.abspath(folder), _file)
                self.digest_single(full_path_single_file)

    def digest_single(self, single, sample='infer', ingester_type='files', extra_data=None):
        """ Digests input data from input directory """

        ingested = self.atomic_ingester(single, self.base_dir, self.out_dir, sample=sample, ingester_type=ingester_type,
                                        extra_data=extra_data)
        ingested.mapping['extra_data'] = extra_data
        self.mapping_list.append(ingested.mapping)
        self.mapping_list = list({v['raw_vcf_file_full_path']: v for v in self.mapping_list}.values())

    def digest_design_file(self, design_df):
        """ Digests input data from csv file """

        vcf_files = []
        design_file_mapping = design_df.set_index('Sample_Names').T.to_dict()
        for sample in design_file_mapping.keys():
            if sample.endswith('.vcf'):
                sample_dir = self.base_dir
                print(os.listdir(sample_dir))
                vcf_file = [i for i in os.listdir(sample_dir) if i.startswith(sample)]
                print(vcf_file)
                if len(vcf_file) > 1:
                    raise NameError('More tha one file found with unique name in design file')
                else:
                    logging.info('Found %i unique vcf files for sample %s' % (len(set(vcf_files)), sample))
                    full_path_single_file = os.path.join(os.path.abspath(sample_dir), vcf_file[0])
                    self.digest_single(full_path_single_file, sample=sample, ingester_type='files',
                                       extra_data=design_file_mapping[sample])
            else:
                sample_dir = os.path.join(self.base_dir, sample)
                sample_dir_preprocessed = os.path.join(self.base_dir, 'sample_' + sample)
                if not os.path.exists(sample_dir) or os.listdir(sample_dir) == []:
                    if not os.path.exists(sample_dir_preprocessed):
                        raise NameError('Could not find directory named %s as provided in design file' % sample)
                    else:
                        sample_dir = sample_dir_preprocessed

                vcf_files = [i for i in os.listdir(sample_dir) if i.endswith('.vcf')]

                for vcf_file in vcf_files:
                    full_path_single_file = os.path.join(os.path.abspath(sample_dir), vcf_file)
                    self.digest_single(full_path_single_file, sample=sample, ingester_type='dirs',
                                       extra_data=design_file_mapping[sample])


class SingleInstance:

    """ Process data of single vcf file, populates dictionary """

    def __init__(self, single_instance, input_dir, out_dir, sample='infer', ingester_type='files', extra_data=None):

        self.instance = single_instance
        self.out_dir = out_dir
        self.input_dir = input_dir
        self.reader = vcf.Reader(open(single_instance, 'r'))
        self.sample = sample
        self.ingester_type = ingester_type
        self.extra_data = extra_data
        self.mapping = {'raw_vcf_file_full_path': self.fill_input_vcf_path(),
                        'vcf_file_basename': self.fill_vcf_file_basename(),
                        'csv_file_basename': self.fill_csv_file_basename(),
                        'sample_names': self.fill_sample_names(),
                        'num_samples_in_csv': len(self.fill_sample_names()),
                        'csv_file_full_path': self.fill_csv_sample_dir(),
                        'vcf_sample_dir': self.fill_vcf_sample_dir()}

        self.create_path()
        # self.move_file()

    def add_key(self, key, value):
        self.mapping[key] = value

    def fill_sample_names(self):

        if self.sample == 'infer':
            return self.reader.samples
        else:
            return [self.sample]

    def sample_dir_name(self):
        if self.ingester_type == 'dirs':
            return self.sample
        else:
            return ''

    def fill_input_vcf_path(self):
        return os.path.abspath(self.instance)

    def fill_vcf_file_basename(self):
        return os.path.basename(self.instance)

    def fill_csv_file_basename(self):
        return os.path.splitext(os.path.basename(self.instance))[0] + '_annotated'

    def fill_csv_sample_dir(self):
        return os.path.join(self.out_dir, self.sample_dir_name())

    def fill_vcf_sample_dir(self):
        return os.path.join(self.input_dir, self.sample_dir_name())

    def create_path(self):
        """ Creates directory named after the sample in a vcf file """
        try:
            os.mkdir(self.mapping['csv_file_full_path'])
        except OSError:
            logging.info('Csv output dir %s for sample exists, moving files there' % self.mapping['csv_file_full_path'])

    def move_file(self):
        """ Moves files to newly created directory """
        try:
            os.rename(self.mapping['raw_vcf_file_full_path'], os.path.join(self.mapping['vcf_sample_dir'],
                                                                           self.mapping['vcf_file_basename']))

            self.mapping['raw_vcf_file_full_path'] = os.path.join(self.mapping['vcf_sample_dir'],
                                                                  self.mapping['vcf_file_basename'])
        except OSError:
            logging.info('Files are organized already')

        print(os.path.isfile(self.mapping['raw_vcf_file_full_path']), self.mapping['raw_vcf_file_full_path'])