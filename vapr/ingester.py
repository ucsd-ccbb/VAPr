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


class Ingester:

    def __init__(self, input_dir, output_dir):

        self.base_dir = input_dir
        self.out_dir = output_dir
        #self.is_base_dir = self.check_contents()
        self.mapping_list = []
        self.atomic_ingester = SingleInstance

    @staticmethod
    def listdir_fullpath(d):
        return [os.path.join(d, f) for f in os.listdir(d)]

    def check_contents(self):
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
        walker = os.walk(self.base_dir)
        for folder, subfolders, files in walker:
            for _file in files:
                print(_file)
                full_path_single_file = os.path.join(os.path.abspath(folder), _file)
                self.digest_single(full_path_single_file)

    def digest_single(self, single):
        ingested = self.atomic_ingester(single, self.base_dir, self.out_dir)
        self.mapping_list.append(ingested.mapping)
        self.mapping_list = list({v['raw_vcf_file_full_path']: v for v in self.mapping_list}.values())

    def digest_design_file(self, design_df):

        design_file_mapping = design_df.set_index('Sample_Names').T.to_dict()
        for sample in design_file_mapping.keys():
            if not os.path.exists(os.path.join(self.base_dir, sample)):
                raise NameError('Could not find directory named %s as provided in design file' % sample)

            vcf_files = [i for i in os.listdir(os.path.join(self.base_dir, sample)) if i.endswith('.vcf')]
            logging.info('Found %i unique vcf files for sample %s' % (len(set(vcf_files)), sample))

            design_file_mapping[sample]['vcf_csv'] = [(vcf_file, os.path.splitext(os.path.basename(vcf_file))[0] +
                                                       '_annotated') for vcf_file in vcf_files]  # God bless listcomps

        return design_file_mapping


class SingleInstance:

    def __init__(self, single_instance, input_dir, out_dir):

        self.instance = single_instance
        self.out_dir = out_dir
        self.input_dir = input_dir
        self.reader = vcf.Reader(open(single_instance, 'r'))

        self.mapping = {'raw_vcf_file_full_path': self.fill_input_vcf_path(),
                        'vcf_file_basename': self.fill_vcf_file_basename(),
                        'csv_file_basename': self.fill_csv_file_basename(),
                        'sample_names': self.fill_sample_names(),
                        'num_samples_in_csv': self.fill_num_samples(),
                        'csv_file_full_path': self.fill_csv_sample_dir(),
                        'vcf_sample_dir': self.fill_vcf_sample_dir()}

        self.create_path()
        self.move_file()

    def add_key(self, key, value):
        self.mapping[key] = value

    def fill_num_samples(self):
        return len(self.reader.samples)

    def fill_sample_names(self):
        return self.reader.samples

    def sample_dir_name(self):
        return '_'.join(['sample'] + self.fill_sample_names())

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
        try:
            os.mkdir(self.mapping['csv_file_full_path'])
        except OSError:
            logging.info('Csv output dir %s for sample exists, moving files there' % self.mapping['csv_file_full_path'])
        try:
            os.mkdir(self.mapping['vcf_sample_dir'])
        except OSError:
            logging.info('Vcf input dir %s for sample exists, moving files there' % self.mapping['vcf_sample_dir'])

    def move_file(self):
        try:
            os.rename(self.mapping['raw_vcf_file_full_path'], os.path.join(self.mapping['vcf_sample_dir'],
                                                                           self.mapping['vcf_file_basename']))

            self.mapping['raw_vcf_file_full_path'] = os.path.join(self.mapping['vcf_sample_dir'],
                                                                  self.mapping['vcf_file_basename'])
        except OSError:
            logging.info('Files are organized already')

        print(os.path.isfile(self.mapping['raw_vcf_file_full_path']), self.mapping['raw_vcf_file_full_path'])