from __future__ import division, print_function

import myvariant
import os
import sys
from pymongo import MongoClient
import vapr.definitions as definitions
from vapr.models import TxtParser, HgvsParser

from multiprocessing import Pool
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class VariantParsing:

    def __init__(self,
                 input_dir,
                 output_csv_path,
                 annovar_path,
                 project_data,
                 mapping,
                 design_file=None,
                 build_ver=None):

        """ Project data """
        self.input_dir = input_dir
        self.output_csv_path = output_csv_path
        self.annovar = annovar_path
        self.project_data = project_data
        self.design_file = design_file
        self.buildver = build_ver
        self.mapping = mapping

        self.chunksize = definitions.chunk_size
        self.step = 0
        # self.csvs, self.vcfs = self.get_file_names()
        self.collection = project_data['project_name']
        self.db = project_data['db_name']
        # self._buffer_len = 50000
        self._last_round = False
        # self.completed_jobs = dict.fromkeys(list(self.mapping.keys()), 0)
        self.verbose = 0
        self.collection = self.mongo_client_setup()

    def annotate_and_saving(self, buffer_vars=False, verbose=2):

        self.verbose = verbose

        list_tupls = self.get_sample_csv_vcf_tuple()

        for tpl in list_tupls:
            hgvs = HgvsParser(tpl[1])
            test = hgvs.get_variants_from_vcf(1)
            print(len(test))
            csv_parsing = TxtParser(tpl[2], samples=hgvs.samples)
            variant_buffer = []
            self.step = 0
            num_parsed = 0

            while csv_parsing.num_lines > self.step * self.chunksize:

                print('STEP: %i, CHNKSIZE: %i, MULT: %i' % (self.step, self.chunksize, self.step * self.chunksize))
                list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
                all_ids = hgvs.get_all_variants_from_vcf()
                print('Len List HGVC ID: %i, ALL IDS: %i' % (len(list_hgvs_ids), len(all_ids)),
                      'Num Lines: %i' % csv_parsing.num_lines)

                myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, verbose, tpl[0])
                num_parsed += len(list_hgvs_ids)

                offset = len(list_hgvs_ids) - self.chunksize
                csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver,
                                                                 offset=offset)

                print('CSV_VARIANTS: %i' % len(csv_variants))
                if self.verbose >= 1:
                    logger.info('Gathered %i variants so far for sample %s, vcf file %s' % (num_parsed,
                                                                                            tpl[0],
                                                                                            tpl[1]))
                merged_list = []
                for i, _ in enumerate(myvariants_variants):
                    for dict_from_sample in csv_variants[i]:
                        merged_list.append(self.merge_dict_lists(myvariants_variants[i], dict_from_sample))

                variant_buffer.extend(merged_list)
                print('VARBUFFLEN: %i' % len(variant_buffer))
                logging.info('Parsing Buffer...')
                if len(variant_buffer) == 0:
                    print('VarBuff Len == 0, no clue why')
                else:
                    self.collection.insert_many(variant_buffer, ordered=False)
                variant_buffer = []
                self.step += 1

                if len(list_hgvs_ids) < self.chunksize:
                    self._last_round = True

                if self._last_round:    # or (len(variant_buffer) > self._buffer_len):
                    logging.info('Done parsing variants for file pair %s, %s' % (tpl[1], tpl[2]))
                    return 'Done'

        return 'Done'
    """
            if not csv:

                while hgvs.num_lines > self.step * self.chunksize:

                    list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
                    myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, sample_id)

                    if len(myvariants_variants) < self.chunksize:
                        self._last_round = True

                    if self._last_round:
                        return 'Done'
                    else:
                        self.export(myvariants_variants)
                        self.step += 1

                return 'Done'

            else:
    """

    def mongo_client_setup(self):
        """
        Setup MongoDB client
        :return: null
        """
        client = MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)
        db = getattr(client, self.db)
        collection = getattr(db, self.collection)

        return collection

    @staticmethod
    def merge_dict_lists(*dict_args):
        """
        Given any number of dicts, shallow copy and merge into a new dict,
        precedence goes to key value pairs in latter dicts.
        """
        result = {}
        for dictionary in dict_args:
            result.update(dictionary)
        return result

    def get_dict_mv_parallel(self, mapping):
        """
        Function designated to place the queries on myvariant.info servers.
        """
        sample_id = mapping[0]
        ids = mapping[1]
        mv = myvariant.MyVariantInfo()
        # This will retrieve a list of dictionaries
        variant_data = mv.getvariants(ids, as_dataframe=False)
        variant_data = self.remove_id_key(variant_data, sample_id)

        return variant_data

    def get_dict_myvariant(self, variant_list, verbose, sample_id):
        """
        Function designated to place the queries on myvariant.info servers.

        :param variant_list: list of HGVS variant ID's. Usually retrived beforehand using the method
        get_variants_from_vcf from the class VariantParsing.
        :param verbose:
        :return: list of dictionaries. Each dictionary contains data about a single variant.
        """

        if verbose >= 2:
            verbose = True
        else:
            verbose = False

        mv = myvariant.MyVariantInfo()
        # This will retrieve a list of dictionaries
        variant_data = mv.getvariants(variant_list, verbose=verbose, as_dataframe=False)
        variant_data = self.remove_id_key(variant_data, sample_id)

        return variant_data

    @staticmethod
    def remove_id_key(variant_data, sample_id):

        for dic in variant_data:
            dic['hgvs_id'] = dic.pop("_id", None)
            dic['hgvs_id'] = dic.pop("query", None)
            dic['sample_id'] = sample_id

        return variant_data

    def get_sample_csv_vcf_tuple(self):
        """ Locate files associated with a specific sample """
        list_tupls = []

        for _map in self.mapping:
            base_name = os.path.splitext(os.path.basename(_map['vcf_file']))[0]

            matching_csv = [i for i in os.listdir(_map['annotation_dir'])
                            if i.startswith(base_name + '_annotated') and i.endswith('txt')]

            matching_vcf = [i for i in os.listdir(_map['annotation_dir'])
                            if i.startswith(base_name + '_annotated') and i.endswith('vcf')]

            print(base_name, matching_vcf, matching_csv)

            if len(matching_csv) > 1 or len(matching_vcf) > 1:
                raise ValueError('Too many matching csvs')
            elif len(matching_csv) == 0 or len(matching_vcf) == 0:
                raise ValueError('Csv not found')
            else:
                csv_path = os.path.join(_map['annotation_dir'], matching_csv[0])
                vcf_path = os.path.join(_map['annotation_dir'], matching_vcf[0])
                list_tupls.append((_map['samples'], vcf_path, csv_path))

        return list_tupls

    def parallel_annotation(self, n_processes, verbose=1):

        """
        Set up variant parsing scheme. Since a functional programming style is required for parallel
        processing, the input to the Pool.map function from the multiprocessing library must be
        immutable. We use tuples of the type:

            (samples, vcf_path, csv_path)

        :param n_processes: number of cores to be used
        :param verbose: verbosity level [0,1,2,3]
        :return: None
        """

        self.verbose = verbose
        samples = self.mapping.keys()

        for sample in samples:
            list_tupls = self.get_sample_csv_vcf_tuple(sample)
            self.pooling(n_processes, list_tupls)
            logger.info('Completed annotation and parsing for variants in sample %s' % sample)

    def _variant_parsing(self, maps):

        """
        Given a list of tuples, parse the hgvs IDS using HGVS class, and parse csv data through the
        TxtParser classes. Query variants from MyVarinant.info and aggregate results in a dictionary,
        to be exported to MongoDB

        """

        hgvs = HgvsParser(maps[1])
        csv_parsing = TxtParser(maps[2], samples=hgvs.samples)

        variant_buffer = []
        n_vars = 0

        while csv_parsing.num_lines > self.step * self.chunksize:

            list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
            myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, self.verbose, maps[0])

            offset = len(list_hgvs_ids) - self.chunksize
            csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver, offset=offset)

            merged_list = []
            for i, _ in enumerate(myvariants_variants):

                for dict_from_sample in csv_variants[i]:
                    merged_list.append(self.merge_dict_lists(myvariants_variants[i], dict_from_sample))

            variant_buffer.extend(merged_list)
            n_vars += len(merged_list)
            if self.verbose >= 1:
                logger.info('Gathered %i variants so far for sample %s, vcf file %s' % (n_vars, maps[0], maps[1]))
            self.step += 1

            if len(merged_list) < self.chunksize:
                self._last_round = True

            if (len(variant_buffer) > self._buffer_len) or self._last_round:
                logging.info('Parsing Buffer...')
                self.export(variant_buffer)
                variant_buffer = []

                if self._last_round:
                    self.completed_jobs[maps[0]] += 1
                    return self.completed_jobs

        return self.completed_jobs

    def __call__(self, x):
        return self._variant_parsing(x)

    def pooling(self, n_processes, input_2):
        """ Temporary hack to make instance method _variant_parsing 'pickleable' """
        Pool(n_processes).map(self, input_2)

if __name__ == '__main__':
    # Directory of input files to be annotated
    IN_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_multisample/"

    # Output file directory
    OUT_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/"
    # Location of your annovar dowload. The folder should contain the following files/directories:
    ANNOVAR_PATH = '/Volumes/Carlo_HD1/CCBB/annovar/'

    # Design File (optional)
    # design_file = '/Volumes/Carlo_HD1/CCBB/VAPr_files/guorong_single_sample.csv'

    # Databse and Collection names (optional)
    proj_data = {'db_name': 'VariantDBMultiBenchmark',
                 'project_name': 'collect'}
    mapping = [{'annotation_dir': '/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_multisample/samples_BC001_BC001_MOM_BC001_SIS',
                  'csv_file': 'BC001_family.g_annotated',
                  'path_to_dir': '/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_multisample/samples_BC001_BC001_MOM_BC001_SIS',
                  'samples': ['BC001', 'BC001_MOM', 'BC001_SIS'],
                  'vcf_file': '/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_multisample/samples_BC001_BC001_MOM_BC001_SIS/BC001_family.g.vcf'}]

    annotator =  VariantParsing(IN_PATH,
                                OUT_PATH,
                                ANNOVAR_PATH,
                                proj_data,
                                mapping)
    annotator.annotate_and_saving(verbose=2)
