from __future__ import division, print_function

import myvariant
import os
import sys
from pymongo import MongoClient
from base import AnnotationProject
from models import TxtParser, HgvsParser
from multiprocessing import Pool
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass

def _variant_parsing_parallel(varparse):
    varparse._variant_parsing()


class VariantParsing(AnnotationProject):

    def __init__(self,
                 input_dir,
                 output_csv_path,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None):

        super(AnnotationProject, self).__init__(input_dir,
                                                output_csv_path,
                                                annovar_path,
                                                project_data,
                                                design_file=design_file,
                                                build_ver=build_ver)

        self.chunksize = 950
        self.step = 0
        # self.csvs, self.vcfs = self.get_file_names()
        self.collection = project_data['project_name']
        self.db = project_data['db_name']
        self._buffer_len = 40000
        self._last_round = False
        self.completed_jobs = {}
        self.verbose = 0

    def annotate_and_save(self, buffer=False):

        for csv, vcf in list(zip(self.csvs, self.vcfs)):

            csv_parsing = TxtParser(csv)
            hgvs = HgvsParser(vcf)
            sample_id = os.path.splitext(os.path.basename(vcf))[0]

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

                if buffer:
                    variant_buffer = []
                    while csv_parsing.num_lines > self.step * self.chunksize:

                        list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
                        myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, sample_id)
                        offset = len(list_hgvs_ids) - self.chunksize
                        csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver, offset=offset)

                        merged_list = []
                        for i, _ in enumerate(myvariants_variants):
                            merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

                        variant_buffer.extend(merged_list)
                        self.step += 1

                        if len(merged_list) < self.chunksize:
                            self._last_round = True

                        if (len(variant_buffer) > self._buffer_len) or self._last_round:
                            logging.info('Parsing Buffer...')
                            self.export(variant_buffer)
                            variant_buffer = []

                            if self._last_round:
                                return 'Done2'

                    return 'Done1'

                else:

                    while csv_parsing.num_lines > self.step*self.chunksize:

                        list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
                        myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, sample_id)
                        offset = len(list_hgvs_ids) - self.chunksize
                        csv_variants = csv_parsing.open_and_parse_chunks(self.step, offset=offset)

                        merged_list = []
                        for i, _ in enumerate(myvariants_variants):
                            merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

                        if len(merged_list) < self.chunksize:
                            self._last_round = True

                        if self._last_round:
                            return 'Done'
                        else:
                            self.export(merged_list)
                            self.step += 1

                return 'Done'

    def check_ver(self, build_ver):
        """ Checking user input validity for genome build version """

        if not build_ver:
            return 'hg19'  # Default genome build vesion

        if build_ver not in self.supported_build_vers:
            raise ValueError('Build version must not recognized. Supported builds are'
                             ' %s, %s, %s' % (self.supported_build_vers[0],
                                              self.supported_build_vers[1],
                                              self.supported_build_vers[2]))
        else:
            return build_ver

    def export(self, list_docs):
        """
        Export data do a MongoDB server
        :param list_docs: list of dictionaries containing variant information
        :return: null
        """
        client = MongoClient()
        db = getattr(client, self.db)
        collection = getattr(db, self.collection)
        collection.insert_many(list_docs, ordered=False)

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
        get_variants_from_vcf
        from the class VariantParsing.
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

    def get_sample_csv_vcf_tuple(self, sample):
        """ Locate files associated with a specific sample """

        list_tupls = []
        vcfs = [i for i in os.listdir(os.path.join(self.input_dir, sample)) if i.endswith('vcf')]
        for vcf in vcfs:
            base_name = os.path.splitext(vcf)[0]
            matching_csv = [i for i in os.listdir(os.path.join(self.output_csv_path, sample)) if
                            i.startswith(base_name + '_annotated') and i.endswith('txt')]
            if len(matching_csv) > 1:
                raise ValueError('Too many matching csvs')
            if len(matching_csv) == 0:
                raise ValueError('Csv not found')
            else:
                csv_path = os.path.join(self.output_csv_path, sample, matching_csv[0])
                vcf_path = os.path.join(self.input_dir, sample, vcf)
                list_tupls.append((sample, vcf_path, csv_path))

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
        csv_parsing = TxtParser(maps[2])

        variant_buffer = []
        n_vars = 0

        while csv_parsing.num_lines > self.step * self.chunksize:

            list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
            myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, self.verbose, maps[0])

            offset = len(list_hgvs_ids) - self.chunksize
            csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver, offset=offset)

            merged_list = []
            for i, _ in enumerate(myvariants_variants):
                merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

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





