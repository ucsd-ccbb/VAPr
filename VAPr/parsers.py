from __future__ import division, print_function

import myvariant
import os
import sys
import tqdm
from multiprocessing import Pool
from multiprocessing.dummy import Pool as ThreadPool
from pymongo import MongoClient
from src.base import AnnotationProject
from src.models import TxtParser, HgvsParser
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.handlers[0].stream = sys.stdout


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
        self.csvs, self.vcfs = self.get_file_names()
        self.collection = project_data['project_name']
        self.db = project_data['db_name']
        self._buffer_len = 40000
        self._last_round = False

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

    def get_file_names(self):
        vcfs = [i[0] for i in self.mapping]
        inits = [os.path.basename(i[1]) for i in self.mapping]
        csvs = []
        all_csvs = os.listdir(self.output_csv_path)
        for csv in all_csvs:
            for csv_init in inits:
                if csv.startswith(csv_init):
                    csvs.append(os.path.join(self.output_csv_path, csv))

        return vcfs, csvs

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

    def get_dict_myvariant(self, variant_list, sample_id):
        """
        Function designated to place the queries on myvariant.info servers.

        :param variant_list: list of HGVS variant ID's. Usually retrived beforehand using the method
        get_variants_from_vcf
        from the class VariantParsing.
        :return: list of dictionaries. Each dictionary contains data about a single variant.
        """

        mv = myvariant.MyVariantInfo()
        # This will retrieve a list of dictionaries
        variant_data = mv.getvariants(variant_list, as_dataframe=False)
        variant_data = self.remove_id_key(variant_data, sample_id)

        return variant_data

    def remove_id_key(self, variant_data, sample_id):

        for dic in variant_data:
            dic['hgvs_id'] = dic.pop("_id", None)
            dic['hgvs_id'] = dic.pop("query", None)
            dic['sample_id'] = sample_id

        return variant_data

    def add_csv_to_mapping(self):

        inits = self.mapping.keys()
        all_csvs = os.listdir(self.output_csv_path)
        list_of_dicts = []
        for csv in all_csvs:
            for init in inits:
                if all([csv.startswith(init), csv.endswith('txt')]):
                    self.mapping[init]['csv'] = os.path.join(self.output_csv_path, csv)
                    # list_of_dicts.append({init: self.mapping.values()})
        for key in self.mapping.keys():
            list_of_dicts.append({key: self.mapping[key]})
        return list_of_dicts

    def parallel_annotation(self, n_processes):

        list_tupls = []
        mapping = self.add_csv_to_mapping()
        for k in mapping:
            key = list(k.keys())[0]
            list_tupls.append((key, k[key]['vcf'], k[key]['csv']))

        self.pooling(n_processes, self._variant_parsing, list_tupls)

    def _variant_parsing(self, maps):

        hgvs = HgvsParser(maps[1])
        csv_parsing = TxtParser(maps[2])

        variant_buffer = []
        while csv_parsing.num_lines > self.step * self.chunksize:
            # print('LINE EXIT COND: %i, %i' % (csv_parsing.num_lines, self.step * self.chunksize))
            # logging.info('Parsing %i/%i variants from file %s' % ((self.step + 1) * self.chunksize,
            #                                                      csv_parsing.num_lines,
            #                                                      os.path.basename(maps[1])))

            list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
            # tpl = (maps[0], list_hgvs_ids)

            myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, maps[0]) #self.threading(8, self.get_dict_mv_parallel, tpl)
            offset = len(list_hgvs_ids) - self.chunksize
            csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver, offset=offset)

            merged_list = []
            for i, _ in enumerate(myvariants_variants):
                merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

            variant_buffer.extend(merged_list)
            self.step += 1
            # print(self._last_round)

            if len(merged_list) < self.chunksize:
                self._last_round = True

            # print(self._last_round)
            if (len(variant_buffer) > self._buffer_len) or self._last_round:
                logging.info('Parsing Buffer...')
                self.export(variant_buffer)
                variant_buffer = []

                if self._last_round:
                    return 'Done2'
            # print('LINE EXIT COND: %i, %i' % (csv_parsing.num_lines, self.step * self.chunksize))

        return 'Done1'

    @staticmethod
    def threading(n_threads, input_1, input_2):

        pool = ThreadPool(n_threads)
        results = pool.map(input_1, input_2)
        pool.close()
        pool.join()
        return results

    @staticmethod
    def pooling(n_processes, input_1, input_2):

        pool = Pool(processes=n_processes)
        # rs = pool.imap(input_1, input_2)
        for _ in tqdm.tqdm(pool.imap_unordered(input_1, input_2), total=len(input_2)):
            pass
