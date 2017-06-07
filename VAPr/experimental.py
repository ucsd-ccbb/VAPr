


#### DEPRECATED #####

import VAPr.definitions as definitions
from VAPr.models import TxtParser, HgvsParser
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
        self.collection = project_data['project_name']
        self.db = project_data['db_name']
        self._last_round = False
        # self.completed_jobs = dict.fromkeys(list(self.mapping.keys()), 0)
        self.verbose = 0
        self.mongo_client = self.mongo_client_setup()
        self.n_vars = 0


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
    def remove_id_key(variant_data, sample_id):
        for dic in variant_data:
            dic['hgvs_id'] = dic.pop("_id", None)
            dic['hgvs_id'] = dic.pop("query", None)
            dic['sample_id'] = sample_id

        return variant_data


    def annotate_and_saving(self, buffer_vars=False, verbose=2):
        self.verbose = verbose

        list_tupls = self.get_sample_csv_vcf_tuple()
        print(list_tupls)

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

                if self._last_round:  # or (len(variant_buffer) > self._buffer_len):
                    logging.info('Done parsing variants for file pair %s, %s' % (tpl[1], tpl[2]))
                    return 'Done'

        return 'Done'



    def parallel_annotation_efficient(self, n_processes, verbose=1):
        self.verbose = verbose
        samples = self.mapping.keys()
        process_mapping = dict.fromkeys(['process_%i' % i for i in range(n_processes)], 0)

        for sample in samples:
            list_tupls = self.get_sample_csv_vcf_tuple(sample)
            for tpl in list_tupls:
                hgvs = HgvsParser(tpl[1])
                csv_parsing = TxtParser(tpl[2])

            self.pooling(n_processes, list_tupls)
            logger.info('Completed annotation and parsing for variants in sample %s' % sample)


    def _variant_parsing_efficitent(self, maps):
        hgvs = HgvsParser(maps[1])
        csv_parsing = TxtParser(maps[2])

        variant_buffer = []
        n_vars = 0

        # for i in n_processes:
        #    chunk_n =

        while csv_parsing.num_lines > self.step * self.chunksize:

            list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
            myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, self.verbose, maps[0])

            offset = len(list_hgvs_ids) - self.chunksize
            csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver, offset=offset)

            merged_list = []
            for i, _ in enumerate(myvariants_variants):
                merged_list.append(
                    self.merge_dict_lists(myvariants_variants[i], csv_variants[i][j]) for j in csv_variants[i])

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