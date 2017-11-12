from __future__ import division, print_function

# built-in libraries
import itertools
import os
import sys
from multiprocessing import Pool
import logging
import subprocess
import shlex
import time
import tqdm

# third-party libraries
import myvariant
from pymongo import MongoClient
import vcf

# project libraries
# import VAPr.definitions as definitions
from VAPr.annovar_output_parsing import AnnovarTxtParser
from VAPr.writes import Writer
from VAPr.queries import Filters

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'

# TODO: Understand, vet this logging set-up
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


def _get_num_lines_in_file(file_path):
    return sum(1 for _ in open(file_path))


class VariantParsing:
    CHUNK_INDEX_INDEX = 0
    FILE_PATH_INDEX = 1
    CHUNK_SIZE_INDEX = 2
    DB_NAME_INDEX = 3
    COLLECTION_NAME_INDEX = 4
    GENOME_BUILD_VERSION_INDEX = 5
    VERBOSE_LEVEL_INDEX = 6
    SAMPLE_LIST_INDEX = 7

    # def __init__(self, input_dir, output_csv_path, annovar_path, mongo_db_and_collection_names_dict,
    #              vcf_mapping_dict, design_file=None, build_ver=None, mongod_cmd=None):
    #
    #     # Project data
    #     self.input_dir = input_dir
    #     self.output_csv_path = output_csv_path
    #     self.annovar_path = annovar_path
    #     self.project_data = mongo_db_and_collection_names_dict
    #     self.design_file = design_file
    #     self.genome_build_version = build_ver
    #     self._vcf_mapping_dict = vcf_mapping_dict
    #     self.chunksize = definitions.chunk_size
    #     self.step = 0
    #
    #     # TODO: refactor these string keys into symbolic constants
    #     self.collection = mongo_db_and_collection_names_dict['collection_name']
    #     self.db = mongo_db_and_collection_names_dict['db_name']
    #     self._last_round = False
    #     self.verbose = 0
    #     self.n_vars = 0
    #     self.mongod = mongod_cmd
    #
    # def parallel_annotation(self, num_processes, verbose=1):
    #     """
    #     Set up variant parsing scheme. Since a functional programming style is required for parallel
    #     processing, the input to the Pool.map function from the multiprocessing library must be
    #     immutable. We use tuples of the format:
    #
    #         (samples, vcf_path, csv_path, db_name, collection_name, step)
    #
    #     The parallel annotation will split the vcf, csv file in N steps where N equals the number of lines
    #     (variants) to be annotated divided a pre-defined 'chunk size'. For instance, a vcf, csv pair containing
    #     250 000 variants and a chunksize of 5 000 will yield 50 steps. Given 8 parallel procrsses, the annotation
    #     will happen over those 50 steps with 8 cores working cuncurrently until the 50 jobs are fully annotated.
    #
    #     :param num_processes: number of cores to be used
    #     :param verbose: verbosity level [0,1,2,3]
    #     :return: None
    #     """
    #
    #     self.verbose = verbose
    #     list_tuples = self._get_sample_csv_vcf_tuple()
    #     map_job = []
    #     for tpl in list_tuples:
    #         # hgvs = HgvsParser(tpl[1])
    #         num_lines = _get_num_lines_in_file(tpl[2])  # num lines in csv_path
    #         n_steps = int(num_lines/self.chunksize) + 1
    #         map_job.extend(self._make_full_tuple_list(tpl, n_steps, extra_data=tpl[5], mongod_cmd=self.mongod))
    #
    #     pool = Pool(num_processes)
    #     for _ in tqdm.tqdm(pool.imap_unordered(_parse_by_step, map_job), total=len(map_job)):
    #         pass
    #     pool.close()
    #     pool.join()
    #
    #     logger.info('Completed annotation and parsing for variants in sample %s' % tpl[0])

    # TODO: someday: extra_data from design file needs to come back in here ...
    @staticmethod
    def merge_annotations_and_save(db_name, collection_name, build_version, file_path, num_processes,
                                   chunk_size, sample=None, verbose_level=1):

        shared_job_params = [file_path,
                             chunk_size,
                             db_name,
                             collection_name,
                             build_version,
                             verbose_level]

        if sample is not None:
            shared_job_params.append(sample)

        jobs_params_tuples_list = VariantParsing._make_jobs_params_tuples_list(file_path, chunk_size, shared_job_params)

        pool = Pool(num_processes)
        for _ in tqdm.tqdm(pool.imap_unordered(_parse_by_step, jobs_params_tuples_list),
                           total=len(jobs_params_tuples_list)):
            pass
        pool.close()
        pool.join()

        # logger.info('Completed annotation and parsing for variants in sample %s' % tpl[0])

    # def quick_annotate_and_save(self, n_processes=8):
    #     """ Annotation that doesn't require annovar """
    #
    #     list_tuples = []
    #     simple_map = {k: v for k, v in self._vcf_mapping_dict.items() if k not in ['csv_file_basename', 'sample_names']}
    #     list_tuples.append((simple_map['raw_vcf_file_full_path'], self.db, self.collection))
    #
    #     for _tuple in list_tuples:
    #         # hgvs = HgvsParser(_tuple[0])
    #         num_lines = _get_num_lines_in_file(_tuple[0])
    #         n_steps = int(num_lines / self.chunksize) + 1
    #         map_job = self._make_simple_tuple_list(_tuple, n_steps)
    #         pool = Pool(n_processes)
    #         for _ in tqdm.tqdm(pool.imap_unordered(_parallel_get_dict_mv, map_job), total=len(map_job)):
    #             pass
    #         pool.close()
    #         pool.join()

    @staticmethod
    def generate_output_files_by_sample(db_name, collection_name, output_path):
        client = MongoClient()
        db = getattr(client, db_name)
        collection = getattr(db, collection_name)

        # Get all distinct sample_ids
        samples = collection.distinct('sample_id')

        fwriter = Writer(db_name, collection_name)
        filt = Filters(db_name, collection_name)

        # Generate files for each sample
        for sample in samples:
            q = filt.variants_from_sample(sample)
            list_docs = list(q)
            out_path_by_sample = os.path.join(output_path, 'csv_by_sample')
            if not os.path.exists(out_path_by_sample):
                os.makedirs(out_path_by_sample)
            fname = os.path.join(out_path_by_sample, 'annotated_csv_' + sample + '_all_vars.csv')
            fwriter.generate_annotated_csv(list_docs, fname)

            # def _get_csv_path(self):
            #     """ Locate files associated with a specific sample """
            #
            #     matching_csv = [i for i in os.listdir(self._vcf_mapping_dict['csv_file_full_path']) if
            #                     i.startswith(self._vcf_mapping_dict['csv_file_basename']) and i.endswith('txt')]
            #
            #     if len(matching_csv) > 1:
            #         raise ValueError('Too many matching csvs')
            #     elif len(matching_csv) == 0:
            #         raise ValueError('Csv not found')
            #
            #     csv_path = os.path.join(self._vcf_mapping_dict['csv_file_full_path'], matching_csv[0])
            #     return csv_path
            #
            # def _get_sample_csv_vcf_tuple(self):
            #     """ Locate files associated with a specific sample """
            #
            #     list_tuples = []
            #
            #     matching_csv = [i for i in os.listdir(self._vcf_mapping_dict['csv_file_full_path']) if
            #                     i.startswith(self._vcf_mapping_dict['csv_file_basename']) and i.endswith('txt')]
            #
            #     matching_vcf = [i for i in os.listdir(self._vcf_mapping_dict['csv_file_full_path'])
            #                     if i.startswith(self._vcf_mapping_dict['csv_file_basename']) and i.endswith('vcf')]
            #
            #     if len(matching_csv) > 1 or len(matching_vcf) > 1:
            #         raise ValueError('Too many matching csvs')
            #     elif len(matching_csv) == 0 or len(matching_vcf) == 0:
            #         raise ValueError('Csv not found')
            #     else:
            #         csv_path = os.path.join(self._vcf_mapping_dict['csv_file_full_path'], matching_csv[0])
            #         vcf_path = os.path.join(self._vcf_mapping_dict['csv_file_full_path'], matching_vcf[0])
            #
            #         list_tuples.append(
            #             (self._vcf_mapping_dict['sample_names'],
            #              vcf_path,
            #              csv_path,
            #              self.db,
            #              self.collection,
            #              self._vcf_mapping_dict['extra_data'])
            #         )
            #
            #     return list_tuples
            #
            # def _make_full_tuple_list(self, _tuple, n_steps, extra_data=None, mongod_cmd=None):
            #     """ Assign step number to each tuple to be consumed by parsing function """
            #     new_tuple_list = []
            #     if mongod_cmd:
            #         mongod = mongod_cmd
            #     else:
            #         mongod = None
            #
            #     for i in range(n_steps):
            #         sample = _tuple[0]
            #         vcf_file = _tuple[1]
            #         csv_file = _tuple[2]
            #         db_name = _tuple[3]
            #         collection_name = _tuple[4]
            #         extra_data = extra_data
            #         step = i
            #         new_tuple_list.append((sample, vcf_file, csv_file, extra_data, db_name, collection_name, step, mongod,
            #                                definitions.myvariant_fields,
            #                                self.genome_build_version))
            #     return new_tuple_list
            #
            # def _make_simple_tuple_list(self, _tuple, n_steps):
            #     new_tuple_list = []
            #     for i in range(n_steps):
            #         step = i
            #         vcf_path = _tuple[0]
            #         db_name = _tuple[1]
            #         collec_name = _tuple[2]
            #         new_tuple_list.append((vcf_path,
            #                                db_name,
            #                                collec_name,
            #                                step,
            #                                definitions.myvariant_fields,
            #                                self.genome_build_version))
            #     return new_tuple_list

    @staticmethod
    def _make_jobs_params_tuples_list(file_path, chunk_size, shared_params_list):
        num_lines = _get_num_lines_in_file(file_path)
        num_steps = int(num_lines / chunk_size) + 1

        jobs_params_tuples_list = []
        for curr_chunk_index in range(num_steps):
            curr_job_params_tuple = tuple([curr_chunk_index] + shared_params_list)
            jobs_params_tuples_list.append(curr_job_params_tuple)

        return jobs_params_tuples_list


#
#
# def _parallel_get_dict_mv(maps):
#     vcf_path = maps[0]
#     db_name = maps[1]
#     collection_name = maps[2]
#     step = maps[3]
#     fields = maps[4]
#     genome_build_version = maps[5]
#
#     client = MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)
#     db = getattr(client, db_name)
#     collection = getattr(db, collection_name)
#
#     list_hgvs_ids = _get_hgvs_ids_from_vcf(vcf_path, step, definitions.chunk_size)
#     myvariants_variants = _get_myvariantinfo_annotations_dict(list_hgvs_ids, genome_build_version, 1)
#
#     # logging.info('Parsing Buffer...')
#     collection.insert_many(myvariants_variants, ordered=False)
#     client.close()
#     return None


def _parse_by_step(job_params_tuple):
    chunk_index = job_params_tuple[VariantParsing.CHUNK_INDEX_INDEX]
    file_path = job_params_tuple[VariantParsing.FILE_PATH_INDEX]
    chunk_size = job_params_tuple[VariantParsing.CHUNK_SIZE_INDEX]
    db_name = job_params_tuple[VariantParsing.DB_NAME_INDEX]
    collection_name = job_params_tuple[VariantParsing.COLLECTION_NAME_INDEX]
    genome_build_version = job_params_tuple[VariantParsing.GENOME_BUILD_VERSION_INDEX]
    verbose_level = job_params_tuple[VariantParsing.VERBOSE_LEVEL_INDEX]

    if len(job_params_tuple) > VariantParsing.SAMPLE_LIST_INDEX:
        merge_variants = True
        sample_names_list = job_params_tuple[VariantParsing.SAMPLE_LIST_INDEX]
        annovar_variants, hgvs_ids_list = AnnovarTxtParser.read_chunk_of_annotations_to_dicts_list(
            file_path, sample_names_list, chunk_index, chunk_size)
    else:
        merge_variants = False
        annovar_variants = None
        hgvs_ids_list = _get_hgvs_ids_from_vcf(file_path, chunk_index, chunk_size)

    myvariants_variants = _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version, verbose_level)

    variant_dicts_to_store = myvariants_variants
    if merge_variants:
        variant_dicts_to_store = []
        for i in range(0, len(hgvs_ids_list)):
            variant_dicts_to_store.append(
                _merge_annovar_and_myvariant_dicts(myvariants_variants[i], annovar_variants[i]))

    return _insert_handler(variant_dicts_to_store, db_name, collection_name)


def _get_hgvs_ids_from_vcf(vcf_file_path, step_number, chunk_size):
    reader = vcf.Reader(open(vcf_file_path, 'r'))
    hgvs_ids = []

    for record in itertools.islice(reader, step_number * chunk_size, (step_number + 1) * chunk_size):
        hgvs_id = myvariant.format_hgvs(record.CHROM, record.POS, record.REF, str(record.ALT[0]))
        normed_hgvs_id = _complete_chromosome(hgvs_id)
        hgvs_ids.append(normed_hgvs_id)

    return hgvs_ids


def _complete_chromosome(hgvs_id):
    """ Ensuring syntax consistency """

    result = hgvs_id
    if 'M' in hgvs_id:
        one = hgvs_id.split(':')[0]
        two = hgvs_id.split(':')[1]
        if 'MT' not in one:
            one = 'chrMT'
        result = "".join([one, ':', two])
    return result


def _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version, verbose):
    myvariant_fields = [
        'cadd.1000g',
        'cadd.esp',
        'cadd.phred',
        'cadd.gerp',
        'cadd.polyphen',
        'cadd.sift',
        'dbsnp.rsid',
        'cosmic.cosmic_id',
        'cosmic.tumor_site',
        'clinvar.rcv.accession',
        'clinvar.rcv.clinical_significance',
        'clinvar.rcv.conditions',
        'civic.description',
        'civic.evidence_items',
        'cgi',
        'gwassnps',
        'wellderly.alleles.freq'
    ]

    """ Retrieve variants from MyVariant.info"""
    if verbose >= 2:
        verbose = True
    else:
        verbose = False

    mv = myvariant.MyVariantInfo()
    # This will retrieve a list of dictionaries
    try:
        getvariants = mv.getvariants(hgvs_ids_list, verbose=1, as_dataframe=False, fields=myvariant_fields,
                                     assembly=genome_build_version)
        variant_data = getvariants
    except Exception as error:
        logging.info('Error: ' + str(error) + 'while fetching from MyVariant, retrying...')
        time.sleep(5)
        variant_data = _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version, verbose)

    variant_data = _remove_id_key(variant_data)
    return variant_data


def _remove_id_key(variant_data):
    """ Let mongo create an _id key to prevent insert attempts of documents with same key """

    for dic in variant_data:
        dic['hgvs_id'] = dic.pop("_id", None)
        dic['hgvs_id'] = dic.pop("query", None)
    return variant_data


def _merge_annovar_and_myvariant_dicts(myvariant_dict, annovar_dict):
    """
    Merge myvariant_dict with annovar_dict
    """
    if myvariant_dict['hgvs_id'] != annovar_dict['hgvs_id']:
        raise ValueError("myvariant hgvs_id {0} not equal to annovar hgvs_id {1}".format(myvariant_dict['hgvs_id'],
                                                                                         annovar_dict['hgvs_id']))
    annovar_dict.update(myvariant_dict)
    return annovar_dict


def _insert_handler(merged_list, db_name, collection_name, mongod_cmd=None):
    client = MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)
    db = getattr(client, db_name)
    collection = getattr(db, collection_name)

    logging.info('Parsing Buffer...')
    if len(merged_list) == 0:
        logging.info('Empty list of documents trying to be parsed, skipping and continuing operation...')
        return None
    else:
        try:
            collection.insert_many(merged_list, ordered=False)
        except Exception as error:
            if "Connection refused" in str(error):
                if mongod_cmd:
                    logging.info('MongoDB Server seems to be off. Attempting restart...')
                    args = shlex.split(mongod_cmd)
                    subprocess.Popen(args, stdout=subprocess.PIPE)
                else:
                    logging.error('Error connecting to mongodb.' + str(error))

                logging.info('Retrying to parse...')
                time.sleep(4)

                # TODO: Think about whether to keep this; is it necessary/feasible if moved client in here?
                _insert_handler(merged_list, collection, mongod_cmd)  # Recurse!
            else:
                logging.error('Unrecoverable error: ' + str(error))

    client.close()

# def _get_dict_myvariant(variant_list, verbose, fields, genome_build_version):
#     """ Retrieve variants from MyVariant.info"""
#     if verbose >= 2:
#         verbose = True
#     else:
#         verbose = False
#
#     mv = myvariant.MyVariantInfo()
#     # This will retrieve a list of dictionaries
#     try:
#         getvariants = mv.getvariants(variant_list, verbose=1, as_dataframe=False, fields=fields,
#                                      assembly=genome_build_version)
#         variant_data = getvariants
#     except Exception as error:
#         logging.info('Error: ' + str(error) + 'while fetching from MyVariant, retrying...')
#         time.sleep(5)
#         variant_data = _get_dict_myvariant(variant_list, verbose, fields, genome_build_version)
#
#     variant_data = _remove_id_key(variant_data)
#     return variant_data
