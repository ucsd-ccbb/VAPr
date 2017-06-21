from __future__ import division, print_function

import myvariant
import os
import sys
from pymongo import MongoClient
import VAPr.definitions as definitions
from VAPr.models import TxtParser, HgvsParser
from VAPr.writes import Writer
from VAPr.queries import Filters
import tqdm
from multiprocessing import Pool
import logging
import subprocess
import shlex
import time
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'


class VariantParsing:

    def __init__(self,
                 input_dir,
                 output_csv_path,
                 annovar_path,
                 project_data,
                 mapping,
                 design_file=None,
                 build_ver=None,
                 mongod_cmd=None):

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
        self.n_vars = 0
        self.mongod = mongod_cmd

    def get_sample_csv_vcf_tuple(self):
        """ Locate files associated with a specific sample """

        list_tupls = []

        for _map in self.mapping:

            matching_csv = [i for i in os.listdir(_map['csv_file_full_path']) if i.startswith(_map['csv_file_basename'])
                            and i.endswith('txt')]

            matching_vcf = [i for i in os.listdir(_map['csv_file_full_path']) if i.startswith(_map['csv_file_basename'])
                            and i.endswith('vcf')]

            if len(matching_csv) > 1 or len(matching_vcf) > 1:
                raise ValueError('Too many matching csvs')
            elif len(matching_csv) == 0 or len(matching_vcf) == 0:
                raise ValueError('Csv not found')
            else:
                csv_path = os.path.join(_map['csv_file_full_path'], matching_csv[0])
                vcf_path = os.path.join(_map['csv_file_full_path'], matching_vcf[0])
                list_tupls.append((_map['sample_names'],
                                   vcf_path,
                                   csv_path,
                                   self.db,
                                   self.collection,
                                   _map['extra_data']))

        return list_tupls

    def parallel_annotation(self, n_processes, verbose=1):

        """
        Set up variant parsing scheme. Since a functional programming style is required for parallel
        processing, the input to the Pool.map function from the multiprocessing library must be
        immutable. We use tuples of the type:

            (samples, vcf_path, csv_path, db_name, collection_name, step)

        The parallel annotation will split the vcf, csv file in N steps where N equals the number of lines
        (variants) to be annotated divided a pre-defined 'chunk size'. For instance, a vcf, csv pair containing
        250 000 variants and a chunksize of 5 000 will yield 50 steps. Given 8 parallel procrsses, the annotation
        will happen over those 50 steps with 8 cores working cuncurrently until the 50 jobs are fully annotated.

        :param n_processes: number of cores to be used
        :param verbose: verbosity level [0,1,2,3]
        :return: None
        """

        self.verbose = verbose
        list_tupls = self.get_sample_csv_vcf_tuple()
        map_job = []
        for tpl in list_tupls:
            hgvs = HgvsParser(tpl[1])
            csv_parsing = TxtParser(tpl[2], samples=hgvs.samples, extra_data=tpl[5])
            num_lines = csv_parsing.num_lines
            n_steps = int(num_lines/self.chunksize) + 1
            map_job.extend(self.parallel_annotator_mapper(tpl, n_steps, extra_data=tpl[5], mongod_cmd=self.mongod))

        pool = Pool(n_processes)
        for _ in tqdm.tqdm(pool.imap_unordered(parse_by_step, map_job), total=len(map_job)):
            pass
        pool.close()
        pool.join()
        # logger.info('Completed annotation and parsing for variants in sample %s' % tpl[0])

    @staticmethod
    def parallel_annotator_mapper(_tuple, n_steps, extra_data=None, mongod_cmd=None):
        """ Assign step number to each tuple to be consumed by parsing function """
        new_tuple_list = []
        if mongod_cmd:
            mongod = mongod_cmd
        else:
            mongod = None

        for i in range(n_steps):
            sample = _tuple[0]
            vcf_file = _tuple[1]
            csv_file = _tuple[2]
            db_name = _tuple[3]
            collection_name = _tuple[4]
            extra_data = extra_data
            step = i
            new_tuple_list.append((sample,
                                   vcf_file,
                                   csv_file,
                                   extra_data,
                                   db_name,
                                   collection_name,
                                   step,
                                   mongod))
        return new_tuple_list

    def quick_annotate_and_save(self, n_processes=8):
        """ Annotation that doesn't require annovar """

        list_tuples = []
        for _map in self.mapping:
            simple_map = {k: v for k, v in _map.items() if k not in ['csv_file_basename', 'sample_names']}
            list_tuples.append((simple_map['raw_vcf_file_full_path'], self.db, self.collection))

        for _tuple in list_tuples:

            hgvs = HgvsParser(_tuple[0])
            num_lines = hgvs.get_num_lines()
            n_steps = int(num_lines / self.chunksize) + 1
            map_job = self.quick_annotate_mapper(_tuple, n_steps)
            pool = Pool(n_processes)
            for _ in tqdm.tqdm(pool.imap_unordered(parallel_get_dict_mv, map_job), total=len(map_job)):
                pass

    @staticmethod
    def quick_annotate_mapper(_tuple, n_steps):
        new_tuple_list = []
        for i in range(n_steps):
            step = i
            vcf_path = _tuple[0]
            db_name = _tuple[1]
            collec_name = _tuple[2]
            new_tuple_list.append((vcf_path,
                                   db_name,
                                   collec_name,
                                   step))
        return new_tuple_list

    def generate_output_files_by_sample(self):

        client = MongoClient()
        db = getattr(client, self.project_data['db_name'])
        collection = getattr(db, self.project_data['collection_name'])

        # Get all distinct sample_ids
        samples = collection.distinct('sample_id')

        fwriter = Writer(self.project_data['db_name'], self.project_data['collection_name'])
        filt = Filters(self.project_data['db_name'], self.project_data['collection_name'])

        # Generate files for each sample
        for sample in samples:
            q = filt.variants_from_sample(sample)
            list_docs = list(q)
            out_path_by_sample = os.path.join(self.output_csv_path, 'csv_by_sample')
            if not os.path.exists(out_path_by_sample):
                os.makedirs(out_path_by_sample)
            fname = os.path.join(out_path_by_sample, 'annotated_csv_' + sample + '_all_vars.csv')
            fwriter.generate_annotated_csv(list_docs, fname)


def parse_by_step(maps):
    """ The function that implements the parsing """

    sample = maps[0]
    vcf_file = maps[1]
    csv_file = maps[2]
    extra_data = maps[3]
    db_name = maps[4]
    collection_name = maps[5]
    step = maps[6]
    mongod_cmd = maps[7]

    client = MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)
    db = getattr(client, db_name)
    collection = getattr(db, collection_name)
    hgvs = HgvsParser(vcf_file)
    csv_parsing = TxtParser(csv_file, samples=sample, extra_data=extra_data)
    list_hgvs_ids = hgvs.get_variants_from_vcf(step)
    myvariants_variants = get_dict_myvariant(list_hgvs_ids, 1, sample)

    csv_variants = csv_parsing.open_and_parse_chunks(step, build_ver='hg19')

    merged_list = []
    for i, _ in enumerate(myvariants_variants):
        for dict_from_sample in csv_variants[i]:
            merged_list.append(merge_dict_lists(myvariants_variants[i], dict_from_sample))

    return insert_handler(merged_list, collection, mongod_cmd)


def insert_handler(merged_list, collection, mongod_cmd):

    logging.info('Parsing Buffer...')
    if len(merged_list) == 0:
        logging.info('Empty list of documents trying to be parsed, skipping and continuing operation...')
        return None
    else:
        try:
            collection.insert_many(merged_list, ordered=False)
        except Exception, error:
            if "Connection refused" in str(error):
                if mongod_cmd:
                    logging.info('MongoDB Server seems to be off. Attempting restart...')
                    print(mongod_cmd)
                    args = shlex.split(mongod_cmd)
                    subprocess.Popen(args, stdout=subprocess.PIPE)
                else:
                    logging.error('Error connecting to mongodb.' + str(error))
                logging.info('Retrying to parse...')
                time.sleep(4)
                insert_handler(merged_list, collection, mongod_cmd)  # Recurse!
            else:
                logging.error('Unrecoverable error: ' + str(error))

    return None


def merge_dict_lists(*dict_args):
    """
    Given any number of dicts, shallow copy and merge into a new dict,
    precedence goes to key value pairs in latter dicts.
    """
    result = {}
    for dictionary in dict_args:
        result.update(dictionary)
    return result


def parallel_get_dict_mv(maps):

    vcf_path = maps[0]
    db_name = maps[1]
    collection_name = maps[2]
    step = maps[3]

    client = MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)
    db = getattr(client, db_name)
    collection = getattr(db, collection_name)
    hgvs = HgvsParser(vcf_path)

    list_hgvs_ids = hgvs.get_variants_from_vcf(step)
    myvariants_variants = get_dict_myvariant(list_hgvs_ids, 1, hgvs.samples)

    # logging.info('Parsing Buffer...')
    collection.insert_many(myvariants_variants, ordered=False)
    client.close()
    return None


def get_dict_myvariant(variant_list, verbose, sample_id):
    """ Retrieve variants from MyVariant.info"""

    if verbose >= 2:
        verbose = True
    else:
        verbose = False

    mv = myvariant.MyVariantInfo()
    # This will retrieve a list of dictionaries
    try:
        variant_data = mv.getvariants(variant_list, verbose=1, as_dataframe=False)
    except Exception, error:
        logging.info('Error: ' + str(error) + 'while fetching from MyVariant, retrying...')
        time.sleep(5)
        variant_data = get_dict_myvariant(variant_list, verbose, sample_id)

    variant_data = remove_id_key(variant_data, sample_id)
    return variant_data


def remove_id_key(variant_data, sample_id):
    """ Let mongo create an _id key to prevent insert attempts of documents with same key """

    for dic in variant_data:
        dic['hgvs_id'] = dic.pop("_id", None)
        dic['hgvs_id'] = dic.pop("query", None)
        dic['sample_id'] = sample_id

    return variant_data
