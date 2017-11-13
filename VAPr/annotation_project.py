from __future__ import division, print_function

# built-in libraries
import itertools
import logging
import multiprocessing
import sys
import time
import tqdm

# third-party libraries
import myvariant
import vcf

# project libraries
from VAPr.vcf_merge import MergeVcfs
from VAPr.annovar_runner import AnnovarWrapper
import VAPr.vcf_mappings_maker
from VAPr.annovar_output_parsing import AnnovarTxtParser
import VAPr.queries

# TODO: Understand, vet this logging set-up
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class AnnotationProject:
    @staticmethod
    def _make_jobs_params_tuples_list(file_path, chunk_size, db_name, collection_name, genome_build_version,
                                      sample_names_list=None, verbose_level=1):
        shared_job_params = []
        shared_job_params[AnnotationJobParamsIndices.CHUNK_SIZE_INDEX] = chunk_size
        shared_job_params[AnnotationJobParamsIndices.FILE_PATH_INDEX] = file_path
        shared_job_params[AnnotationJobParamsIndices.DB_NAME_INDEX] = db_name
        shared_job_params[AnnotationJobParamsIndices.COLLECTION_NAME_INDEX] = collection_name
        shared_job_params[AnnotationJobParamsIndices.GENOME_BUILD_VERSION_INDEX] = genome_build_version
        shared_job_params[AnnotationJobParamsIndices.VERBOSE_LEVEL_INDEX] = verbose_level

        if sample_names_list is not None:
            shared_job_params[AnnotationJobParamsIndices.SAMPLE_LIST_INDEX] = sample_names_list

        num_lines = _get_num_lines_in_file(file_path)
        num_steps = int(num_lines / chunk_size) + 1

        jobs_params_tuples_list = []
        for curr_chunk_index in range(num_steps):
            curr_job_params_tuple = tuple([curr_chunk_index] + shared_job_params)
            jobs_params_tuples_list.append(curr_job_params_tuple)

        return jobs_params_tuples_list

    def __init__(self, input_dir, output_dir, analysis_name, annovar_path, vcf_file_extension, mongo_db_and_collection_names_dict,


    def __init__(self, input_dir, output_dir, analysis_name, annovar_path, mongo_db_name, mongo_collection_name,
                 design_file=None, build_ver=None, mongod_cmd=None, split_vcf=False):

        self._input_dir = input_dir
        self._output_dir = output_dir
        self._analysis_name = analysis_name
        self._design_file = design_file
        self._path_to_annovar_install = annovar_path
        self._mongo_db_name = mongo_db_name
        self._mongo_collection_name = mongo_collection_name
        self._genome_build_version = AnnovarWrapper.get_validated_genome_version(build_ver)
        # self.times_called = 0
        # self.split = split_vcf
        # self.mongod = mongod_cmd
        self._vcf_file_extension = vcf_file_extension
        # return vcf mapping dict
        self.vcf_mapping_dict = MergeVcfs(self._input_dir,
                                          self._output_dir,
                                          self._analysis_name,
                                          self._design_file,
                                          self._vcf_file_extension).merge_vcfs()

        self.annovar_wrapper = AnnovarWrapper(self._input_dir, self._output_dir, self._path_to_annovar_install,
                                              self._make_mongo_db_and_collection_names_dict(), self.vcf_mapping_dict,
                                              design_file=self._design_file,
                                              genome_build_version=self._genome_build_version)

    def download_annovar_databases(self):
        """Run ANNOVAR to download its databases."""
        self.annovar_wrapper.download_dbs()

    def gather_basic_annotations(self, num_processes=8, chunk_size=2000, verbose_level=1):
        vcf_fp = self.vcf_mapping_dict[VAPr.vcf_mappings_maker.SingleVcfFileMappingMaker.VCF_FILE_PATH_KEY]

        self._collect_annotations_and_store(vcf_fp, chunk_size, num_processes, sample_names_list=None,
                                            verbose_level=verbose_level)

    def gather_detailed_annotations(self, num_processes=4, chunk_size=2000, verbose_level=1, multisample=False):
        sample_names_list = self.vcf_mapping_dict[
            VAPr.vcf_mappings_maker.SingleVcfFileMappingMaker.SAMPLE_NAMES_LIST_KEY]

        annovar_output_fp = self._run_annovar_annotation(multisample)
        self._collect_annotations_and_store(annovar_output_fp, chunk_size, num_processes,
                                            sample_names_list=sample_names_list, verbose_level=verbose_level)

    def write_output_files_by_sample(self):
        # TODO: finish refactoring this functionality
        raise NotImplementedError("function has not been refactored yet")
        generate_output_files_by_sample(self._mongo_db_name, self._mongo_collection_name,  self._output_dir)

    # TODO: someday: update AnnovarWrapper init to take db_name and collection_name as individual arguments
    # However, for now, to avoid breaking existing interface, continue to send in as a dict.
    def _make_mongo_db_and_collection_names_dict(self):
        return {'db_name': self._mongo_db_name, 'collection_name': self._mongo_collection_name}

    def _run_annovar_annotation(self, multisample=False):
        """Run ANNOVAR to annotate variants in a vcf file."""
        return self.annovar_wrapper.run_annovar(vcf_is_multisample=multisample)

    # TODO: someday: extra_data from design file needs to come back in here ...
    def _collect_annotations_and_store(self, file_path, chunk_size, num_processes, sample_names_list=None,
                                      verbose_level=1):
        jobs_params_tuples_list = self._make_jobs_params_tuples_list(
            file_path, chunk_size, self._mongo_db_name, self._mongo_collection_name, self._genome_build_version,
            sample_names_list, verbose_level)

        pool = multiprocessing.Pool(num_processes)
        for _ in tqdm.tqdm(pool.imap_unordered(_collect_chunk_annotations_and_store, jobs_params_tuples_list),
                           total=len(jobs_params_tuples_list)):
            pass
        pool.close()
        pool.join()

        logger.info('Done collecting and saving annotations')


class AnnotationJobParamsIndices:
    CHUNK_INDEX_INDEX = 0
    FILE_PATH_INDEX = 1
    CHUNK_SIZE_INDEX = 2
    DB_NAME_INDEX = 3
    COLLECTION_NAME_INDEX = 4
    GENOME_BUILD_VERSION_INDEX = 5
    VERBOSE_LEVEL_INDEX = 6
    SAMPLE_LIST_INDEX = 7


def _get_num_lines_in_file(file_path):
    return sum(1 for _ in open(file_path))


def _collect_chunk_annotations_and_store(job_params_tuple):
    chunk_index = job_params_tuple[AnnotationJobParamsIndices.CHUNK_INDEX_INDEX]
    chunk_size = job_params_tuple[AnnotationJobParamsIndices.CHUNK_SIZE_INDEX]
    file_path = job_params_tuple[AnnotationJobParamsIndices.FILE_PATH_INDEX]
    db_name = job_params_tuple[AnnotationJobParamsIndices.DB_NAME_INDEX]
    collection_name = job_params_tuple[AnnotationJobParamsIndices.COLLECTION_NAME_INDEX]
    genome_build_version = job_params_tuple[AnnotationJobParamsIndices.GENOME_BUILD_VERSION_INDEX]
    verbose_level = job_params_tuple[AnnotationJobParamsIndices.VERBOSE_LEVEL_INDEX]

    if len(job_params_tuple) > AnnotationJobParamsIndices.SAMPLE_LIST_INDEX:
        merge_variants = True
        sample_names_list = job_params_tuple[AnnotationJobParamsIndices.SAMPLE_LIST_INDEX]
        annovar_variants, hgvs_ids_list = AnnovarTxtParser.read_chunk_of_annotations_to_dicts_list(
            file_path, sample_names_list, chunk_index, chunk_size)
    else:
        merge_variants = False
        annovar_variants = None
        hgvs_ids_list = _get_hgvs_ids_from_vcf(file_path, chunk_index, chunk_size)

    myvariants_variants = _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version,
                                                              verbose_level)

    variant_dicts_to_store = myvariants_variants
    if merge_variants:
        variant_dicts_to_store = []
        for i in range(0, len(hgvs_ids_list)):
            variant_dicts_to_store.append(
                _merge_annovar_and_myvariant_dicts(myvariants_variants[i], annovar_variants[i]))

    return VAPr.queries.store_annotations_to_db(variant_dicts_to_store, db_name, collection_name)


def _get_hgvs_ids_from_vcf(vcf_file_path, chunk_index, chunk_size):
    reader = vcf.Reader(open(vcf_file_path, 'r'))
    hgvs_ids = []

    for record in itertools.islice(reader, chunk_index * chunk_size, (chunk_index + 1) * chunk_size):
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
    """ Retrieve variants from MyVariant.info"""

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
        raise ValueError(
            "myvariant hgvs_id {0} not equal to annovar hgvs_id {1}".format(myvariant_dict['hgvs_id'],
                                                                            annovar_dict['hgvs_id']))
    annovar_dict.update(myvariant_dict)
    return annovar_dict


    # def _find_annovar_output_file_path(self):
    #     # Get the annovar output file
    #     annovar_output_dir = self.vcf_mapping_dict[VAPr.vcf_mappings_maker.SingleVcfFileMappingMaker.ANNOVAR_OUTPUT_DIR]
    #     annovar_output_basename = self.vcf_mapping_dict[
    #         VAPr.vcf_mappings_maker.SingleVcfFileMappingMaker.ANNOVAR_OUTPUT_BASENAME_KEY]
    #     annovar_output_file_names = [curr_file_name for curr_file_name in os.listdir(annovar_output_dir) if
    #                                  curr_file_name.startswith(annovar_output_basename) and
    #                                  curr_file_name.endswith('txt')]
    #
    #     if len(annovar_output_file_names) > 1:
    #         raise ValueError('Too many matching annovar output files found')
    #     elif len(annovar_output_file_names) == 0:
    #         raise ValueError('No matching annovar output files found')
    #
    #     result = os.path.join(annovar_output_dir, annovar_output_file_names[0])
    #     return result
