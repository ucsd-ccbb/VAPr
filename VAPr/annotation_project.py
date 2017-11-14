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
import VAPr.vcf_merge
from VAPr.annovar_runner import AnnovarWrapper
from VAPr.annovar_output_parsing import AnnovarTxtParser, AnnovarAnnotatedVariant
import VAPr.queries

# TODO: Understand, vet this logging set-up
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class AnnotationProject:
    HG19_VERSION = "hg19"
    HG38_VERSION = "hg38"
    DEFAULT_GENOME_VERSION = HG19_VERSION
    SUPPORTED_GENOME_BUILD_VERSIONS = [HG19_VERSION, HG38_VERSION]

    @staticmethod
    def _make_jobs_params_tuples_list(file_path, num_file_lines, chunk_size, db_name, collection_name,
                                      genome_build_version, sample_names_list=None, verbose_level=1):

        num_params = AnnotationJobParamsIndices.get_num_possible_indices()
        if sample_names_list is not None:
            shared_job_params = [None] * num_params
            shared_job_params[AnnotationJobParamsIndices.SAMPLE_LIST_INDEX] = sample_names_list
        else:
            shared_job_params = [None] * (num_params - 1)

        shared_job_params[AnnotationJobParamsIndices.CHUNK_SIZE_INDEX] = chunk_size
        shared_job_params[AnnotationJobParamsIndices.FILE_PATH_INDEX] = file_path
        shared_job_params[AnnotationJobParamsIndices.DB_NAME_INDEX] = db_name
        shared_job_params[AnnotationJobParamsIndices.COLLECTION_NAME_INDEX] = collection_name
        shared_job_params[AnnotationJobParamsIndices.GENOME_BUILD_VERSION_INDEX] = genome_build_version
        shared_job_params[AnnotationJobParamsIndices.VERBOSE_LEVEL_INDEX] = verbose_level

        jobs_params_tuples_list = []
        num_steps = int(num_file_lines / chunk_size) + 1
        for curr_chunk_index in range(num_steps):
            shared_job_params[AnnotationJobParamsIndices.CHUNK_INDEX_INDEX] = curr_chunk_index
            curr_job_params_tuple = tuple(shared_job_params)
            jobs_params_tuples_list.append(curr_job_params_tuple)

        return jobs_params_tuples_list

    @classmethod
    def _get_validated_genome_version(cls, input_genome_build_version):
        """ Make sure genome version is acceptable """

        if input_genome_build_version is None:
            result = cls.DEFAULT_GENOME_VERSION
        elif input_genome_build_version not in cls.SUPPORTED_GENOME_BUILD_VERSIONS:
            str_of_acceptable_versions = ", ".join(cls.SUPPORTED_GENOME_BUILD_VERSIONS)
            raise ValueError("Input genome build version '{0}' is not recognized. Supported builds are {1}".format(
                input_genome_build_version, str_of_acceptable_versions))
        else:
            result = input_genome_build_version

        return result

    def __init__(self, input_dir, output_dir, analysis_name, path_to_annovar_install, vcf_file_extension, mongo_db_name,
                 mongo_collection_name, design_file=None, build_ver=None, mongod_cmd=None):

        self._input_dir = input_dir
        self._output_dir = output_dir
        self._analysis_name = analysis_name
        self._design_file = design_file
        self._path_to_annovar_install = path_to_annovar_install
        self._vcf_file_extension = vcf_file_extension
        self._mongo_db_name = mongo_db_name
        self._mongo_collection_name = mongo_collection_name
        self._genome_build_version = self._get_validated_genome_version(build_ver)
        # self.mongod = mongod_cmd

        self._single_vcf_path, self._annovar_output_basename, self._sample_names_list = VAPr.vcf_merge.merge_vcfs(
            self._input_dir, self._output_dir, self._design_file, self._analysis_name, self._vcf_file_extension)

        self.annovar_wrapper = AnnovarWrapper(self._input_dir, self._output_dir, self._path_to_annovar_install,
                                              self._single_vcf_path, self._annovar_output_basename,
                                              design_file=self._design_file,
                                              genome_build_version=self._genome_build_version)

    def download_annovar_databases(self):
        """Run ANNOVAR to download its databases."""
        self.annovar_wrapper.download_dbs()

    def gather_basic_annotations(self, num_processes=8, chunk_size=2000, verbose_level=1):
        self._collect_annotations_and_store(self._single_vcf_path, chunk_size, num_processes, sample_names_list=None,
                                            verbose_level=verbose_level)

    def gather_detailed_annotations(self, num_processes=4, chunk_size=2000, verbose_level=1, multisample=False):
        annovar_output_fp = self._run_annovar_annotation(multisample)
        self._collect_annotations_and_store(annovar_output_fp, chunk_size, num_processes,
                                            sample_names_list=self._sample_names_list, verbose_level=verbose_level)

    def write_output_files_by_sample(self):
        # TODO: finish refactoring this functionality
        raise NotImplementedError("function has not been refactored yet")
        # generate_output_files_by_sample(self._mongo_db_name, self._mongo_collection_name,  self._output_dir)

    def _run_annovar_annotation(self, multisample=False):
        """Run ANNOVAR to annotate variants in a vcf file."""
        return self.annovar_wrapper.run_annovar(vcf_is_multisample=multisample)

    # TODO: someday: extra_data from design file needs to come back in here ...
    def _collect_annotations_and_store(self, file_path, chunk_size, num_processes, sample_names_list=None,
                                       verbose_level=1):
        num_file_lines = _get_num_lines_in_file(file_path)
        jobs_params_tuples_list = self._make_jobs_params_tuples_list(
            file_path, num_file_lines, chunk_size, self._mongo_db_name, self._mongo_collection_name,
            self._genome_build_version, sample_names_list, verbose_level)

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

    # TODO: someday: refactor so one doesn't have to remember to add new indices to the below function
    @classmethod
    def get_num_possible_indices(cls):
        max_index = max(cls.CHUNK_INDEX_INDEX, cls.FILE_PATH_INDEX, cls.CHUNK_SIZE_INDEX, cls.DB_NAME_INDEX,
                        cls.COLLECTION_NAME_INDEX, cls.GENOME_BUILD_VERSION_INDEX, cls.VERBOSE_LEVEL_INDEX,
                        cls.SAMPLE_LIST_INDEX)
        return max_index+1


def _get_num_lines_in_file(file_path):
    with open(file_path) as file_obj:
        result = sum(1 for _ in file_obj)
    return result


def _collect_chunk_annotations_and_store(job_params_tuple):
    db_name = job_params_tuple[AnnotationJobParamsIndices.DB_NAME_INDEX]
    collection_name = job_params_tuple[AnnotationJobParamsIndices.COLLECTION_NAME_INDEX]
    variant_dicts_to_store = _collect_chunk_annotations(job_params_tuple)
    VAPr.queries.store_annotations_to_db(variant_dicts_to_store, db_name, collection_name)


def _collect_chunk_annotations(job_params_tuple):
    chunk_index = job_params_tuple[AnnotationJobParamsIndices.CHUNK_INDEX_INDEX]
    chunk_size = job_params_tuple[AnnotationJobParamsIndices.CHUNK_SIZE_INDEX]
    file_path = job_params_tuple[AnnotationJobParamsIndices.FILE_PATH_INDEX]
    genome_build_version = job_params_tuple[AnnotationJobParamsIndices.GENOME_BUILD_VERSION_INDEX]
    verbose_level = job_params_tuple[AnnotationJobParamsIndices.VERBOSE_LEVEL_INDEX]

    with open(file_path, 'r') as input_file_obj:
        if len(job_params_tuple) > AnnotationJobParamsIndices.SAMPLE_LIST_INDEX:
            merge_variants = True
            sample_names_list = job_params_tuple[AnnotationJobParamsIndices.SAMPLE_LIST_INDEX]
            hgvs_ids_list, annovar_variants = AnnovarTxtParser.read_chunk_of_annotations_to_dicts_list(
                input_file_obj, sample_names_list, chunk_index, chunk_size)
        else:
            merge_variants = False
            annovar_variants = None
            hgvs_ids_list = _get_hgvs_ids_from_vcf(input_file_obj, chunk_index, chunk_size)

    myvariants_variants = _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version,
                                                              verbose_level)

    result = myvariants_variants
    if merge_variants:
        result = []
        for i in range(0, len(hgvs_ids_list)):
            result.append(_merge_annovar_and_myvariant_dicts(myvariants_variants[i], annovar_variants[i]))

    return result


def _get_hgvs_ids_from_vcf(vcf_file_obj, chunk_index, chunk_size):
    reader = vcf.Reader(vcf_file_obj)
    hgvs_ids = []

    for record in itertools.islice(reader, chunk_index * chunk_size, (chunk_index + 1) * chunk_size):
        hgvs_id = myvariant.format_hgvs(record.CHROM, record.POS, record.REF, str(record.ALT[0]))

        # ensure syntax consistency for chromosome M variants
        if 'M' in hgvs_id:
            one = hgvs_id.split(':')[0]
            two = hgvs_id.split(':')[1]
            if 'MT' not in one:
                one = 'chrMT'
                hgvs_id = "".join([one, ':', two])

        hgvs_ids.append(hgvs_id)

    return hgvs_ids


def _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version, verbose_level, num_failed_attempts=0):
    """ Retrieve variants from MyVariant.info"""

    max_failed_attempts = 5
    # TODO: Ask Adam who picked these fields; e.g., why only getting freq from wellderly.alleles but not allele??
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

    be_verbose = verbose_level >= 2
    mv = myvariant.MyVariantInfo()
    try:
        myvariantinfo_dicts_list = mv.getvariants(hgvs_ids_list, verbose=int(be_verbose), as_dataframe=False,
                                                  fields=myvariant_fields, assembly=genome_build_version)
    except Exception as error:
        # TODO: Discuss w/Adam what kinds of errors really should be retried; some def shouldn't
        logging.info('Error: ' + str(error) + 'while fetching from MyVariant')
        num_failed_attempts += 1
        if num_failed_attempts < max_failed_attempts:
            time.sleep(5)
            logging.info("Retrying MyVariant.info fetch")
            myvariantinfo_dicts_list = _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version,
                                                                           verbose_level, num_failed_attempts)
        else:
            # give up and raise error
            raise error

    myvariantinfo_dicts_list = _remove_unwanted_keys(myvariantinfo_dicts_list)
    return myvariantinfo_dicts_list


def _remove_unwanted_keys(myvariantinfo_dicts_list):
    for curr_myvariantinfo_dict in myvariantinfo_dicts_list:
        # Put the contents of the query key into a new hgvs_id key; NB, use AnnovarAnnotatedVariant.HGVS_ID_KEY
        # as the value of this key because later, when we combine this dict with the annovar-generated dict (if any),
        # we want to make sure the keys are duplicates and thus collapse into one.
        curr_myvariantinfo_dict[AnnovarAnnotatedVariant.HGVS_ID_KEY] = curr_myvariantinfo_dict.pop("query", None)
        # Also, just drop the _id key--we want mongo db to create its own _id key.
        curr_myvariantinfo_dict.pop("_id", None)
    return myvariantinfo_dicts_list


def _merge_annovar_and_myvariant_dicts(myvariantinfo_annotations_dict, annovar_annotations_dict):
    hgvs_id_key = AnnovarAnnotatedVariant.HGVS_ID_KEY
    if myvariantinfo_annotations_dict[hgvs_id_key] != annovar_annotations_dict[hgvs_id_key]:
        raise ValueError(
            "myvariant HGVS id '{0}' not equal to annovar HGVS id '{1}'".format(
                myvariantinfo_annotations_dict[hgvs_id_key], annovar_annotations_dict[hgvs_id_key]))

    annovar_annotations_dict.update(myvariantinfo_annotations_dict)
    return annovar_annotations_dict


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
