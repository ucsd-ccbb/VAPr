from __future__ import division, print_function

# built-in libraries
import itertools
import logging
import time

# third-party libraries
import myvariant
import pymongo
import vcf

# project libraries
from VAPr.annovar_output_parsing import AnnovarTxtParser, AnnovarAnnotatedVariant


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


def collect_chunk_annotations_and_store(job_params_tuple):
    db_name = job_params_tuple[AnnotationJobParamsIndices.DB_NAME_INDEX]
    collection_name = job_params_tuple[AnnotationJobParamsIndices.COLLECTION_NAME_INDEX]
    variant_dicts_to_store = _collect_chunk_annotations(job_params_tuple)
    _store_annotations_to_db(variant_dicts_to_store, db_name, collection_name)


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


# TODO: someday: refactor myvariant fields into external file so easy to modify which are pulled
def _get_myvariantinfo_annotations_dict(hgvs_ids_list, genome_build_version, verbose_level, num_failed_attempts=0):
    """ Retrieve variants from MyVariant.info"""

    max_failed_attempts = 5
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
        'wellderly.alleles'
    ]

    be_verbose = verbose_level >= 2
    mv = myvariant.MyVariantInfo()
    try:
        myvariantinfo_dicts_list = mv.getvariants(hgvs_ids_list, verbose=int(be_verbose), as_dataframe=False,
                                                  fields=myvariant_fields, assembly=genome_build_version)
    except ValueError as unrecoverable_error:
        # If myvariant.info returned a value error, recalling with the same values won't help so error out now
        raise unrecoverable_error
    except Exception as error:
        # If we got something other than a value error, problem may be with internet connection or myvariant.info
        # availability, so try again a couple of times just in case we can recover
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


def _store_annotations_to_db(annotation_dicts_list, db_name, collection_name, client=None, num_failed_attempts=0):
    max_failed_attempts = 5

    if client is None:
        client = pymongo.MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)

    db = getattr(client, db_name)
    collection = getattr(db, collection_name)

    logging.info('Parsing Buffer...')
    if len(annotation_dicts_list) == 0:
        logging.info('List of annotations to store is empty; continuing.')
        return

    try:
        collection.insert_many(annotation_dicts_list, ordered=False)
    except Exception as error:
        if "Connection refused" in str(error) and num_failed_attempts < max_failed_attempts:
            num_failed_attempts += 1
            logging.error('Error connecting to mongodb: ' + str(error))
            logging.info('Retrying connection to mongodb')
            time.sleep(4)

            _store_annotations_to_db(annotation_dicts_list, db_name, collection_name, client=client,
                                     num_failed_attempts=num_failed_attempts)  # Recurse!
        else:
            error_msg = "Encountered error '{0}' when attempting to store the following annotations " \
                        "to mongo db collection '{1}': '{2}'".format(str(error), collection.full_name,
                                                                     annotation_dicts_list)
            raise RuntimeError(error_msg)

    try:
        client.close()
    except:
        pass  # if the client is already closed, just move along
