from __future__ import division, print_function

# built-in libraries
import logging
import multiprocessing
import os
import pymongo
import tqdm
import warnings

# third-party libraries
import pandas
import vcf

# project libraries
import VAPr.vcf_merging
import VAPr.annovar_running
import VAPr.filtering
import VAPr.chunk_processing


class VaprDataset(object):
    def __init__(self, mongo_db_name, mongo_collection_name, merged_vcf_path=None):
        self._mongo_db_name = mongo_db_name
        self._mongo_collection_name = mongo_collection_name
        self._merged_vcf_path = merged_vcf_path

        self._mongo_client = pymongo.MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)
        self._mongo_db = getattr(self._mongo_client, self._mongo_db_name)
        self._mongo_db_collection = getattr(self._mongo_db, self._mongo_collection_name)

    @property
    def full_name(self):
        return self._mongo_db_collection.full_name

    def get_rare_deleterious_variants(self, sample_names_list=None):
        return self._get_filtered_variants_by_sample(VAPr.filtering.make_rare_deleterious_variants_filter,
                                                     sample_names_list)

    def get_known_disease_variants(self, sample_names_list=None):
        return self._get_filtered_variants_by_sample(VAPr.filtering.make_known_disease_variants_filter,
                                                     sample_names_list)

    def get_deleterious_compound_heterozygous_variants(self, sample_names_list=None):
        return self._get_filtered_variants_by_sample(
            VAPr.filtering.make_deleterious_compound_heterozygous_variants_filter, sample_names_list)

    def get_de_novo_variants(self, proband, ancestor1, ancestor2):
        filter_dict = VAPr.filtering.make_de_novo_variants_filter(proband, ancestor1, ancestor2)
        return self.get_custom_filtered_variants(filter_dict)

    def get_custom_filtered_variants(self, filter_dictionary):
        return list(self._mongo_db_collection.find(filter_dictionary))

    def get_distinct_sample_ids(self):
        result = self._mongo_db_collection.distinct(VAPr.filtering.SAMPLE_ID_SELECTOR)
        return result

    def get_all_variants(self):
        return self.get_custom_filtered_variants({})

    def get_variants_for_sample(self, sample_name):
        filter_dict = VAPr.filtering.get_sample_id_filter(sample_name)
        return self.get_custom_filtered_variants(filter_dict)

    def get_variants_for_samples(self, sample_names_list):
        filter_dict = VAPr.filtering.get_any_of_sample_ids_filter(sample_names_list)
        return self.get_custom_filtered_variants(filter_dict)

    def get_variants_as_dataframe(self, filtered_variants=None):
        if filtered_variants is None:
            filtered_variants = self.get_all_variants()
        return pandas.DataFrame(filtered_variants)

    def write_unfiltered_annotated_csv(self, output_fp):
        all_variants = self.get_all_variants()
        self._write_annotated_csv("write_unfiltered_annotated_csv", all_variants, output_fp)

    def write_filtered_annotated_csv(self, filtered_variants, output_fp):
        self._write_annotated_csv("write_filtered_annotated_csv", filtered_variants, output_fp)

    def write_unfiltered_annotated_vcf(self, vcf_output_path, info_out=True):
        filtered_variants = self.get_all_variants()
        self._write_annotated_vcf(filtered_variants, vcf_output_path, info_out=info_out)

    def write_filtered_annotated_vcf(self, filtered_variants, vcf_output_path, info_out=True):
        self._write_annotated_vcf(filtered_variants, vcf_output_path, info_out=info_out)

    def write_unfiltered_annotated_csvs_per_sample(self, output_dir):
        sample_ids_list = self.get_distinct_sample_ids()

        for curr_sample_id in sample_ids_list:
            variant_dicts_list = self.get_variants_for_sample(curr_sample_id)
            curr_output_fp = os.path.join(output_dir, curr_sample_id + 'unfiltered_annotated_variants.csv')
            self.write_filtered_annotated_csv(variant_dicts_list, curr_output_fp)

        self._warn_if_no_output("write_unfiltered_annotated_csvs_per_sample", sample_ids_list)

    def _construct_sample_ids_list(self, sample_names):
        result = sample_names
        if not sample_names:
            result = self.get_distinct_sample_ids()
        elif not isinstance(sample_names, list):
            result = [sample_names]

        if len(result) == 0:
            warnings.warn("No sample ids found.")
        return result

    def _write_annotated_csv(self, func_name, filtered_variants, output_fp):
        no_output = self._warn_if_no_output(func_name, filtered_variants)
        if not no_output:
            dataframe = self.get_variants_as_dataframe(filtered_variants)
            dataframe.to_csv(output_fp)

    def _get_filtered_variants_by_sample(self, filter_builder_func, sample_names=None):
        sample_ids_list = self._construct_sample_ids_list(sample_names)
        filter_dict = filter_builder_func(sample_ids_list)
        return self.get_custom_filtered_variants(filter_dict)

    # TODO: I'd like to do a bit more work on this one; not sure it is in the right place
    def _write_annotated_vcf(self, filtered_variants_dicts_list, vcf_output_path, info_out=True):
        """
        :param vcf_input_path: template vcf file (initial vcf from which a new one will be created)
        :param vcf_output_path: name and filepath to where new vcf file will be written
        :param filtered_variants_dicts_list: list of dictionaries (one per variant) containing annotations
        :param info_out: if set to true (Default), will write all annotation data to INFO column, else, it won't.
        """

        # TODO: check if merged vcf path is None; if so, throw error
        # TODO: Must bgzip and index merged vcf path  before writing; use method from vcf_merging

        chr_vars = []
        location_vars_ant = []
        location_vars_pos = []

        for i in range(0, len(filtered_variants_dicts_list)):
            if filtered_variants_dicts_list[i]['chr'] == 'chrMT':
                chr_vars.append('chrM')
            else:
                chr_vars.append(filtered_variants_dicts_list[i]['chr'])
            location_vars_ant.append(filtered_variants_dicts_list[i]['start'] + 1)
            location_vars_pos.append(filtered_variants_dicts_list[i]['start'] - 1)

        vcf_reader = vcf.Reader(filename=self._merged_vcf_path)
        vcf_writer = vcf.Writer(open(vcf_output_path, 'w'), vcf_reader)

        for i in range(0, len(chr_vars)):
            for record in vcf_reader.fetch(chr_vars[i], location_vars_pos[i], location_vars_ant[i]):
                if info_out is True:
                    record.INFO.update(filtered_variants_dicts_list[i])
                    vcf_writer.write_record(record)
                else:
                    vcf_writer.write_record(record)

        self._warn_if_no_output("write_unfiltered_annotated_csvs_per_sample", filtered_variants_dicts_list)

    def _warn_if_no_output(self, output_func_name, items_list):
        no_output = False
        if len(items_list) == 0:
            no_output = True
            warnings.warn("{0} wrote no file(s) because no relevant samples were found in dataset '{1}'.".format(
                output_func_name, self._mongo_db_collection.full_name))

        return no_output


class VaprAnnotator(object):
    SAMPLE_NAMES_KEY = "Sample_Names"
    HG19_VERSION = "hg19"
    HG38_VERSION = "hg38"
    DEFAULT_GENOME_VERSION = HG19_VERSION
    SUPPORTED_GENOME_BUILD_VERSIONS = [HG19_VERSION, HG38_VERSION]

    @staticmethod
    def _get_num_lines_in_file(file_path):
        with open(file_path) as file_obj:
            result = sum(1 for _ in file_obj)
        return result

    @staticmethod
    def _make_jobs_params_tuples_list(file_path, num_file_lines, chunk_size, db_name, collection_name,
                                      genome_build_version, sample_names_list=None, verbose_level=1):

        num_params = VAPr.chunk_processing.AnnotationJobParamsIndices.get_num_possible_indices()
        if sample_names_list is not None:
            shared_job_params = [None] * num_params
            shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.SAMPLE_LIST_INDEX] = sample_names_list
        else:
            shared_job_params = [None] * (num_params - 1)

        shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.CHUNK_SIZE_INDEX] = chunk_size
        shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.FILE_PATH_INDEX] = file_path
        shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.DB_NAME_INDEX] = db_name
        shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.COLLECTION_NAME_INDEX] = collection_name
        shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.GENOME_BUILD_VERSION_INDEX] = \
            genome_build_version
        shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.VERBOSE_LEVEL_INDEX] = verbose_level

        jobs_params_tuples_list = []
        num_steps = int(num_file_lines / chunk_size) + 1
        for curr_chunk_index in range(num_steps):
            shared_job_params[VAPr.chunk_processing.AnnotationJobParamsIndices.CHUNK_INDEX_INDEX] = curr_chunk_index
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

    # TODO: Decide on what tests to write for top-level functions
    def __init__(self, input_dir, output_dir, mongo_db_name, mongo_collection_name, annovar_install_path=None,
                 design_file=None, build_ver=None, vcfs_gzipped=False):

        self._input_dir = input_dir
        self._output_dir = output_dir
        self._mongo_db_name = mongo_db_name
        self._mongo_collection_name = mongo_collection_name
        self._analysis_name = mongo_db_name + "_" + mongo_collection_name
        self._path_to_annovar_install = annovar_install_path
        self._design_file = design_file
        self._vcfs_gzipped = vcfs_gzipped

        self._genome_build_version = self._get_validated_genome_version(build_ver)

        self._single_vcf_path = self._make_merged_vcf_fp()
        self._output_basename = os.path.splitext(os.path.basename(self._single_vcf_path))[0]
        self._sample_names_list = vcf.Reader(open(self._single_vcf_path, 'r')).samples

        # TODO: someday: put back the functionality for custom annovar dbs?
        self._annovar_wrapper = None
        if self._path_to_annovar_install is not None:
            self._annovar_wrapper = VAPr.annovar_running.AnnovarWrapper(
                self._path_to_annovar_install, genome_build_version=self._genome_build_version,
                custom_annovar_dbs_to_use=None)

        try:
            os.mkdir(output_dir)
        except OSError:
            logging.info('Output directory %s for analysis already exists; using existing directory' % output_dir)

    def download_annovar_databases(self):
        """Run ANNOVAR to download its databases."""
        if self._path_to_annovar_install is None:
            raise ValueError("No ANNOVAR install path provided.")

        self._annovar_wrapper.download_databases()

    def annotate_lite(self, num_processes=8, chunk_size=2000, verbose_level=1, allow_adds=False):
        result = self._make_dataset_for_results("annotate_lite", allow_adds)
        self._collect_annotations_and_store(self._single_vcf_path, chunk_size, num_processes, sample_names_list=None,
                                            verbose_level=verbose_level)
        return result

    def annotate(self, num_processes=4, chunk_size=2000, verbose_level=1, allow_adds=False):
        if self._path_to_annovar_install is None:
            raise ValueError("No ANNOVAR install path provided.")

        result = self._make_dataset_for_results("annotate", allow_adds)
        annovar_output_fp = self._annovar_wrapper.run_annotation(self._single_vcf_path, self._output_basename,
                                                                 self._output_dir)
        self._collect_annotations_and_store(annovar_output_fp, chunk_size, num_processes,
                                            sample_names_list=self._sample_names_list, verbose_level=verbose_level)
        return result

    def _make_merged_vcf_fp(self):
        vcf_file_paths_list = None
        if self._design_file is not None:
            design_df = pandas.read_csv(self._design_file)
            vcf_file_paths_list = design_df[self.SAMPLE_NAMES_KEY].tolist()

        result = VAPr.vcf_merging.merge_vcfs(self._input_dir, self._output_dir, self._analysis_name,
                                             vcf_file_paths_list, self._vcfs_gzipped)
        return result

    def _make_dataset_for_results(self, func_name, allow_adds):
        result = VaprDataset(self._mongo_db_name, self._mongo_collection_name)

        distinct_ids_list = result.get_distinct_sample_ids()
        if len(distinct_ids_list) > 0:
            msg_prefix = "Dataset '{0}' already contains {1} records".format(result.full_name, len(distinct_ids_list))
            if allow_adds:
                logging.info("{0}; adding to this dataset.".format(msg_prefix))
            else:
                error_msg = "{0}, which is disallowed by default.  Either create a VaprAnnotator with a new " \
                            "collection name, clear your existing collection manually, or (if you definitely wish to " \
                            "add to an existing dataset), rerun {1} with the 'allow_adds' parameter set to " \
                            "True".format(msg_prefix, func_name)
                raise ValueError(error_msg)

        return result


    # TODO: someday: extra_data from design file needs to come back in here
    def _collect_annotations_and_store(self, file_path, chunk_size, num_processes, sample_names_list=None,
                                       verbose_level=1):
        num_file_lines = self._get_num_lines_in_file(file_path)
        jobs_params_tuples_list = self._make_jobs_params_tuples_list(
            file_path, num_file_lines, chunk_size, self._mongo_db_name, self._mongo_collection_name,
            self._genome_build_version, sample_names_list, verbose_level)

        pool = multiprocessing.Pool(num_processes)
        for _ in tqdm.tqdm(
                pool.imap_unordered(
                    VAPr.chunk_processing.collect_chunk_annotations_and_store, jobs_params_tuples_list),
                total=len(jobs_params_tuples_list)):
            pass
        pool.close()
        pool.join()
