# standard libraries
import tempfile
import unittest
import warnings

# project-specific libraries
import VAPr.vapr_core as ns_test


class TestVaprDataset(unittest.TestCase):
    # This effectively created a test database with the documents in mongo required by the other tests,
    # not much else is needed to be tested here.

    @classmethod
    def setUpClass(cls):
        cls.var1 = {"hgvs_id": "chr1:g.1000A>C",
                     "cadd":
                         {"esp" : {"af": 0.05},
                          "phred": 11},
                     "func_knowngene": "exonic",
                     "1000g2015aug_all": 0.05,
                     "exonicfunc_knowngene": "nonsynonymous SNV",
                     "clinvar" :
                         {"rcv":
                              {"accession": "ABC123",
                               "clinical_significance": "Pathogenic"}},
                     "cosmic" : {"cosmic_id": "XYZ789"},
                     "samples": {"sample_id": "sample1"}}

        cls.var2 = {"hgvs_id": "chr1:g.2000G>T",
                     "cadd":
                         {"esp" : {"af": 0.05}},
                     "func_knowngene": "intronic",
                     "1000g2015aug_all": 0.06,
                     "exonicfunc_knowngene": "nonsynonymous SNV",
                     "cosmic": {"cosmic_id": "XYZ789"},
                     "samples": {"sample_id": "sample2"}}

        cls.var3 = {"hgvs_id": "chr1:g.3000T>A",
                     "cadd":
                         {"esp" : {"af": 0.95},
                          "phred": 40},
                     "func_knowngene": "exonic",
                     "1000g2015aug_all": 0.05,
                     "exonicfunc_knowngene": "synonymous SNV",
                     "genotype_subclass_by_class": {"heterozygous": "compound"},
                     "samples": {"sample_id": "sample3"}}

        cls._db_name = "queries_test"
        cls._collection_name = "collect"

    def test__write_annotated_vcf(self):
        self.fail("test not implemented")

    def test__get_filtered_variants_by_sample(self):
        self.fail("test not implemented")

    def test__write_annotated_csv(self):
        self.fail("test not implemented")

    def test__warn_if_no_output_true(self):
        self.fail("test not implemented")

    def test__warn_if_no_output_false(self):
        self.fail("test not implemented")

    def test_de_novo_variants(self):
        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        dnv = test_dataset.get_de_novo_variants("sample1", "sample2", "sample3")
        self.assertListEqual([var['hgvs_id'] for var in dnv], [self.var1['hgvs_id']])

    def test_deleterious_compound_heterozygote_variants_all_samples(self):
        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        dch = test_dataset.get_deleterious_compound_heterozygote_variants()
        self.assertListEqual([var['hgvs_id'] for var in dch], [self.var3['hgvs_id']])

    def test_deleterious_compound_heterozygote_variants_specific_samples(self):
        self.fail("test not implemented")

        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        dch = test_dataset.get_deleterious_compound_heterozygote_variants()
        self.assertListEqual([var['hgvs_id'] for var in dch], [self.var3['hgvs_id']])

    def test_known_disease_variants_all_samples(self):
        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        kdv = test_dataset.get_known_disease_variants()
        self.assertListEqual([var['hgvs_id'] for var in kdv], [self.var1['hgvs_id'], self.var2['hgvs_id']])

    def test_known_disease_variants_specific_samples(self):
        self.fail("test not implemented")

        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        kdv = test_dataset.get_known_disease_variants()
        self.assertListEqual([var['hgvs_id'] for var in kdv], [self.var1['hgvs_id'], self.var2['hgvs_id']])

    def test_rare_deleterious_variants_all_samples(self):
        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        rdv = test_dataset.get_rare_deleterious_variants()
        self.assertEqual(rdv[0]['samples']['sample_id'], self.var1['samples']['sample_id'])

    def test_rare_deleterious_variants_specific_samples(self):
        self.fail("test not implemented")

        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        rdv = test_dataset.get_rare_deleterious_variants()
        self.assertEqual(rdv[0]['samples']['sample_id'], self.var1['samples']['sample_id'])

    def test_get_custom_filtered_variants(self):
        self.fail("test not implemented")

    def test_get_distinct_sample_ids(self):
        self.fail("test not implemented")

    def test_get_all_variants(self):
        self.fail("test not implemented")

    def test_get_variants_for_sample(self):
        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        sample1_var = test_dataset.get_variants_for_sample("sample1")
        self.assertTrue(len(sample1_var) == 1)
        self.assertEqual(sample1_var[0]['hgvs_id'], self.var1['hgvs_id'])

    def test_get_variants_for_samples(self):
        test_dataset = ns_test.VaprDataset(self._db_name, self._collection_name)
        test_dataset._mongo_db_collection.delete_many({})

        test_dataset._mongo_db_collection.insert_many([self.var1, self.var2, self.var3])
        sample_var = test_dataset.get_variants_for_samples(["sample1", "sample2"])
        self.assertTrue(len(sample_var) == 2)
        self.assertListEqual([var['hgvs_id'] for var in sample_var], [self.var1['hgvs_id'], self.var2['hgvs_id']])

    def test_get_variants_as_dataframe(self):
        self.fail("test not implemented")

    def test_write_unfiltered_annotated_csv(self):
        self.fail("test not implemented")

    def test_write_filtered_annotated_csv(self):
        self.fail("test not implemented")

    def test_write_unfiltered_annotated_vcf(self):
        self.fail("test not implemented")

    def test_write_filtered_annotated_vcf(self):
        self.fail("test not implemented")

    def test_write_unfiltered_annotated_csvs_per_sample(self):
        self.fail("test not implemented")



class TestVaprAnnotator(unittest.TestCase):
    def test__get_num_lines_in_file(self):
        # create a temporary file with 10,000 lines and ensure that is how many lines we get back
        num_lines = 10000
        temp_file = tempfile.NamedTemporaryFile(delete=False)
        for _ in range(0, num_lines):
            temp_file.write(b'test line\n')
        temp_file.close()  # but don't delete yet, as delete=False

        real_output = ns_test.VaprAnnotator._get_num_lines_in_file(temp_file.name)
        self.assertEqual(num_lines, real_output)

    # region _make_jobs_params_tuples_list tests
    def test__make_jobs_params_tuples_list_no_samples_default_verbose(self):
        input_file_path = "my/path/to/file.txt"
        input_num_file_lines = 21
        input_chunk_size = 10
        input_db_name = "mydb"
        input_collection_name = "mycol"
        input_build_version = "hg19"
        default_verbose_level = 1

        expected_output = [(0, input_file_path, input_chunk_size, input_db_name, input_collection_name,
                            input_build_version, default_verbose_level),
                           (1, input_file_path, input_chunk_size, input_db_name, input_collection_name,
                            input_build_version, default_verbose_level),
                           (2, input_file_path, input_chunk_size, input_db_name, input_collection_name,
                            input_build_version, default_verbose_level)]

        real_output = ns_test.VaprAnnotator._make_jobs_params_tuples_list(
            input_file_path, input_num_file_lines, input_chunk_size, input_db_name, input_collection_name,
            input_build_version)

        self.assertListEqual(expected_output, real_output)

    def test__make_jobs_params_tuples_list_no_samples_default_verbose_less_than_chunk(self):
        input_file_path = "my/path/to/file.txt"
        input_num_file_lines = 2
        input_chunk_size = 10
        input_db_name = "mydb"
        input_collection_name = "mycol"
        input_build_version = "hg19"
        default_verbose_level = 1

        expected_output = [(0, input_file_path, input_chunk_size, input_db_name, input_collection_name,
                            input_build_version, default_verbose_level)]

        real_output = ns_test.VaprAnnotator._make_jobs_params_tuples_list(
            input_file_path, input_num_file_lines, input_chunk_size, input_db_name, input_collection_name,
            input_build_version)

        self.assertListEqual(expected_output, real_output)

    def test__make_jobs_params_tuples_list_samples_with_verbose(self):
        input_file_path = "my/path/to/file.txt"
        input_num_file_lines = 21
        input_chunk_size = 10
        input_db_name = "mydb"
        input_collection_name = "mycol"
        input_build_version = "hg19"
        input_verbose_level = 2
        input_sample_names_list = ["sample_1", "sample_2"]

        expected_output = [(0, input_file_path, input_chunk_size, input_db_name, input_collection_name,
                            input_build_version, input_verbose_level, input_sample_names_list),
                           (1, input_file_path, input_chunk_size, input_db_name, input_collection_name,
                            input_build_version, input_verbose_level, input_sample_names_list),
                           (2, input_file_path, input_chunk_size, input_db_name, input_collection_name,
                            input_build_version, input_verbose_level, input_sample_names_list)]

        real_output = ns_test.VaprAnnotator._make_jobs_params_tuples_list(
            input_file_path, input_num_file_lines, input_chunk_size, input_db_name, input_collection_name,
            input_build_version, sample_names_list=input_sample_names_list, verbose_level=input_verbose_level)

        self.assertListEqual(expected_output, real_output)

    # endregion

    # region _get_validated_genome_version tests
    def test__get_validated_genome_version_default(self):
        real_output = ns_test.VaprAnnotator._get_validated_genome_version(None)
        self.assertEqual(ns_test.VaprAnnotator.DEFAULT_GENOME_VERSION, real_output)

    def test__get_validated_genome_version_error(self):
        with self.assertRaises(ValueError):
            ns_test.VaprAnnotator._get_validated_genome_version("blue")

    def test__get_validated_genome_version(self):
        real_output = ns_test.VaprAnnotator._get_validated_genome_version(ns_test.VaprAnnotator.HG38_VERSION)
        self.assertEqual(ns_test.VaprAnnotator.HG38_VERSION, real_output)

    # endregion

    def test__make_merged_vcf_fp(self):
        self.fail("test not implemented")

    def test__make_dataset_for_results(self):
        self.fail("test not implemented")

    def test__collect_annotations_and_store(self):
        self.fail("test not implemented")

    def test_annotate(self):
        self.fail("test not implemented")

    def test_annotate_lite(self):
        self.fail("test not implemented")

    def test_download_annovar_databases(self):
        self.fail("test not implemented")

    def test___init__(self):
        self.fail("test not implemented")
