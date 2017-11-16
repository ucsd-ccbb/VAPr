# standard libraries
import unittest

# project-specific libraries
import VAPr.filtering as ns_test


class TestFunctions(unittest.TestCase):
    def test_get_sample_id_filter(self):
        expected_output = {'samples.sample_id': "testname"}
        real_output = ns_test.get_sample_id_filter("testname")
        self.assertEqual(expected_output, real_output)

    def test_get_any_of_sample_ids_filter(self):
        expected_output = {'samples.sample_id': {'$in': ["testname1", "testname2"]}}
        real_output = ns_test.get_any_of_sample_ids_filter(["testname1", "testname2"])
        self.assertEqual(expected_output, real_output)

    def test_make_rare_deleterious_variants_filter_w_samples(self):
        expected_output = {
            "$and":
                [
                    {
                        "$or":
                            [
                                {"cadd.esp.af": {"$lt": 0.051}},
                                {"cadd.esp.af": {"$exists": False}}
                            ]
                    },
                    {
                        "$or":
                            [
                                {"func_knowngene": "exonic"},
                                {"func_knowngene": "splicing"}
                            ]
                    },
                    {"exonicfunc_knowngene": {"$ne": "synonymous SNV"}},
                    {"1000g2015aug_all": {"$lt": 0.051}},
                    {'samples.sample_id': {"$in":["testname1", "testname2"]}}

                ]
        }
        real_output = ns_test.make_rare_deleterious_variants_filter(["testname1", "testname2"])
        self.assertEqual(expected_output, real_output)

    def test_make_rare_deleterious_variants_filter_wo_samples(self):
        expected_output = {
            "$and":
                [
                    {
                        "$or":
                            [
                                {"cadd.esp.af": {"$lt": 0.051}},
                                {"cadd.esp.af": {"$exists": False}}
                            ]
                    },
                    {
                        "$or":
                            [
                                {"func_knowngene": "exonic"},
                                {"func_knowngene": "splicing"}
                            ]
                    },
                    {"exonicfunc_knowngene": {"$ne": "synonymous SNV"}},
                    {"1000g2015aug_all": {"$lt": 0.051}}

                ]
        }
        real_output = ns_test.make_rare_deleterious_variants_filter()
        self.assertEqual(expected_output, real_output)

    def test_make_known_disease_variants_filter_w_samples(self):
        expected_output = {
            "$and":
                [
                    {"$or":
                        [
                            {
                                "$and":
                                    [
                                        {"clinvar.rcv.accession": {"$exists": True}},
                                        {"clinvar.rcv.clinical_significance": {"$nin": ["Benign", "Likely benign"]}}
                                    ]
                            },
                            {"cosmic.cosmic_id": {"$exists": True}}
                        ]},
                    {'samples.sample_id': {"$in": ["testname1", "testname2"]}}
                ]
        }
        real_output = ns_test.make_known_disease_variants_filter(["testname1", "testname2"])
        self.assertEqual(expected_output, real_output)

    def test_make_known_disease_variants_filter_wo_samples(self):
        expected_output = {
            "$or":
                    [
                        {
                            "$and":
                                [
                                    {"clinvar.rcv.accession": {"$exists": True}},
                                    {"clinvar.rcv.clinical_significance": {"$nin": ["Benign", "Likely benign"]}}
                                ]
                        },
                        {"cosmic.cosmic_id": {"$exists": True}}
                    ]
        }
        real_output = ns_test.make_known_disease_variants_filter()
        self.assertEqual(expected_output, real_output)

    def test_make_deleterious_compound_heterozygote_variants_filter_w_samples(self):
        expected_output = {
            "$and":
                [
                    {"genotype_subclass_by_class.heterozygous": "compound"},
                    {"cadd.phred": {"$gte": 10}},
                    {'samples.sample_id': {"$in": ["testname1", "testname2"]}}
                ]
        }
        real_output = ns_test.make_deleterious_compound_heterozygous_variants_filter(["testname1", "testname2"])
        self.assertEqual(expected_output, real_output)

    def test_make_deleterious_compound_heterozygote_variants_filter_wo_samples(self):
        expected_output = {
            "$and":
                [
                    {"genotype_subclass_by_class.heterozygous": "compound"},
                    {"cadd.phred": {"$gte": 10}}
                ]
        }
        real_output = ns_test.make_deleterious_compound_heterozygote_variants_filter()
        self.assertEqual(expected_output, real_output)

    def test_make_de_novo_variants_filter(self):
        expected_output = {
            "$and":
                    [
                        {'samples.sample_id': "sampleA"},
                        {
                            "$and":
                                [
                                    {'samples.sample_id': {"$ne": "sampleB"}},
                                    {'samples.sample_id': {"$ne": "sampleC"}}
                                ]
                        }
                    ]
            }
        real_output = ns_test.make_de_novo_variants_filter("sampleA", "sampleB", "sampleC")
        self.assertEqual(expected_output, real_output)

    def test__append_sample_id_constraint_if_needed_is_needed(self):
        input_list = [
                        {"genotype_subclass_by_class.heterozygous": "compound"},
                        {"cadd.phred": {"$gte": 10}}
                     ]
        expected_output = {
            "$and":
                [
                    {"genotype_subclass_by_class.heterozygous": "compound"},
                    {"cadd.phred": {"$gte": 10}},
                    {'samples.sample_id': {"$in": ["testname1", "testname2"]}}
                ]
        }
        real_output = ns_test._append_sample_id_constraint_if_needed(input_list, ["testname1", "testname2"])
        self.assertDictEqual(expected_output, real_output)

    def test__append_sample_id_constraint_if_needed_is_not_needed(self):
        input_list = [
                        {"genotype_subclass_by_class.heterozygous": "compound"},
                        {"cadd.phred": {"$gte": 10}}
                     ]
        expected_output = {
            "$and":
                [
                    {"genotype_subclass_by_class.heterozygous": "compound"},
                    {"cadd.phred": {"$gte": 10}}
                ]
        }
        real_output = ns_test._append_sample_id_constraint_if_needed(input_list, None)
        self.assertDictEqual(expected_output, real_output)
