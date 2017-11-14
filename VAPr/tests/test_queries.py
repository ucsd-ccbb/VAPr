# project-specific libraries
from VAPr.queries import Filters
from pymongo import MongoClient
import unittest


__author__ = 'Adam Mark'


class TestQueries(unittest.TestCase):

    # This effectively created a test database with the documents in mongo required by the other tests,
    # not much else is needed to be tested here.

    def setUp(self):

        self.var1 = {"hgvs_id": "chr1:g.1000A>C",
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

        self.var2 = {"hgvs_id": "chr1:g.2000G>T",
                     "cadd":
                         {"esp" : {"af": 0.05}},
                     "func_knowngene": "intronic",
                     "1000g2015aug_all": 0.06,
                     "exonicfunc_knowngene": "nonsynonymous SNV",
                     "cosmic": {"cosmic_id": "XYZ789"},
                     "samples": {"sample_id": "sample2"}}

        self.var3 = {"hgvs_id": "chr1:g.3000T>A",
                     "cadd":
                         {"esp" : {"af": 0.95},
                          "phred": 40},
                     "func_knowngene": "exonic",
                     "1000g2015aug_all": 0.05,
                     "exonicfunc_knowngene": "synonymous SNV",
                     "genotype_subclass_by_class": {"heterozygous": "compound"},
                     "samples": {"sample_id": "sample3"}}

        self.client = MongoClient()
        self.db_name = "queries_test"
        self.db = getattr(self.client, self.db_name)
        self.collection_name = "collect"
        self.collection = getattr(self.db, self.collection_name)

    def test_variants_from_sample(self):
        self.collection.delete_many({})
        self.collection.insert_many([self.var1, self.var2, self.var3])
        filter_collection = Filters(self.db_name, self.collection_name)
        sample1_var = filter_collection.variants_from_sample("sample1")
        self.assertTrue(len(sample1_var) == 1)
        self.assertEqual(sample1_var[0]['hgvs_id'], self.var1['hgvs_id'])

    def test_variants_from_samples(self):
        self.collection.delete_many({})
        self.collection.insert_many([self.var1, self.var2, self.var3])
        filter_collection = Filters(self.db_name, self.collection_name)
        sample_var = filter_collection.variants_from_samples(["sample1", "sample2"])
        self.assertTrue(len(sample_var) == 2)
        self.assertListEqual([var['hgvs_id'] for var in sample_var], [self.var1['hgvs_id'], self.var2['hgvs_id']])

    def test_rare_deleterious_variants(self):
        self.collection.delete_many({})
        self.collection.insert_many([self.var1, self.var2, self.var3])
        filter_collection = Filters(self.db_name, self.collection_name)
        rdv = filter_collection.rare_deleterious_variants()
        self.assertEqual(rdv[0]['samples']['sample_id'], self.var1['samples']['sample_id'])

    def test_known_disease_variants(self):
        self.collection.delete_many({})
        self.collection.insert_many([self.var1, self.var2, self.var3])
        filter_collection = Filters(self.db_name, self.collection_name)
        kdv = filter_collection.known_disease_variants()
        self.assertListEqual([var['hgvs_id'] for var in kdv], [self.var1['hgvs_id'], self.var2['hgvs_id']])

    def test_deleterious_compound_heterozygote_variants(self):
        self.collection.delete_many({})
        self.collection.insert_many([self.var1, self.var2, self.var3])
        filter_collection = Filters(self.db_name, self.collection_name)
        dch = filter_collection.deleterious_compound_heterozygote_variants()
        self.assertListEqual([var['hgvs_id'] for var in dch], [self.var3['hgvs_id']])

    def test_de_novo_variants(self):
        self.collection.delete_many({})
        self.collection.insert_many([self.var1, self.var2, self.var3])
        filter_collection = Filters(self.db_name, self.collection_name)
        dnv = filter_collection.de_novo_variants("sample1", "sample2", "sample3")
        self.assertListEqual([var['hgvs_id'] for var in dnv], [self.var1['hgvs_id']])




