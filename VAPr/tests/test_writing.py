# project-specific libraries

from VAPr.writes import Writer
from VAPr.queries import Filters
import logging
import unittest
import os
from pymongo import MongoClient
logger = logging.getLogger()
logger.setLevel(logging.INFO)


__author__ = 'Mazzaferro'


class TestWrites(unittest.TestCase):

    """ These tests work when the test_parallel_annotation script has been run twice, where the first
    run is done with the self.minimapper set to (['RAND'], ...) and the second with
    self.minimapper set to (['RAND1'], ...) . The sample names RAND and RAND1 are parsed to
    mongo """

    def setUp(self):

        self.base_dir = os.getcwd()
        self.fwriter_output_dir = os.path.join(self.base_dir, 'test_files/fwriter_output')
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}

    def test_query_and_write_rand1(self):
        fwriter = Writer(self.project_data['collection_name'], self.project_data['db_name'])
        filt = Filters(self.project_data['db_name'], self.project_data['collection_name'])
        filtered = filt.rare_cancer_variant(samples=['RAND'])
        self.assertEqual(len(filtered), 0)

        # Test a stub function that is more lenient
        filtered = self.rare_cancer_variant_stub(samples=['RAND'])
        self.assertEqual(len(filtered), 2)
        fname = os.path.join(self.fwriter_output_dir, 'annotated_csv_rand1.csv')
        fwriter.generate_annotated_csv(filtered, fname)
        row_count = sum(1 for _ in open(fname, 'r'))
        self.assertEqual(row_count, len(filtered) + 1)  # to account for the header

    def test_and_query_rand(self):
        fwriter = Writer(self.project_data['db_name'], self.project_data['collection_name'])
        filt = Filters(self.project_data['db_name'], self.project_data['collection_name'])
        filtered = filt.rare_cancer_variant(samples=['RAND'])
        self.assertEqual(len(filtered), 0)

        # Test a stub function that is more lenient
        filtered = self.rare_cancer_variant_stub(samples=['RAND'])
        self.assertEqual(len(filtered), 2)
        fname = os.path.join(self.fwriter_output_dir, 'annotated_csv_rand1.csv')
        fwriter.generate_annotated_csv(filtered, fname)
        row_count = sum(1 for _ in open(fname, 'r'))
        self.assertEqual(row_count, len(filtered) + 1)  # to account for the header

    def test_query_by_sample(self):
        """ Tests the function that generates files for each sample, implemented in
        VAPr.parsers.generate_output_files_by_sample """

        samples = ['RAND', 'RAND1']
        fwriter = Writer(self.project_data['db_name'], self.project_data['collection_name'])
        filt = Filters(self.project_data['db_name'], self.project_data['collection_name'])

        for sample in samples:
            q = filt.variants_from_sample(sample)
            list_docs = list(q)
            for doc in list_docs:
                self.assertEqual(doc['sample_id'], sample)
            fname = os.path.join(self.fwriter_output_dir, 'annotated_csv_' + sample + '_all_vars.csv')
            fwriter.generate_annotated_csv(list_docs, fname)
            row_count = sum(1 for _ in open(fname, 'r'))
            self.assertEqual(row_count, 16000 + 1)

    def rare_cancer_variant_stub(self, samples):
        """ Function for retrieving rare cancer variants, just a bit more lenient to make
         the tests easier to debug, and needing less variants in the stub database created
         by the test_parallel_annotation_test """

        client = MongoClient()
        db = getattr(client, self.project_data['db_name'])
        collection = getattr(db, self.project_data['collection_name'])
        if not samples:
            samples = collection.distinct('sample_id')

        filtered = collection.find(
            {
                "$and":
                    [
                        {
                            "$or":
                                [
                                    {"esp6500siv2_all": {"$lt": 0.05}},
                                    {"esp6500siv2_all": {"$exists": False}}
                                ]
                        },
                        {
                            "$or":
                                [
                                    {"func_knowngene": "exonic"},
                                    {"func_knowngene": "splicing"}
                                ]
                        },
                        {"1000g2015aug_all": {"$lt": 100}},
                        {"sample_id": {"$in": samples}}
                    ]
            }
        )

        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered