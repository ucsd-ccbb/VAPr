# project-specific libraries
from VAPr.annotation_project import AnnotationProject
from VAPr.queries import Filters
import logging
import unittest
import os
logger = logging.getLogger()
logger.setLevel(logging.INFO)


__author__ = 'Mazzaferro'


class TestQueries(unittest.TestCase):

    def setUp(self):

        self.base_dir = '/Volumes/Carlo_HD1/CCBB/VAPr_files/'
        self.input_dir = os.path.join(self.base_dir, 'vcf_files/multi_sample')
        self.design_file_dirs = os.path.join(self.base_dir, 'design_file_three_samples.csv')
        self.out_path = os.path.join(self.base_dir, 'csv_multisample')
        self.annovar = os.path.join(self.base_dir, '../annovar')
        self.project_data = {'db_name': 'VariantDatabase',
                             'collection_name': 'collect'}

        self.project = AnnotationProject(self.input_dir,
                                         self.out_path,
                                         self.annovar,
                                         self.project_data,
                                         design_file=self.design_file_dirs,
                                         build_ver='hg19',
                                         mongod_cmd='/usr/local/bin/mongod --dbpath /Volumes/Carlo_HD1/data/db/ '
                                                    '--storageEngine wiredTiger')

    def test_query_one(self):
        self.project.parallel_annotation_and_saving(n_processes=6)
        """

        :return:

        filt = Filters(self.project_data['db_name'], self.project_data['collection_name'])
        filtered = filt.rare_cancer_variant(samples=['RAND1', 'RAND'])
        self.assertEqual(len(filtered), 20)
        filtered = filt.rare_cancer_variant(samples=['RAND1'])
        self.assertEqual(len(filtered), 13)
        filtered = filt.rare_cancer_variant(samples=['RAND'])
        self.assertEqual(len(filtered), 7)
        """