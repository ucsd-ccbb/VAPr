import unittest
import sys
import os
#quick and dirty way of importing functions
from variantannotation import csv_to_df

sys.path.append('/Users/carlomazzaferro/Documents/Code/variant-annotation/variantannotation')
file_name = os.path.dirname(os.path.realpath('__file__')) + '/test_file.csv'
sample_list = csv_to_df.open_and_parse(file_name)
sample_df = csv_to_df.parse_to_df(sample_list)


class MongoDBQueryTest(unittest.TestCase):

    def setUp(self):
        self.data = file_name