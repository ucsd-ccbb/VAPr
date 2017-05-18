import unittest
import sys
import os
#quick and dirty way of importing functions


sys.path.append('/Users/carlomazzaferro/Documents/Code/variant-annotation/variantannotation')
file_name = os.path.dirname(os.path.realpath('__file__')) + '/test_file.csv'


class MongoDBQueryTest(unittest.TestCase):

    def setUp(self):
        self.data = file_name