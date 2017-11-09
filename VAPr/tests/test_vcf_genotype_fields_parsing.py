# standard libraries
import unittest
import warnings

# project-specific libraries
import VAPr.vcf_genotype_fields_parsing as ns_test

__author__ = 'Birmingham'

# Cause all warnings to always be triggered.
warnings.simplefilter("always")


def _help_get_warn_msg(warn_obj):
    return str(warn_obj[-1].message)


class TestFunctions(unittest.TestCase):
    # region _fill_genotype tests
    def test_fill_genotype(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        ns_test._fill_genotype('0/2', genotype_to_fill)
        self.assertEqual('0/2', genotype_to_fill.genotype)

    def test_fill_genotype_warn(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        with warnings.catch_warnings(record=True) as w:
            genotype_to_fill = ns_test._fill_genotype('1/0/2', genotype_to_fill)
            self.assertEqual("The GT tag value 1/0/2 does not split into exactly two values so genotype information "
                             "could not be captured for the current variant.", _help_get_warn_msg(w))
            self.assertIsNone(genotype_to_fill.genotype)  # no genotype created if GT tag splits wrong

    # endregion

    # region _fill_unfiltered_reads_counts tests
    def test_fill_unfiltered_reads_counts(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        ns_test._fill_unfiltered_reads_counts('0,64,12', genotype_to_fill)
        self.assertEqual(3, len(genotype_to_fill.alleles))
        self.assertEqual(0, genotype_to_fill.alleles[0].unfiltered_read_counts)
        self.assertEqual(64, genotype_to_fill.alleles[1].unfiltered_read_counts)
        self.assertEqual(12, genotype_to_fill.alleles[2].unfiltered_read_counts)

    def test_fill_unfiltered_reads_counts_warn(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        with warnings.catch_warnings(record=True) as w:
            genotype_to_fill = ns_test._fill_unfiltered_reads_counts('0;64', genotype_to_fill)
            self.assertEqual("The AD tag value 0;64 does not split into at least two values so unfiltered allele depth"
                             " information could not be captured for the current variant.", _help_get_warn_msg(w))
            self.assertEqual(0, len(genotype_to_fill.alleles))  # no alleles created if AD tag splits wrong

    # endregion

    # region fill_filtered_reads tests
    def test_fill_filtered_reads_count(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        ns_test._fill_filtered_reads_count('44', genotype_to_fill)
        self.assertEqual(44, genotype_to_fill.filter_passing_reads_count)

    # endregion

    # region _fill_genotype_confidence tests
    def test_fill_genotype_confidence(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        ns_test._fill_genotype_confidence('44.1', genotype_to_fill)
        self.assertEqual(44.1, genotype_to_fill.genotype_confidence)

    # endregion

    # region _fill_genotype_likelihoods tests
    def test_fill_genotype_likelihoods_three_alleles(self):
        expected_values = [(0, 0, 495),
                           (0, 1, 162),
                           (1, 1, 123),
                           (0, 2, 213),
                           (1, 2, 129),
                           (2, 2, 175),
                           (0, 3, 67),
                           (1, 3, 0),
                           (2, 3, 46),
                           (3, 3, 28.1)]
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [ns_test.Allele(10), ns_test.Allele(11), ns_test.Allele(12),
                                    ns_test.Allele(13)]
        ns_test._fill_genotype_likelihoods('495,162,123,213,129,175,67,0,46,28.1', genotype_to_fill)
        self.assertEqual(len(expected_values), len(genotype_to_fill.genotype_likelihoods))
        for index in range(0, len(expected_values)):
            curr_expected_values = expected_values[index]
            real_values = genotype_to_fill.genotype_likelihoods[index]
            self.assertEqual(curr_expected_values[0], real_values.allele1_number)
            self.assertEqual(curr_expected_values[1], real_values.allele2_number)
            self.assertEqual(curr_expected_values[2], real_values.likelihood_neg_exponent)

    def test_fill_genotype_likelihoods_no_alleles(self):
        expected_values = [(0, 0, 495),
                           (0, 1, 162),
                           (1, 1, 123)]
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        self.assertEqual(0, len(genotype_to_fill.alleles))
        ns_test._fill_genotype_likelihoods('495,162,123', genotype_to_fill)
        self.assertEqual(2, len(genotype_to_fill.alleles))  # should have added two
        self.assertEqual(len(expected_values), len(genotype_to_fill.genotype_likelihoods))
        for index in range(0, len(expected_values)):
            curr_expected_values = expected_values[index]
            real_values = genotype_to_fill.genotype_likelihoods[index]
            self.assertEqual(curr_expected_values[0], real_values.allele1_number)
            self.assertEqual(curr_expected_values[1], real_values.allele2_number)
            self.assertEqual(curr_expected_values[2], real_values.likelihood_neg_exponent)

    def test_fill_genotype_likelihoods_warn_too_many_alleles_implied(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [ns_test.Allele(10), ns_test.Allele(11), ns_test.Allele(12)]
        with warnings.catch_warnings(record=True) as w:
            genotype_to_fill = ns_test._fill_genotype_likelihoods('495,162,123,213,129,175,67,0,46,28.1',
                                                                  genotype_to_fill)
            self.assertEqual("The PL tag value 495,162,123,213,129,175,67,0,46,28.1 appears to contain information for "
                             "more alleles than expected so 'normalized' Phred-scaled likelihoods of possible genotypes"
                             " information could not be captured for the current variant.", _help_get_warn_msg(w))
            self.assertEqual(0, len(genotype_to_fill.genotype_likelihoods))  # no likelihoods made if PL splits wrong

    def test_fill_genotype_likelihoods_warn_too_few_likelihoods_1(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [ns_test.Allele(10), ns_test.Allele(11), ns_test.Allele(12),
                                    ns_test.Allele(13)]
        with warnings.catch_warnings(record=True) as w:
            genotype_to_fill = ns_test._fill_genotype_likelihoods('495,162,123,213,129,175', genotype_to_fill)
            self.assertEqual("The PL tag value 495,162,123,213,129,175 appears to contain information for fewer alleles"
                             " than expected so 'normalized' Phred-scaled likelihoods of possible genotypes information"
                             " could not be captured for the current variant.", _help_get_warn_msg(w))
            self.assertEqual(0, len(genotype_to_fill.genotype_likelihoods))  # no likelihoods made if PL splits wrong

    def test_fill_genotype_likelihoods_warn_too_few_likelihoods_2(self):
        genotype_to_fill = ns_test.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [ns_test.Allele(10), ns_test.Allele(11), ns_test.Allele(12),
                                    ns_test.Allele(13)]
        with warnings.catch_warnings(record=True) as w:
            genotype_to_fill = ns_test._fill_genotype_likelihoods('495,162,123,213,129,175,67,0,46', genotype_to_fill)
            self.assertEqual("The PL tag value 495,162,123,213,129,175,67,0,46 appears to contain information for fewer"
                             " alleles than expected so 'normalized' Phred-scaled likelihoods of possible genotypes "
                             "information could not be captured for the current variant.", _help_get_warn_msg(w))
            self.assertEqual(0, len(genotype_to_fill.genotype_likelihoods))  # no likelihoods made if PL splits wrong

    # endregion


class TestVCFGenotypeInfo(unittest.TestCase):
    # No tests of __init__ as it is just setting empty values
    # No explicit tests of getter properties as they're just returning an internal variable

    # region contains_no_genotype_call tests
    def test_contains_no_genotype_call_true(self):
        dummy_vcfgenotypeinfo = ns_test.VCFGenotypeInfo('./.:whatevs')
        self.assertTrue(dummy_vcfgenotypeinfo.contains_no_genotype_call)

    def test_contains_no_genotype_call_false(self):
        dummy_vcfgenotypeinfo = ns_test.VCFGenotypeInfo('1/1:0,34:34:99:1187.2,101,0')
        self.assertFalse(dummy_vcfgenotypeinfo.contains_no_genotype_call)

    # endregion

    # region genotype_confidence setter tests
    def test_genotype_confidence_setter(self):
        dummy_vcfgenotypeinfo = ns_test.VCFGenotypeInfo('')
        dummy_vcfgenotypeinfo.genotype_confidence = "89"
        self.assertEqual(89.0, dummy_vcfgenotypeinfo.genotype_confidence)
        dummy_vcfgenotypeinfo.genotype_confidence = 89
        self.assertEqual(89.0, dummy_vcfgenotypeinfo.genotype_confidence)
        dummy_vcfgenotypeinfo.genotype_confidence = "-89.10"
        self.assertEqual(-89.1, dummy_vcfgenotypeinfo.genotype_confidence)

    def test_genotype_confidence_setter_error(self):
        dummy_vcfgenotypeinfo = ns_test.VCFGenotypeInfo('')
        with self.assertRaises(ValueError):
            dummy_vcfgenotypeinfo.genotype_confidence = "blue"

    # endregion

    # region filter_passing_reads_count setter tests
    def test_filter_passing_reads_count_setter(self):
        dummy_vcfgenotypeinfo = ns_test.VCFGenotypeInfo('')
        dummy_vcfgenotypeinfo.filter_passing_reads_count = "42"
        self.assertEqual(42, dummy_vcfgenotypeinfo.filter_passing_reads_count)
        dummy_vcfgenotypeinfo.filter_passing_reads_count = 0
        self.assertEqual(0, dummy_vcfgenotypeinfo.filter_passing_reads_count)
        dummy_vcfgenotypeinfo.filter_passing_reads_count = "."
        self.assertEqual(None, dummy_vcfgenotypeinfo.filter_passing_reads_count)

    def test_filter_passing_reads_count_setter_error(self):
        dummy_vcfgenotypeinfo = ns_test.VCFGenotypeInfo('')
        with self.assertRaises(ValueError):
            dummy_vcfgenotypeinfo.filter_passing_reads_count = -1
        with self.assertRaises(ValueError):
            dummy_vcfgenotypeinfo.filter_passing_reads_count = 48.5

    # endregion


class TestAllele(unittest.TestCase):
    # No tests of __init__ as it is just setting values
    # No explicit tests of getter properties as they're just returning an internal variable

    def test_read_counts_setter(self):
        dummy_allele = ns_test.Allele(0)
        self.assertEqual(0, dummy_allele.unfiltered_read_counts)
        dummy_allele.unfiltered_read_counts = 94
        self.assertEqual(94, dummy_allele.unfiltered_read_counts)

    def test_read_counts_setter_error(self):
        dummy_allele = ns_test.Allele(0)
        with self.assertRaises(ValueError):
            dummy_allele.unfiltered_read_counts = -94
        with self.assertRaises(ValueError):
            dummy_allele.unfiltered_read_counts = 48.5


class TestGenotypeLikelihood(unittest.TestCase):
    # No tests of init as just calls tested setters
    # No explicit tests of getter properties as they're just returning an internal variable

    # region _validate_allele_relationship tests
    def test__validate_allele_relationship_pass(self):
        ns_test.GenotypeLikelihood._validate_allele_relationship(0, 2)
        # if we got this far, the test passed

    def test__validate_allele_relationship_fail(self):
        with self.assertRaises(ValueError):
            ns_test.GenotypeLikelihood._validate_allele_relationship(2, 0)

    # endregion

    # region allele1_number setter tests
    def test_allele1_number_setter(self):
        temp_likelihood = ns_test.GenotypeLikelihood(0, 5, 0)
        temp_likelihood.allele1_number = 4
        self.assertEqual(4, temp_likelihood.allele1_number)

    def test_allele1_number_setter_error(self):
        with self.assertRaises(ValueError):
            ns_test.GenotypeLikelihood(0, 0, 0).allele1_number = -1
        with self.assertRaises(ValueError):
            ns_test.GenotypeLikelihood(0, 5, 0).allele1_number = 4.5
        with self.assertRaises(ValueError):
            ns_test.GenotypeLikelihood(0, 0, 0).allele1_number = "blue"

        temp_likelihood = ns_test.GenotypeLikelihood(0, 0, 0)
        temp_likelihood.allele2_number = 1
        with self.assertRaises(ValueError):
            temp_likelihood.allele1_number = 4  # can't be greater than allele 2

    # endregion

    # region allele2_number setter tests
    def test_allele2_number_setter(self):
        temp_likelihood = ns_test.GenotypeLikelihood(0, 0, 0)
        temp_likelihood.allele2_number = 4
        self.assertEqual(4, temp_likelihood.allele2_number)

    def test_allele2_number_setter_error(self):
        with self.assertRaises(ValueError):
            ns_test.GenotypeLikelihood(0, 0, 0).allele2_number = -1
        with self.assertRaises(ValueError):
            ns_test.GenotypeLikelihood(0, 0, 0).allele2_number = 4.5
        with self.assertRaises(ValueError):
            ns_test.GenotypeLikelihood(0, 0, 0).allele2_number = "blue"

        temp_likelihood = ns_test.GenotypeLikelihood(4, 4, 0)
        with self.assertRaises(ValueError):
            temp_likelihood.allele2_number = 2  # can't be smaller than allele 1

    # endregion

    # region likelihood_neg_exponent setter tests
    def test_likelihood_neg_exponent_setter(self):
        temp_likelihood = ns_test.GenotypeLikelihood(0, 0, 0)
        temp_likelihood.likelihood_neg_exponent = "89"
        self.assertEqual(89.0, temp_likelihood.likelihood_neg_exponent)
        temp_likelihood.likelihood_neg_exponent = 89
        self.assertEqual(89.0, temp_likelihood.likelihood_neg_exponent)
        temp_likelihood.likelihood_neg_exponent = "-89.10"
        self.assertEqual(-89.1, temp_likelihood.likelihood_neg_exponent)

    def test_genotype_confidence_setter_error(self):
        temp_likelihood = ns_test.GenotypeLikelihood(0, 0, 0)
        with self.assertRaises(ValueError):
            temp_likelihood.likelihood_neg_exponent = "blue"

    # endregion


class TestVCFGenotypeString(unittest.TestCase):
    def test_parse_GT_GQ_PL(self):
        format_string = 'GT:GQ:PL'
        info_string = '1/1:99:1187.2,101,0'
        parser = ns_test.VCFGenotypeParser()
        genotype_to_fill = parser.parse(format_string, info_string)
        self.assertEqual('1/1', genotype_to_fill.genotype)
        self.assertIsNone(genotype_to_fill.filter_passing_reads_count)
        self.assertEqual(99, genotype_to_fill.genotype_confidence)
        self.assertEqual(2, len(genotype_to_fill.alleles))
        self.assertEqual(3, len(genotype_to_fill.genotype_likelihoods))
        self.assertEqual(1187.2, genotype_to_fill.genotype_likelihoods[0].likelihood_neg_exponent)
        self.assertEqual(101, genotype_to_fill.genotype_likelihoods[1].likelihood_neg_exponent)
        self.assertEqual(0, genotype_to_fill.genotype_likelihoods[2].likelihood_neg_exponent)

    def test_parse_GT_AD_DP_GQ_PL(self):
        format_string = 'GT:AD:DP:GQ:PL'
        info_string = '1/1:0,34:34:99:1187.2,101,0'
        parser = ns_test.VCFGenotypeParser()
        genotype_to_fill = parser.parse(format_string, info_string)

        self.assertEqual('1/1', genotype_to_fill.genotype)
        self.assertEqual(34, genotype_to_fill.filter_passing_reads_count)
        self.assertEqual(99, genotype_to_fill.genotype_confidence)
        self.assertEqual(2, len(genotype_to_fill.alleles))
        self.assertEqual(3, len(genotype_to_fill.genotype_likelihoods))
        self.assertEqual(0, genotype_to_fill.alleles[0].unfiltered_read_counts)
        self.assertEqual(1187.2, genotype_to_fill.genotype_likelihoods[0].likelihood_neg_exponent)
        self.assertEqual(34, genotype_to_fill.alleles[1].unfiltered_read_counts)
        self.assertEqual(101, genotype_to_fill.genotype_likelihoods[1].likelihood_neg_exponent)
        self.assertEqual(0, genotype_to_fill.genotype_likelihoods[2].likelihood_neg_exponent)

    def test_parse_error_reported_as_warn(self):
        format_string = 'GT:AD:GQ:PL'
        info_string = 'blue'
        parser = ns_test.VCFGenotypeParser()
        with warnings.catch_warnings(record=True) as w:
            genotype_to_fill = parser.parse(format_string, info_string)
            self.assertEqual("Encountered error 'list index out of range' so genotype fields information could not be "
                             "captured for the current variant.", _help_get_warn_msg(w))
            self.assertIsNone(genotype_to_fill)  # warn and return None but don't error out if any one parse fails
