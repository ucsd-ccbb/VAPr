# standard libraries
import unittest

# personal libraries
import utilities.database

# project-specific libraries
import vcf_parsing

__author__ = 'Birmingham'


class TestFunctions(unittest.TestCase):
    # No test for ignore_pid because it just logs a string
    # No test for ignore_pgt because it just logs a string
    # No test for ignore_field because it just logs a string

    # region fill_genotype tests
    def test_fill_genotype(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        vcf_parsing.fill_genotype('0/2', genotype_to_fill)
        self.assertEqual('0/2',genotype_to_fill.genotype)

    def test_fill_genotype_error(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        with self.assertRaises(ValueError):
            vcf_parsing.fill_genotype('1/0/2', genotype_to_fill)
    # endregion

    # region fill_unfiltered_reads_counts tests
    def test_fill_unfiltered_reads_counts(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        vcf_parsing.fill_unfiltered_reads_counts('0,64,12', genotype_to_fill)
        self.assertEqual(3, len(genotype_to_fill.alleles))
        self.assertEqual(0, genotype_to_fill.alleles[0].read_counts)
        self.assertEqual(64, genotype_to_fill.alleles[1].read_counts)
        self.assertEqual(12, genotype_to_fill.alleles[2].read_counts)

    def test_fill_unfiltered_reads_counts_error(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        with self.assertRaises(ValueError):
            vcf_parsing.fill_unfiltered_reads_counts('0;64', genotype_to_fill)
    # endregion

    # region fill_filtered_reads tests
    def test_fill_filtered_reads_count(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        vcf_parsing.fill_filtered_reads_count('44', genotype_to_fill)
        self.assertEqual(44,genotype_to_fill.filter_passing_reads_count)
    # endregion

    # region fill_genotype_confidence tests
    def test_fill_genotype_confidence(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        vcf_parsing.fill_genotype_confidence('44.1', genotype_to_fill)
        self.assertEqual(44.1,genotype_to_fill.genotype_confidence)
    # endregion

    # region fill_genotype_likelihoods tests
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
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [vcf_parsing.Allele(10), vcf_parsing.Allele(11), vcf_parsing.Allele(12),
                                    vcf_parsing.Allele(13)]
        vcf_parsing.fill_genotype_likelihoods('495,162,123,213,129,175,67,0,46,28.1', genotype_to_fill)
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
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        self.assertEqual(0, len(genotype_to_fill.alleles))
        vcf_parsing.fill_genotype_likelihoods('495,162,123', genotype_to_fill)
        self.assertEqual(2, len(genotype_to_fill.alleles))  # should have added two
        self.assertEqual(len(expected_values), len(genotype_to_fill.genotype_likelihoods))
        for index in range(0, len(expected_values)):
            curr_expected_values = expected_values[index]
            real_values = genotype_to_fill.genotype_likelihoods[index]
            self.assertEqual(curr_expected_values[0], real_values.allele1_number)
            self.assertEqual(curr_expected_values[1], real_values.allele2_number)
            self.assertEqual(curr_expected_values[2], real_values.likelihood_neg_exponent)

    def test_fill_genotype_likelihoods_error_too_few_alleles(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [vcf_parsing.Allele(10), vcf_parsing.Allele(11), vcf_parsing.Allele(12)]
        with self.assertRaises(ValueError):
            vcf_parsing.fill_genotype_likelihoods('495,162,123,213,129,175,67,0,46,28.1', genotype_to_fill)

    def test_fill_genotype_likelihoods_error_too_few_likelihoods_1(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [vcf_parsing.Allele(10), vcf_parsing.Allele(11), vcf_parsing.Allele(12),
                                    vcf_parsing.Allele(13)]
        with self.assertRaises(ValueError):
            vcf_parsing.fill_genotype_likelihoods('495,162,123,213,129,175', genotype_to_fill)

    def test_fill_genotype_likelihoods_error_too_few_likelihoods_2(self):
        genotype_to_fill = vcf_parsing.VCFGenotypeInfo('')
        genotype_to_fill.alleles = [vcf_parsing.Allele(10), vcf_parsing.Allele(11), vcf_parsing.Allele(12),
                                    vcf_parsing.Allele(13)]
        with self.assertRaises(ValueError):
            vcf_parsing.fill_genotype_likelihoods('495,162,123,213,129,175,67,0,46', genotype_to_fill)

    # endregion


class TestVCFGenotypeInfo(unittest.TestCase):
    # No tests of __init__ as it is just setting empty values
    # No explicit tests of getter properties as they're just returning an internal variable

    # region is_null_call tests
    def test_is_null_call_true(self):
        self.fail("fill in")

    def test_is_null_call_false(self):
        self.fail("fill in")

    # endregion

    # region genotype_confidence setter tests
    def test_genotype_confidence_setter(self):
        dummy_vcfgenotypeinfo = vcf_parsing.VCFGenotypeInfo('')
        dummy_vcfgenotypeinfo.genotype_confidence = "89"
        self.assertEqual(89.0, dummy_vcfgenotypeinfo.genotype_confidence)
        dummy_vcfgenotypeinfo.genotype_confidence = 89
        self.assertEqual(89.0, dummy_vcfgenotypeinfo.genotype_confidence)
        dummy_vcfgenotypeinfo.genotype_confidence = "-89.10"
        self.assertEqual(-89.1, dummy_vcfgenotypeinfo.genotype_confidence)

    def test_genotype_confidence_setter_error(self):
        dummy_vcfgenotypeinfo = vcf_parsing.VCFGenotypeInfo('')
        with self.assertRaises(ValueError):
            dummy_vcfgenotypeinfo.genotype_confidence = "blue"

    # endregion

    # region filter_passing_reads_count setter tests
    def test_filter_passing_reads_count_setter(self):
        dummy_vcfgenotypeinfo = vcf_parsing.VCFGenotypeInfo('')
        dummy_vcfgenotypeinfo.filter_passing_reads_count = "42"
        self.assertEqual(42, dummy_vcfgenotypeinfo.filter_passing_reads_count)
        dummy_vcfgenotypeinfo.filter_passing_reads_count = 0
        self.assertEqual(0, dummy_vcfgenotypeinfo.filter_passing_reads_count)
        dummy_vcfgenotypeinfo.filter_passing_reads_count = "."
        self.assertEqual(utilities.database.NULL, dummy_vcfgenotypeinfo.filter_passing_reads_count)

    def test_filter_passing_reads_count_setter_error(self):
        dummy_vcfgenotypeinfo = vcf_parsing.VCFGenotypeInfo('')
        with self.assertRaises(ValueError):
            dummy_vcfgenotypeinfo.filter_passing_reads_count = -1
        with self.assertRaises(ValueError):
            dummy_vcfgenotypeinfo.filter_passing_reads_count = 48.5

    # endregion


class TestAllele(unittest.TestCase):
    # No tests of __init__ as it is just setting values
    # No explicit tests of getter properties as they're just returning an internal variable

    def test_read_counts_setter(self):
        dummy_allele = vcf_parsing.Allele(0)
        self.assertEqual(0,dummy_allele.read_counts)
        dummy_allele.read_counts = 94
        self.assertEqual(94, dummy_allele.read_counts)

    def test_read_counts_setter_error(self):
        dummy_allele = vcf_parsing.Allele(0)
        with self.assertRaises(ValueError):
            dummy_allele.read_counts = -94
        with self.assertRaises(ValueError):
            dummy_allele.read_counts = 48.5


class TestGenotypeLikelihood(unittest.TestCase):
    # No tests of init as just calls tested setters
    # No explicit tests of getter properties as they're just returning an internal variable

    # region _validate_allele_relationship tests
    def test__validate_allele_relationship_pass(self):
        vcf_parsing.GenotypeLikelihood._validate_allele_relationship(0,2)
        # if we got this far, the test passed

    def test__validate_allele_relationship_fail(self):
        with self.assertRaises(ValueError):
            vcf_parsing.GenotypeLikelihood._validate_allele_relationship(2,0)

    # endregion

    # region allele1_number setter tests
    def test_allele1_number_setter(self):
        temp_likelihood = vcf_parsing.GenotypeLikelihood(0, 5, 0)
        temp_likelihood.allele1_number = 4
        self.assertEqual(4, temp_likelihood.allele1_number)

    def test_allele1_number_setter_error(self):
        with self.assertRaises(ValueError):
            vcf_parsing.GenotypeLikelihood(0, 0, 0).allele1_number = -1
        with self.assertRaises(ValueError):
            vcf_parsing.GenotypeLikelihood(0, 5, 0).allele1_number = 4.5
        with self.assertRaises(ValueError):
            vcf_parsing.GenotypeLikelihood(0, 0, 0).allele1_number = "blue"

        temp_likelihood = vcf_parsing.GenotypeLikelihood(0, 0, 0)
        temp_likelihood.allele2_number = 1
        with self.assertRaises(ValueError):
            temp_likelihood.allele1_number = 4  # can't be greater than allele 2

    # endregion

    # region allele2_number setter tests
    def test_allele2_number_setter(self):
        temp_likelihood = vcf_parsing.GenotypeLikelihood(0, 0, 0)
        temp_likelihood.allele2_number = 4
        self.assertEqual(4, temp_likelihood.allele2_number)

    def test_allele2_number_setter_error(self):
        with self.assertRaises(ValueError):
            vcf_parsing.GenotypeLikelihood(0, 0, 0).allele2_number = -1
        with self.assertRaises(ValueError):
            vcf_parsing.GenotypeLikelihood(0, 0, 0).allele2_number = 4.5
        with self.assertRaises(ValueError):
            vcf_parsing.GenotypeLikelihood(0, 0, 0).allele2_number = "blue"

        temp_likelihood = vcf_parsing.GenotypeLikelihood(4, 4, 0)
        with self.assertRaises(ValueError):
            temp_likelihood.allele2_number = 2  # can't be smaller than allele 1

    # endregion

    # region likelihood_neg_exponent setter tests
    def test_likelihood_neg_exponent_setter(self):
        temp_likelihood = vcf_parsing.GenotypeLikelihood(0, 0, 0)
        temp_likelihood.likelihood_neg_exponent = "89"
        self.assertEqual(89.0, temp_likelihood.likelihood_neg_exponent)
        temp_likelihood.likelihood_neg_exponent = 89
        self.assertEqual(89.0, temp_likelihood.likelihood_neg_exponent)
        temp_likelihood.likelihood_neg_exponent = "-89.10"
        self.assertEqual(-89.1, temp_likelihood.likelihood_neg_exponent)

    def test_genotype_confidence_setter_error(self):
        temp_likelihood = vcf_parsing.GenotypeLikelihood(0, 0, 0)
        with self.assertRaises(ValueError):
            temp_likelihood.likelihood_neg_exponent = "blue"

    # endregion


class TestVCFGenotypeString(unittest.TestCase):
    def test_parse_GT_GQ_PL(self):
        format_string = 'GT:GQ:PL'
        info_string = '1/1:99:1187.2,101,0'
        parser = vcf_parsing.VCFGenotypeStrings()
        genotype_to_fill = parser.parse(format_string, info_string)
        self.assertEqual('1/1', genotype_to_fill.genotype)
        self.assertEqual("NULL", genotype_to_fill.filter_passing_reads_count)
        self.assertEqual(99, genotype_to_fill.genotype_confidence)
        self.assertEqual(2, len(genotype_to_fill.alleles))
        self.assertEqual(3, len(genotype_to_fill.genotype_likelihoods))
        self.assertEqual(1187.2, genotype_to_fill.genotype_likelihoods[0].likelihood_neg_exponent)
        self.assertEqual(101, genotype_to_fill.genotype_likelihoods[1].likelihood_neg_exponent)
        self.assertEqual(0, genotype_to_fill.genotype_likelihoods[2].likelihood_neg_exponent)

    def test_parse_GT_AD_GQ_PL(self):
        self.fail("fill in")

    def test_parse_GT_AD_DP_GQ_PL(self):
        format_string = 'GT:AD:DP:GQ:PL'
        info_string = '1/1:0,34:34:99:1187.2,101,0'
        parser = vcf_parsing.VCFGenotypeStrings()
        genotype_to_fill = parser.parse(format_string, info_string)

        self.assertEqual('1/1', genotype_to_fill.genotype)
        self.assertEqual(34, genotype_to_fill.filter_passing_reads_count)
        self.assertEqual(99, genotype_to_fill.genotype_confidence)
        self.assertEqual(2, len(genotype_to_fill.alleles))
        self.assertEqual(3, len(genotype_to_fill.genotype_likelihoods))
        self.assertEqual(0, genotype_to_fill.alleles[0].read_counts)
        self.assertEqual(1187.2, genotype_to_fill.genotype_likelihoods[0].likelihood_neg_exponent)
        self.assertEqual(34, genotype_to_fill.alleles[1].read_counts)
        self.assertEqual(101, genotype_to_fill.genotype_likelihoods[1].likelihood_neg_exponent)
        self.assertEqual(0, genotype_to_fill.genotype_likelihoods[2].likelihood_neg_exponent)

    def test_parse_GT_AD_DP_GQ_PGT_PID_PL(self):
        self.fail("fill in")

    def test_parse_error(self):
        self.fail("update test")

        format_string = 'GT:AD:GQ:PL'
        info_string = '1/1:0,34:99:1187.2,101,0'
        parser = vcf_parsing.VCFGenotypeStrings()
        with self.assertRaises(ValueError):
            parser.parse(format_string, info_string)

