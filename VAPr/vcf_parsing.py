# standard libraries
import logging

import VAPr.validation

__author__ = 'Birmingham'


# TODO: rewrite with lambdas or partials
def ignore_pid(info_value, genotype_info_to_fill):
    return ignore_field(info_value, genotype_info_to_fill, 'PID')


# TODO: rewrite with lambdas or partials
def ignore_pgt(info_value, genotype_info_to_fill):
    return ignore_field(info_value, genotype_info_to_fill, 'PGT')


def ignore_field(info_value, genotype_info_to_fill, subkey):
    logging.warning("Ignored subkey '{0}' with value '{1}'".format(subkey, info_value))
    return genotype_info_to_fill


def fill_genotype(info_value, genotype_info_to_fill):
    alleles = info_value.split('/')
    if len(alleles) != 2:
        raise ValueError("Did not detect exactly two alleles in genotype '{0}'".format(info_value))
    genotype_info_to_fill.genotype = info_value
    return genotype_info_to_fill


def fill_unfiltered_reads_counts(info_value, genotype_info_to_fill):
    delimiter = ','
    counts = info_value.split(delimiter)
    if len(counts) < 2:
        raise ValueError("Found fewer than 2 alleles with unfiltered read counts in '{0}'".format(info_value))
    for curr_count in counts:
        new_allele = Allele(curr_count)
        genotype_info_to_fill.alleles.append(new_allele)
    return genotype_info_to_fill


def fill_filtered_reads_count(info_value, genotype_info_to_fill):
    genotype_info_to_fill.filter_passing_reads_count = info_value
    return genotype_info_to_fill


def fill_genotype_confidence(info_value, genotype_info_to_fill):
    genotype_info_to_fill.genotype_confidence = info_value
    return genotype_info_to_fill


def fill_genotype_likelihoods(info_value, genotype_info_to_fill):
    generate_alleles = False
    delimiter = ','
    likelihoods = info_value.split(delimiter)

    num_expected_alleles = len(genotype_info_to_fill.alleles)
    if num_expected_alleles == 0:
        generate_alleles = True
        genotype_info_to_fill.alleles.append(Allele(None))

    allele_number = 0
    likelihood_number = 0
    for index in range(len(likelihoods)):
        if likelihood_number > allele_number:
            if generate_alleles:
                genotype_info_to_fill.alleles.append(Allele(None))
            allele_number += 1
            likelihood_number = 0
            if allele_number >= len(genotype_info_to_fill.alleles):
                raise ValueError("Found {0} likelihoods but only {1} alleles".format(len(likelihoods),
                                                                                     num_expected_alleles))

        new_likelihood = GenotypeLikelihood(likelihood_number, allele_number, likelihoods[index])
        genotype_info_to_fill.genotype_likelihoods.append(new_likelihood)
        likelihood_number += 1

    if allele_number < (num_expected_alleles-1) or likelihood_number < num_expected_alleles:

        raise ValueError("Found {0} alleles but only {1} likelihoods".format(num_expected_alleles, len(likelihoods)))

    return genotype_info_to_fill


class VCFGenotypeInfo:

    def __init__(self, raw_string):
        self.db_id = None
        self.genotype = None
        self._genotype_confidence = None
        self._filter_passing_reads_count = 'NULL'
        self.raw_string = raw_string
        self.alleles = []  # 0 for ref, 1 for first alt, etc
        self.genotype_likelihoods = []

    @property
    def is_null_call(self):
        result = True if self.raw_string.startswith('./.:') else False
        return result

    @property
    def genotype_confidence(self):
        return self._genotype_confidence

    @genotype_confidence.setter
    def genotype_confidence(self, value):
        # TODO: Determine if genotype confidence value is limited to being a positive or non-negative number
        self._genotype_confidence = VAPr.validation.convert_to_nullable(value, float)

    @property
    def filter_passing_reads_count(self):
        return self._filter_passing_reads_count

    @filter_passing_reads_count.setter
    def filter_passing_reads_count(self, value):
        self._filter_passing_reads_count = VAPr.validation.convert_to_nonneg_int(value, nullable=True)


class Allele:
    def __init__(self, read_counts=None):
        self.db_id = None
        self._read_counts = None
        if read_counts is not None:
            self.read_counts = read_counts
        else:
            self._read_counts = 'NULL'

    @property
    def read_counts(self):
        return self._read_counts

    @read_counts.setter
    def read_counts(self, value):
        self._read_counts = VAPr.validation.convert_to_nonneg_int(value, nullable=True)


class GenotypeLikelihood:
    @staticmethod
    def _validate_allele_relationship(allele1_number, allele2_number):
        if allele1_number > allele2_number:
            raise ValueError("VCF-format genotypes must have allele 2 number ({0}) "
                             "greater than or equal to allele 1 number ({1})".format(allele2_number, allele1_number))

    def __init__(self, allele1_number, allele2_number, likelihood_neg_exponent):

        self.db_id = None
        self._allele1_number = None
        self._allele2_number = None
        self._likelihood_neg_exponent = None
        self.allele1_number = allele1_number
        self.allele2_number = allele2_number
        self.likelihood_neg_exponent = likelihood_neg_exponent

    @property
    def allele1_number(self):
        return self._allele1_number

    @allele1_number.setter
    def allele1_number(self, value):
        int_value = VAPr.validation.convert_to_nonneg_int(value)
        if self.allele2_number is not None:
            self._validate_allele_relationship(int_value, self.allele2_number)
        self._allele1_number = int_value

    @property
    def allele2_number(self):
        return self._allele2_number

    @allele2_number.setter
    def allele2_number(self, value):
        int_value = VAPr.validation.convert_to_nonneg_int(value)
        if self.allele1_number is not None:
            self._validate_allele_relationship(self.allele1_number, int_value)
        self._allele2_number = int_value

    @property
    def likelihood_neg_exponent(self):
        return self._likelihood_neg_exponent

    @likelihood_neg_exponent.setter
    def likelihood_neg_exponent(self, value):
        self._likelihood_neg_exponent = VAPr.validation.convert_to_nullable(value, float)


class VCFGenotypeStrings:
    _DELIMITER = ':'
    _PARSER_FUNCS = {'GT': fill_genotype,
                     'AD': fill_unfiltered_reads_counts,
                     'DP': fill_filtered_reads_count,
                     'GQ': fill_genotype_confidence,
                     'PL': fill_genotype_likelihoods,
                     'PID': ignore_pid,
                     'PGT': ignore_pgt}

    @classmethod
    def parse(cls, format_string, info_string):
        result = VCFGenotypeInfo(info_string)

        if format_string not in ('GT:GQ:PL', 'GT:AD:GQ:PL', 'GT:AD:DP:GQ:PL', 'GT:AD:DP:GQ:PGT:PID:PL'):
            raise ValueError("Unrecognized format string: {0}".format(format_string))

        if not result.is_null_call:
            info_subkeys = format_string.split(cls._DELIMITER)
            info_values = info_string.split(cls._DELIMITER)
            for index in range(0, len(info_subkeys)):
                curr_key = info_subkeys[index]
                curr_value = info_values[index]
                parse_func = cls._PARSER_FUNCS[curr_key]
                result = parse_func(curr_value, result)

        return result
