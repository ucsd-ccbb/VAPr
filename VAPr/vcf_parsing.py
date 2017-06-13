# standard libraries
import logging
import sys
from VAPr import validation

logger = logging.getLogger()
logger.setLevel(logging.INFO)
#try:
#    logger.handlers[0].stream = sys.stdout
#except:
#    pass

GLOBAL_ERROR_COUNT = 0

__author__ = 'Birmingham'

# TODO: Forgive mal-formatted vcf files. Log the error but dont break the annotation.


# TODO: rewrite with lambdas or partials
def ignore_pid(info_value, genotype_info_to_fill):
    return ignore_field(info_value, genotype_info_to_fill, 'PID')


# TODO: rewrite with lambdas or partials
def ignore_pgt(info_value, genotype_info_to_fill):
    return ignore_field(info_value, genotype_info_to_fill, 'PGT')


def ignore_field(info_value, genotype_info_to_fill, subkey):
    return genotype_info_to_fill


def fill_genotype_class(alleles, genotype_info_to_fill):
    genotype_class = "homozygous"
    genotype_subclass = "reference"
    alt_subclass_name = "alt"

    if alleles[0] != alleles[1]:
        genotype_class = "heterozygous"
        alt_subclass_name = "compound"

    if "0" not in alleles:
        genotype_subclass = alt_subclass_name

    result = {genotype_class: genotype_subclass}
    genotype_info_to_fill.genotype_subclass_by_class = result
    return result


def fill_genotype(info_value, genotype_info_to_fill):
    global GLOBAL_ERROR_COUNT
    alleles = info_value.split('/')
    if len(alleles) != 2:
        GLOBAL_ERROR_COUNT += 1
    genotype_info_to_fill.genotype = info_value
    fill_genotype_class(alleles, genotype_info_to_fill)
    return genotype_info_to_fill


def fill_unfiltered_reads_counts(info_value, genotype_info_to_fill):
    global GLOBAL_ERROR_COUNT
    delimiter = ','
    counts = info_value.split(delimiter)
    if len(counts) < 2:
        GLOBAL_ERROR_COUNT += 1
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
    global GLOBAL_ERROR_COUNT
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
                GLOBAL_ERROR_COUNT += 1
        new_likelihood = GenotypeLikelihood(likelihood_number, allele_number, likelihoods[index])
        genotype_info_to_fill.genotype_likelihoods.append(new_likelihood)
        likelihood_number += 1

    if allele_number < (num_expected_alleles-1) or likelihood_number < num_expected_alleles:
        GLOBAL_ERROR_COUNT += 1

    return genotype_info_to_fill


class VCFGenotypeInfo:

    def __init__(self, raw_string):
        self.db_id = None
        self.genotype = None
        self._genotype_confidence = None
        self._filter_passing_reads_count = None  # NULL
        self.raw_string = raw_string
        self.alleles = []  # 0 for ref, 1 for first alt, etc
        self.genotype_likelihoods = []
        self.errors = GLOBAL_ERROR_COUNT

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
        self._genotype_confidence = validation.convert_to_nullable(value, float)

    @property
    def filter_passing_reads_count(self):
        return self._filter_passing_reads_count

    @filter_passing_reads_count.setter
    def filter_passing_reads_count(self, value):
        self._filter_passing_reads_count = validation.convert_to_nonneg_int(value, nullable=True)[0]
        self.errors += validation.convert_to_nonneg_int(value, nullable=True)[1]


class Allele:
    def __init__(self, read_counts=None):
        self.db_id = None
        self._read_counts = None
        if read_counts is not None:
            self.read_counts = read_counts
        else:
            self._read_counts = 'NULL'
        self.errors = GLOBAL_ERROR_COUNT

    @property
    def read_counts(self):
        return self._read_counts

    @read_counts.setter
    def read_counts(self, value):
        self._read_counts = validation.convert_to_nonneg_int(value, nullable=True)[0]
        self.errors += validation.convert_to_nonneg_int(value, nullable=True)[1]


class GenotypeLikelihood:
    @staticmethod
    def _validate_allele_relationship(allele1_number, allele2_number):
        global GLOBAL_ERROR_COUNT
        if allele1_number > allele2_number:
            GLOBAL_ERROR_COUNT += 1

    def __init__(self, allele1_number, allele2_number, likelihood_neg_exponent):

        self.db_id = None
        self._allele1_number = None
        self._allele2_number = None
        self._likelihood_neg_exponent = None
        self.allele1_number = allele1_number
        self.allele2_number = allele2_number
        self.likelihood_neg_exponent = likelihood_neg_exponent
        self.genotype_subclass_by_class = None
        self.errors = GLOBAL_ERROR_COUNT

    @property
    def allele1_number(self):
        return self._allele1_number

    @allele1_number.setter
    def allele1_number(self, value):
        int_value = validation.convert_to_nonneg_int(value, nullable=True)[0]
        self.errors += validation.convert_to_nonneg_int(value, nullable=True)[1]
        if self.allele2_number is not None:
            self._validate_allele_relationship(int_value, self.allele2_number)
        self._allele1_number = int_value

    @property
    def allele2_number(self):
        return self._allele2_number

    @allele2_number.setter
    def allele2_number(self, value):
        int_value = validation.convert_to_nonneg_int(value, nullable=True)[0]
        self.errors += validation.convert_to_nonneg_int(value, nullable=True)[1]
        if self.allele1_number is not None:
            self._validate_allele_relationship(self.allele1_number, int_value)
        self._allele2_number = int_value

    @property
    def likelihood_neg_exponent(self):
        return self._likelihood_neg_exponent

    @likelihood_neg_exponent.setter
    def likelihood_neg_exponent(self, value):
        self._likelihood_neg_exponent = validation.convert_to_nullable(value, float)

    # GT:AD:BQ:DP:FA

# AD BQ DP FA GQ GT PL SS

class VCFGenotypeStrings:
    #TODO: Parser for Mutect files

    _DELIMITER = ':'
    _PARSER_FUNCS = {'GT': fill_genotype,
                     'AD': fill_unfiltered_reads_counts,
                     # 'BQ': fill_base_quality,
                     'DP': fill_filtered_reads_count,
                     # 'FA': fill_allele_fraction,
                     'GQ': fill_genotype_confidence,
                     'PL': fill_genotype_likelihoods,
                     # 'SS': fill_variant_status,
                     'PID': ignore_pid,
                     'PGT': ignore_pgt}

    @classmethod
    def parse(cls, format_string, info_string):

        mutect_formats = ['BQ', 'FA', 'SS']
        possible_string_formats = ['GT:GQ:PL', 'GT:AD:GQ:PL', 'GT:AD:DP:GQ:PL', 'GT:AD:DP:GQ:PGT:PID:PL']
        result = VCFGenotypeInfo(info_string)

        if not result.is_null_call:
            info_subkeys = format_string.split(cls._DELIMITER)
            info_values = info_string.split(cls._DELIMITER)

            for index, value in enumerate(info_subkeys):
                if value in mutect_formats:
                    continue
                if value not in possible_string_formats[-1]:
                    continue
                curr_key = info_subkeys[index]
                curr_value = info_values[index]
                parse_func = cls._PARSER_FUNCS[curr_key]
                result = parse_func(curr_value, result)

        return result

"""
if __name__ == '__main__':

    for i in range(0,10000):
        if i%2:
            dictionary = {'otherinfo': ['GT:AD:BQ:DP:FA', '0:0,3:.:3:1.00']}
        else:
            dictionary = {'otherinfo': ['GT:AD:GQ:PL', '1/1:0,43:99:1039,129,0']}
        parser = VCFGenotypeStrings()

        genotype_to_fill = parser.parse(dictionary['otherinfo'][0], dictionary['otherinfo'][1])
        gen_dic = {'genotype': genotype_to_fill.genotype}

        genotype = genotype_to_fill.genotype
        filter_passing_reads_count = genotype_to_fill.filter_passing_reads_count
        genotype_likelihoods = genotype_to_fill.genotype_likelihoods
        alleles = genotype_to_fill.alleles

        try:
            gen_dic['filter_passing_reads_count'] = [int(genotype_to_fill.filter_passing_reads_count)]
        except ValueError:
            print('Read depth returned null value')
        try:
            gen_dic['genotype_likelihoods'] = [float(genotype_to_fill.genotype_likelihoods[0].likelihood_neg_exponent),
                                               float(genotype_to_fill.genotype_likelihoods[1].likelihood_neg_exponent),
                                               float(genotype_to_fill.genotype_likelihoods[2].likelihood_neg_exponent)]
        except IndexError:
            print('Genotype Likelihood not understood')

        try:
            gen_dic['alleles'] = [int(genotype_to_fill.alleles[0].read_counts),
                                  int(genotype_to_fill.alleles[1].read_counts)]
        except IndexError:
            print('Allele read counts  not understood')

        print(gen_dic)
"""