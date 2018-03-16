# coding: utf-8
# standard libraries
import logging
import re
import warnings

# project libraries
import VAPr.validation


def _warn_of_unparseable_format_field(field_desc, field_tag, field_value, failure_desc):
    """Raise a warning for the input field's issue using a standard wording and format.

    Args:
        field_desc (str): A brief human-readable phrase identifying the info indicated by the tag.
        field_tag (str): The two- or three-character key for this field in the VCF format line (e.g., AD, PL, etc).
        field_value (object): The value associated with the field in the VCF format value line; may have been parsed
            to a non-string data type.
        failure_desc (str): A brief human-readable phrase describing the problem identified with the tag_value.
    """
    warn_msg = "The {0} tag value {1} {2} so {3} information could not be captured for the current variant.".format(
        field_tag, field_value, failure_desc, field_desc)
    warnings.warn(warn_msg)


def _capture_unprocessed_field(field_tag, field_value, genotype_info_to_fill):
    """Attempt basic delimiter splitting and/or numeric casting for input field, and store results in VCFGenotypeInfo.

    Args:
        field_tag (str): The two- or three-character key for this field in the VCF format line (e.g., AD, PL, etc).
        field_value (str): The string value associated with the field in the VCF format value line.
        genotype_info_to_fill (VCFGenotypeInfo): A partially filled VCFGenotypeInfo object

    Returns:
        VCFGenotypeInfo: The input VCFGenotypeInfo with additional entries added to unprocessed_info dictionary.

    """
    try:
        genotype_info_to_fill.unprocessed_info[field_tag] = float(field_value)
    except ValueError:  # if the value can't be converted to a float
        split_list = field_value.split(",")
        cast_split_list = []

        for value in split_list:
            try:
                cast_split_list.append(float(value))
            except ValueError:  # if the value can't be converted to a float
                cast_split_list.append(value)

        if len(cast_split_list) > 1:
            genotype_info_to_fill.unprocessed_info[field_tag] = cast_split_list
        else:
            genotype_info_to_fill.unprocessed_info[field_tag] = cast_split_list[0]

    return genotype_info_to_fill


def _fill_genotype_class(alleles, genotype_info_to_fill):
    """Id genotype class (homozygous/heterozygous) and subclass (reference, alt, compound) and store in VCFGenotypeInfo.

    Args:
        alleles (List[str]): Results of splitting the value for the GT (genotype) format tag.
        genotype_info_to_fill (VCFGenotypeInfo): A partially filled VCFGenotypeInfo object.

    Returns:
        VCFGenotypeInfo: The input VCFGenotypeInfo with genotype_subclass_by_class value set.

    """
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


def _fill_genotype(field_value, genotype_info_to_fill):
    """Parse the genotype of this sample at this site and store in the VCFGenotypeInfo.

    Args:
        field_value (str): The value associated with the GT tag in the format value string.
        genotype_info_to_fill (VCFGenotypeInfo): A partially filled VCFGenotypeInfo object.

    Returns:
        VCFGenotypeInfo: the input VCFGenotypeInfo with additional fields filled in.

    From https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it :
    " GT : The genotype of this sample at this site.
        For a diploid organism, the GT field indicates the two alleles carried by the sample, encoded by a 0 for the REF
         allele, 1 for the first ALT allele, 2 for the second ALT allele, etc. When there's a single ALT allele (by far
         the more common case), GT will be either:

        0/0 - the sample is homozygous reference
        0/1 - the sample is heterozygous, carrying 1 copy of each of the REF and ALT alleles
        1/1 - the sample is homozygous alternate

      In the three sites shown in the example above, NA12878 is observed with the allele combinations T/G, G/G, and C/T
      respectively.

      For non-diploids, the same pattern applies; in the haploid case there will be just a single value in GT; for
      polyploids there will be more, e.g. 4 values for a tetraploid organism."
    """

    alleles = field_value.split('/')
    if len(alleles) == 1:
        alleles = field_value.split('|')

    if len(alleles) != 2:
        _warn_of_unparseable_format_field("genotype", VCFGenotypeParser.GENOTYPE_TAG, field_value,
                                          "does not split into exactly two values")
        return genotype_info_to_fill
    genotype_info_to_fill.genotype = field_value
    _fill_genotype_class(alleles, genotype_info_to_fill)
    return genotype_info_to_fill


def _fill_unfiltered_reads_counts(field_value, genotype_info_to_fill):
    """Parse the unfiltered reads counts for this sample at this site and store in the VCFGenotypeInfo.

    Args:
        field_value (str): The value associated with the AD tag in the format value string.
        genotype_info_to_fill (VCFGenotypeInfo): A partially filled VCFGenotypeInfo object.

    Returns:
        VCFGenotypeInfo: The input VCFGenotypeInfo with new entries added to the alleles list.

    From https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it  :
    " AD ... : Allele depth ....
      AD is the unfiltered allele depth, i.e. the number of reads that support each of the reported alleles. All reads
      at the position (including reads that did not pass the variant caller’s filters) are included in this number,
      except reads that were considered uninformative. Reads are considered uninformative when they do not provide
      enough statistical evidence to support one allele over another."
    """
    delimiter = ','
    counts = field_value.split(delimiter)
    if len(counts) < 2:
        _warn_of_unparseable_format_field("unfiltered allele depth", VCFGenotypeParser.UNFILTERED_ALLELE_DEPTH_TAG,
                                          field_value, "does not split into at least two values")
    else:
        for curr_count in counts:
            new_allele = Allele(curr_count)
            genotype_info_to_fill.alleles.append(new_allele)

    return genotype_info_to_fill


def _fill_filtered_reads_count(field_value, genotype_info_to_fill):
    """Parse the filtered depth of coverage of this sample at this site and store in the VCFGenotypeInfo.

    Args:
        field_value (str): The value associated with the DP tag in the format value string.
        genotype_info_to_fill (VCFGenotypeInfo): A partially filled VCFGenotypeInfo object.

    Returns:
        VCFGenotypeInfo: The input VCFGenotypeInfo with filter_passing_reads_count filled in.

    From https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it :
    " DP : ... depth of coverage
      DP is the filtered depth, at the sample level. This gives you the number of filtered reads that support each of
      the reported alleles. You can check the variant caller’s documentation to see which filters are applied by
      default. Only reads that passed the variant caller’s filters are included in this number. However, unlike the AD
      calculation, uninformative reads are included in DP."
    """
    genotype_info_to_fill.filter_passing_reads_count = field_value
    return genotype_info_to_fill


def _fill_genotype_confidence(field_value, genotype_info_to_fill):
    """Parse the genotype quality (confidence) of this sample at this site and store in the VCFGenotypeInfo.

    Args:
        field_value (str): The value associated with the GQ tag in the format value string.
        genotype_info_to_fill (VCFGenotypeInfo): A partially filled VCFGenotypeInfo object.

    Returns:
        VCFGenotypeInfo: The input VCFGenotypeInfo with genotype_confidence filled in.

    From https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it :
    " GQ : Quality of the assigned genotype.
      The Genotype Quality represents the Phred-scaled confidence that the genotype assignment (GT) is correct, derived
      from the genotype PLs. Specifically, the GQ is the difference between the PL of the second most likely genotype,
      and the PL of the most likely genotype. As noted above, the values of the PLs are normalized so that the most
      likely PL is always 0, so the GQ ends up being equal to the second smallest PL, unless that PL is greater than 99.
      In GATK, the value of GQ is capped at 99 because larger values are not more informative, but they take more space
      in the file. So if the second most likely PL is greater than 99, we still assign a GQ of 99.

      Basically the GQ gives you the difference between the likelihoods of the two most likely genotypes. If it is low,
      you can tell there is not much confidence in the genotype, i.e. there was not enough evidence to confidently
      choose one genotype over another. See the FAQ article on the Phred scale to get a sense of what would be
      considered low.

      Not to be confused with the site-level annotation QUAL; see this FAQ article for an explanation of the differences
      in what they mean and how they should be used."
    """
    genotype_info_to_fill.genotype_confidence = field_value
    return genotype_info_to_fill


def _fill_genotype_likelihoods(field_value, genotype_info_to_fill):
    """Parse the "normalized" Phred-scaled likelihoods of possible genotypes of this sample at this site and store.

    Args:
        field_value (str): The value associated with the PL tag in the format value string.
        genotype_info_to_fill (VCFGenotypeInfo): A partially filled VCFGenotypeInfo object.

    Returns:
        VCFGenotypeInfo: The input VCFGenotypeInfo with additional fields filled in.

    Note that this function will MAKE the number of Allele objects implied by the likelihood string if no alleles have
    been filled into the genotype_info_to_fill by the time this function is called. The reason for this is that there
    ARE valid VCF format strings (e.g., 'GT:GQ:PL') that have the PL tag (likelihood) but no AD tag in them, and since
    Allele objects are usually created in the processing of the AD tag, some back-up approach was needed to infer
    Alleles in this situation.  Of course, the Alleles created in this situation will all have None as their
    read_counts, since the read_counts value comes from the AD tag.

    From https://gatkforums.broadinstitute.org/gatk/discussion/1268/what-is-a-vcf-and-how-should-i-interpret-it:
    " PL : "Normalized" Phred-scaled likelihoods of the possible genotypes.
      For the typical case of a monomorphic site (where there is only one ALT allele) in a diploid organism, the PL
      field will contain three numbers, corresponding to the three possible genotypes (0/0, 0/1, and 1/1). The PL values
      are "normalized" so that the PL of the most likely genotype (assigned in the GT field) is 0 in the Phred scale. We
      use "normalized" in quotes because these are not probabilities. We set the most likely genotype PL to 0 for easy
      reading purpose.The other values are scaled relative to this most likely genotype.

      Keep in mind, if you're not familiar with the statistical lingo, that when we say PL is the "Phred-scaled
      likelihood of the genotype", we mean it is "How much less likely that genotype is compared to the best one".
    "
    """
    generate_alleles = False
    delimiter = ','
    likelihoods = field_value.split(delimiter)

    num_expected_alleles = len(genotype_info_to_fill.alleles)
    if num_expected_alleles == 0:
        generate_alleles = True
        genotype_info_to_fill.alleles.append(Allele(None))

    allele_number = 0
    likelihood_number = 0
    for index in range(len(likelihoods)):
        if likelihood_number > allele_number:
            allele_number += 1
            likelihood_number = 0

            if generate_alleles:
                genotype_info_to_fill.alleles.append(Allele(None))

            # NB: The check below looks at the *current* number of alleles in genotype_info_to_fill, which may be
            # different than the value found above and stored in num_expected_alleles due to the alleles.append
            # statements at various points above.
            if allele_number >= len(genotype_info_to_fill.alleles):
                _warn_of_unparseable_format_field("'normalized' Phred-scaled likelihoods of possible genotypes",
                                                  VCFGenotypeParser.NORMALIZED_SCALED_LIKELIHOODS_TAG,
                                                  field_value,
                                                  "appears to contain information for more alleles than expected")

                # in case of warning, clear any already-parsed genotype likelihoods (they probably can't be trusted)
                # and return out of function
                genotype_info_to_fill.genotype_likelihoods = []
                return genotype_info_to_fill
            # end warn
        # end if likelihood_number > allele_number

        new_likelihood = GenotypeLikelihood(likelihood_number, allele_number, likelihoods[index])
        genotype_info_to_fill.genotype_likelihoods.append(new_likelihood)
        likelihood_number += 1

    if allele_number < (num_expected_alleles-1) or likelihood_number < num_expected_alleles:
        _warn_of_unparseable_format_field("'normalized' Phred-scaled likelihoods of possible genotypes",
                                          VCFGenotypeParser.NORMALIZED_SCALED_LIKELIHOODS_TAG,
                                          field_value, "appears to contain information for fewer alleles than expected")

        # in case of warning, clear any already-parsed genotype likelihoods (they probably can't be trusted)
        genotype_info_to_fill.genotype_likelihoods = []

    return genotype_info_to_fill


class VCFGenotypeInfo(object):
    """Store parsed info from VCF genotype fields for a single sample.

    Attributes:
        _raw_string (str): The genotype fields values string from a VCF file (e.g., '0/1:173,141:282:99:255,0,255').
        genotype (Optional[`str`]): The type of each of the sample's two alleles, such as 0/0, 0/1, etc.
        alleles (List[Allele]): One Allele object for each allele detected for this variant
            (this can be across samples, so there can be more than 2 alleles).
        genotype_likelihoods (List[GenotypeLikelihood]): The GenotypeLikelihood object for each allele.
        unprocessed_info (Dict[str, Any]): Dictionary of field tag and value(s) for any fields not stored in
            dedicated attributes of VCFGenotypeInfo.  Values are parsed to lists and/or floats if possible.
        genotype_subclass_by_class (Dict[str, str]): Genotype subclass (reference, alt, compound) keyed by genotype
            class (homozygous/heterozygous).
    """

    def __init__(self, raw_string):
        """Create VCFGenotypeInfo object.

        Args:
            raw_string (str): The genotype fields values string from a VCF file (e.g., '0/1:173,141:282:99:255,0,255').
        """
        self._raw_string = raw_string
        self._genotype_confidence = None
        self._filter_passing_reads_count = None

        # TODO: someday: Probably these should become properties so they are protected from user resetting them, etc.
        self.genotype = None
        self.alleles = []  # 0 is ref, 1 is first alt, etc
        self.genotype_likelihoods = []
        self.unprocessed_info = {}
        self.genotype_subclass_by_class = None

    @property
    def genotype_confidence(self):
        """str: Genotype quality (confidence) of this sample at this site, from the GQ field."""
        return self._genotype_confidence

    @genotype_confidence.setter
    def genotype_confidence(self, value):
        # TODO: someday: Determine if genotype confidence value is limited to being a positive or non-negative number
        self._genotype_confidence = VAPr.validation.convert_to_nullable(value, float)

    @property
    def filter_passing_reads_count(self):
        """int or None: Filtered depth of coverage of this sample at this site from the DP field."""
        return self._filter_passing_reads_count

    @filter_passing_reads_count.setter
    def filter_passing_reads_count(self, value):
        self._filter_passing_reads_count = VAPr.validation.convert_to_nonneg_int(value, nullable=True)


# TODO: someday: refactor to remove Allele object and replace with single value
# One could argue pretty convincingly that there is no longer a need for a whole Allele object now that
# (with the removal of database-related properties) it has only one property--surely that could be represented by
# a single value rather than a class!  But this would be a non-trivial refactor, so I'm leaving it for some misty future
# date.
class Allele(object):
    """Store unfiltered read counts, if any, for a particular allele."""

    def __init__(self, unfiltered_read_counts=None):
        """Create Allele object.

        Args:
            unfiltered_read_counts (Optional[str]): Number of unfiltered reads counts for this sample at this site,
                from AD field.
        """
        self._unfiltered_read_counts = None
        if unfiltered_read_counts is not None:
            self.unfiltered_read_counts = unfiltered_read_counts

    @property
    def unfiltered_read_counts(self):
        """int or None: Number of unfiltered reads counts for this sample at this site, from AD field."""
        return self._unfiltered_read_counts

    @unfiltered_read_counts.setter
    def unfiltered_read_counts(self, value):
        self._unfiltered_read_counts = VAPr.validation.convert_to_nonneg_int(value, nullable=True)


class GenotypeLikelihood(object):
    """Store parsed info from VCF genotype likelihood field for a single sample."""

    @staticmethod
    def _validate_allele_relationship(allele1_number, allele2_number):
        """Ensure that allele1_number is not greater than allele2_number.

        Args:
            allele1_number (int): The allele identifier (0 for reference, 1 for first alternate, etc) for the left-hand
                allele inferred for this genotype likelihood.
            allele2_number (int): The allele identifier (0 for reference, 1 for first alternate, etc) for the right-hand
                allele inferred for this genotype likelihood.

        Raises:
            ValueError: If `allele1_number` is greater than `allele2_number`.

        Genotype likelihood strings contain one likelihood for each possible allele combination for this variant at
        this site.  The likelihood string does not explicity state the allele combination associated with each
        likelihood; rather, the allele combinations are inferred based on the convention that they are always listed
        in ascending order from lowest allele number to highest allele number. For example, for the likelihood string
        '495,162,123,213,129,175,67,0,46,28.1', the implied allele combinations are  0/0, 0/1, 1/1, 0/2, 1/2, 2/2,
        0/3, 1/3, 2/3, 3/3 .  Note that these are COMBINATIONS, not PERMUTATIONS, so each pair of values occurs only
        once and, again by convention, the representation expected is the one where the allele number on the left
        (arbitrarily referred to here as allele1) has a value less than (or equal to) the allele on the right
        (arbitrarily referred to here as allele2).

        IF, in parsing the likelihood string and inferring the allele combinations, the code somehow ended up with
        an inferred combination in which the left-hand allele number is GREATER than the right-hand allele number,
        that means something has gone tragically wrong!
        """

        if allele1_number > allele2_number:
            raise ValueError("VCF-format genotypes must have allele 2 number ({0}) "
                             "greater than or equal to allele 1 number ({1})".format(allele2_number, allele1_number))

    def __init__(self, allele1_number, allele2_number, likelihood_neg_exponent):
        """Create GenotypeLikelihood object.

        Args:
            allele1_number (int or str): The allele id for the right-hand allele inferred for this genotype likelihood.
            allele2_number (int or str): The allele id for the left-hand allele inferred for this genotype likelihood.
            likelihood_neg_exponent (float or str): The "normalized" Phred-scaled likelihood of the genotype represented
                by allele1 and allele2, as a string.
        """
        self._allele1_number = None
        self._allele2_number = None
        self._likelihood_neg_exponent = None

        self.allele1_number = allele1_number
        self.allele2_number = allele2_number
        self.likelihood_neg_exponent = likelihood_neg_exponent

    @property
    def allele1_number(self):
        """int: The allele identifier for the left-hand allele inferred for this genotype likelihood."""
        return self._allele1_number

    @allele1_number.setter
    def allele1_number(self, value):
        int_value = VAPr.validation.convert_to_nonneg_int(value, nullable=True)

        if self.allele2_number is not None:
            self._validate_allele_relationship(int_value, self.allele2_number)
        self._allele1_number = int_value

    @property
    def allele2_number(self):
        """int: The allele identifier for the right-hand allele inferred for this genotype likelihood. """
        return self._allele2_number

    @allele2_number.setter
    def allele2_number(self, value):
        int_value = VAPr.validation.convert_to_nonneg_int(value, nullable=True)

        if self.allele1_number is not None:
            self._validate_allele_relationship(self.allele1_number, int_value)
        self._allele2_number = int_value

    @property
    def likelihood_neg_exponent(self):
        """float: The "normalized" Phred-scaled likelihood of the genotype represented by allele1 and allele2."""
        return self._likelihood_neg_exponent

    @likelihood_neg_exponent.setter
    def likelihood_neg_exponent(self, value):
        self._likelihood_neg_exponent = VAPr.validation.convert_to_nullable(value, float)


class VCFGenotypeParser(object):
    """Mine format string and genotype fields string to create a filled VCFGenotypeInfo object."""

    @staticmethod
    def is_valid_genotype_fields_string(genotype_fields_string):
        """Return true if input has any real genotype fields content, false if is just periods, zeroes, and delimiters.

        Args:
            genotype_fields_string (str): A VCF-style genotype fields string, such as 1/1:0,2:2:6:89,6,0 or ./.:.:.:.:.

        Returns
            bool: true if input has any real genotype fields content, false if is just periods, zeroes, and delimiters.
        """
        result = False
        # this regex means "one or more characters that is not a comma, period, colon, zero, or forward slash"
        content_char_match = re.search(r"[^,.:0\/]+", genotype_fields_string)
        # NB: necessary to ALSO check first character of string, even if no match to above regex is found, because
        # "0/0" should be a valid genotype (even though all the characters it contains could signal null content in
        # other configurations), and a genotype fields string with nothing but a genotype in it should be legal.
        if content_char_match is not None or not genotype_fields_string.startswith("."):
            result = True
        return result

    GENOTYPE_TAG = "GT"  # str: VCF tag for the genotype of this sample at this site.
    UNFILTERED_ALLELE_DEPTH_TAG = "AD"  # str: VCF tag for the unfiltered allele depth of this sample at this site.
    FILTERED_ALLELE_DEPTH_TAG = "DP"  # str: VCF tag for the filtered depth of coverage of this sample at this site.
    GENOTYPE_QUALITY_TAG = "GQ"  # str: VCF tag for the genotype quality of this sample at this site.
    NORMALIZED_SCALED_LIKELIHOODS_TAG = "PL"  # str: VCF tag for the genotype likelihoods of this sample at this site.

    _DELIMITER = ':'  # str: Delimiter between fields in format and genotype fields strings.

    # Dict(str, Callable[str, VCFGenotypeInfo]): Special parsing functions by the VCF tag whose value they parse.
    _PARSER_FUNCS = {GENOTYPE_TAG: _fill_genotype,  # GT
                     UNFILTERED_ALLELE_DEPTH_TAG: _fill_unfiltered_reads_counts,  # AD
                     FILTERED_ALLELE_DEPTH_TAG: _fill_filtered_reads_count,  # DP
                     GENOTYPE_QUALITY_TAG: _fill_genotype_confidence,  # GQ
                     NORMALIZED_SCALED_LIKELIHOODS_TAG: _fill_genotype_likelihoods}  # PL

    @classmethod
    def parse(cls, format_key_string, format_value_string):
        """Parse the input format string and genotype fields string into a filled VCFGenotypeInfo object.

        Args:
            format_key_string (str): The VCF format string (e.g., 'GT:AD:DP:GQ:PL') for this sample at this site.
            format_value_string (str): The VCF genotype fields values string (e.g., '1/1:0,34:34:99:1187.2,101,0')
                corresponding to the format_key_string for this sample at this site.

        Returns:
            VCFGenotypeInfo or None: A filled VCFGenotypeInfo for this sample at this site unless an error was
                encountered, in which case None is returned.
                encountered, in which case None is returned.

        """
        result = None

        try:
            if cls.is_valid_genotype_fields_string(format_value_string):
                result = VCFGenotypeInfo(format_value_string)
                format_subkeys = format_key_string.split(cls._DELIMITER)
                format_values = format_value_string.split(cls._DELIMITER)

                for index, curr_key in enumerate(format_subkeys):
                    curr_value = format_values[index]
                    if curr_key in cls._PARSER_FUNCS:
                        # if this key has a special parsing function associated with it, use that function
                        parse_func = cls._PARSER_FUNCS[curr_key]
                        result = parse_func(curr_value, result)
                    else:
                        # otherwise, capture this key/value with minimal processing to a catch-all dictionary
                        result = _capture_unprocessed_field(curr_key, curr_value, result)
        except Exception as e:
            warn_msg = "Encountered error '{0}' so genotype fields information could not be captured for the " \
                       "current variant.".format(e)
            warnings.warn(warn_msg)
            result = None  # reset result to None, as contents can't be trusted

        return result
