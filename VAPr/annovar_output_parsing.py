# built-in libraries
import csv
import itertools
import logging
import sys

# third-party libraries
import myvariant

# project libraries
from VAPr.vcf_genotype_fields_parsing import VCFGenotypeParser


class AnnovarTxtParser(object):
    """ Class that processes an Annovar-created tab-delimited text file."""

    CHR_HEADER = 'chr'
    START_HEADER = 'start'
    END_HEADER = 'end'
    REF_HEADER = 'ref'
    ALT_HEADER = 'alt'
    OTHERINFO_HEADER = 'otherinfo'
    THOU_G_2015_ALL_HEADER = '1000g2015aug_all'
    ESP6500_ALL_HEADER = 'esp6500siv2_all'
    NCI60_HEADER = 'nci60'
    CYTOBAND_HEADER = 'cytoband'
    GENOMIC_SUPERDUPS_HEADER = 'genomicsuperdups'
    TFBS_CONS_SITES_HEADER = 'tfbsconssites'
    FUNC_KNOWNGENE_HEADER = 'func_knowngene'
    GENE_KNOWNGENE_HEADER = 'gene_knowngene'
    GENEDETAIL_KNOWNGENE_HEADER = 'genedetail_knowngene'
    EXONICFUNC_KNOWNGENE_HEADER = 'exonicfunc_knowngene'
    SCORE_KEY = "Score"
    RAW_CHR_MT_VAL = "chrM"
    RAW_CHR_MT_SUFFIX_VAL = RAW_CHR_MT_VAL.replace(CHR_HEADER, "")
    STANDARDIZED_CHR_MT_VAL = "chrMT"
    STANDARDIZED_CHR_MT_SUFFIX_VAL = STANDARDIZED_CHR_MT_VAL.replace(CHR_HEADER, "")

    _ANNOVAR_OUTPUT_COLS = [CHR_HEADER,
                           START_HEADER,
                           END_HEADER,
                           REF_HEADER,
                           ALT_HEADER,
                           FUNC_KNOWNGENE_HEADER,
                           GENE_KNOWNGENE_HEADER,
                           GENEDETAIL_KNOWNGENE_HEADER,
                           EXONICFUNC_KNOWNGENE_HEADER,
                           # CYTOBAND_HEADER,
                           # GENOMIC_SUPERDUPS_HEADER,
                           THOU_G_2015_ALL_HEADER] #,
                           # ESP6500_ALL_HEADER,
                           # 'cosmic70',
                           # NCI60_HEADER,
                           #OTHERINFO_HEADER]

    @staticmethod
    def _normalize_header(raw_headers_list):
        """Lower-case all header strings and replace periods with underscores.

        Args:
            raw_headers_list (List[str]): A list of the headers from the Annovar txt output file.

        Returns:
            List[str]: A list of the normalized headers, in the same order as the input list.
        """

        normalized_headers_list = []

        for curr_header in raw_headers_list:
            normalized_headers_list.append(curr_header.lower().replace('.', '_'))

        return normalized_headers_list

    @classmethod
    def read_chunk_of_annotations_to_dicts_list(cls, annovar_txt_file_like_obj, sample_names_list, chunk_index,
                                                chunk_size):
        annotations_dict_per_variant_list = []
        hgvsid_list = []

        # read in the header line and normalize the header names
        reader = csv.reader(annovar_txt_file_like_obj, delimiter='\t')
        normed_headers_list = cls._normalize_header(next(reader))

        # for each row in this chunk--which is to say, each variant
        for curr_line_fields_list in itertools.islice(reader, (chunk_index * chunk_size),
                                                      ((chunk_index + 1) * chunk_size)):
            hgvs_id, annotations_dict_for_curr_variant = cls._parse_single_variant_record(
                normed_headers_list, curr_line_fields_list, sample_names_list)
            hgvsid_list.append(hgvs_id)
            annotations_dict_per_variant_list.append(annotations_dict_for_curr_variant)

        return hgvsid_list, annotations_dict_per_variant_list

    @classmethod
    def _parse_single_variant_record(cls, normed_headers_list, curr_line_fields_list, sample_names_list):
        # This code assumes that the VCF-produced format string and the genotype fields string(s) for the sample(s)
        # will be the last fields on every line and that they will NOT all have their own headers--rather, it
        # assumes the last header will indicate that the rest of the fields are "other info".  Here is a simplified
        # example:
        # chr	start	end	ref	alt	func_knowngene	    otherinfo
        # chrM	146     146	T	C	upstream;downstream 1	    61.74	AC=2;AF=1.00;AN=2;DP=2;FS=0.000	GT:AD:DP:GQ:PL	1/1:0,22:22:66:794,66,0	./.:0,0	1/1:0,40:40:99:1494,119,0
        # Note that the content at and after the position of the otherinfo header may list additional information
        # before the format string and genotype fields info (which are required to be at the end of the line);
        # any such extra info is ignored.

        # make a dictionary that pairs every named header *except* the last one with its content in this line
        last_field_index = len(normed_headers_list) - 1
        raw_fields_dict = dict(zip(normed_headers_list[0:last_field_index], curr_line_fields_list[0:last_field_index]))

        # TODO: someday: perhaps stop limiting this to only the fields in ANNOVAR_OUTPUT_COLS instead of all
        # For only a limited subset of columns, look those columns up in raw_fields_dict; if they hold real content,
        # do any clean-up necessary to their values and write them into a new dict
        cleaned_fields_dict = {}
        for curr_header in cls._ANNOVAR_OUTPUT_COLS:
            curr_value = raw_fields_dict[curr_header]
            if curr_value != ".":
                curr_value = cls._rewrite_value_if_special_header(curr_header, curr_value)
                cleaned_fields_dict[curr_header] = curr_value

        # generate the hgvs id for this variant
        hgvs_id = myvariant.format_hgvs(cleaned_fields_dict[cls.CHR_HEADER],
                                        cleaned_fields_dict[cls.START_HEADER],
                                        cleaned_fields_dict[cls.REF_HEADER],
                                        cleaned_fields_dict[cls.ALT_HEADER])

        # now grab the number-of-samples-plus-one-th field from the *end* of the line--this holds the format
        # string--and also grab a list of the number-of-samples fields from the *end* of the line--these are
        # the genotype fields strings for each sample.
        num_samples_plus_one = len(sample_names_list) + 1
        format_string = curr_line_fields_list[-num_samples_plus_one]
        genotype_field_strings_per_sample = curr_line_fields_list[-len(sample_names_list)::]
        genotype_field_strings_by_sample_name = dict(zip(sample_names_list, genotype_field_strings_per_sample))

        # turn the dictionary of annovar fields into a dictionary of annotations for the variant, including
        # nested structures containing sample-specific genotype-related info
        annotations_dict_for_curr_variant = AnnovarAnnotatedVariant.make_per_variant_annotation_dict(
            cleaned_fields_dict, hgvs_id, format_string, genotype_field_strings_by_sample_name)

        return hgvs_id, annotations_dict_for_curr_variant

    @classmethod
    def _rewrite_value_if_special_header(cls, format_header, field_value):
        result = field_value  # by default, assume no special coddling needed

        if format_header == cls.CHR_HEADER:
            if field_value == cls.RAW_CHR_MT_VAL:
                result = cls.STANDARDIZED_CHR_MT_VAL
        elif format_header in [cls.THOU_G_2015_ALL_HEADER, cls.ESP6500_ALL_HEADER, cls.NCI60_HEADER]:
            result = float(field_value)
        elif format_header in [cls.START_HEADER, cls.END_HEADER]:
            result = int(field_value)
        # elif curr_header == cls.CYTOBAND_HEADER:
        #    cytoband_data = CytoBand(curr_value)
        #    result = cytoband_data.fill()
        elif format_header in [cls.GENOMIC_SUPERDUPS_HEADER, cls.TFBS_CONS_SITES_HEADER]:
            result = cls._parse_to_dict_with_score_key(field_value)
        # end checking if this is a field that needs special coddling

        return result

    @classmethod
    def _parse_to_dict_with_score_key(cls, delimited_str):
        """Parse delimited string of key/value pairs to a dictionary, with Score key's value cast to a float.

        Args:
            delimited_str (str): A string of key/value pairs, with each pair delimited from the next by a ';' and,
                with each pair, the key delimited from the value by an '='.  A key named 'Score' must be present.

        Returns:
            dict(str, Any): Dictionary of key/value pairs from the input string; the value for the 'Score' key is cast
                to a float.
        """
        key_val_pairs_as_str = delimited_str.split(";")
        result = dict(curr_key_val_pair_str.split("=") for curr_key_val_pair_str in key_val_pairs_as_str)
        result[cls.SCORE_KEY] = float(result[cls.SCORE_KEY])
        return result


class AnnovarAnnotatedVariant(object):
    GENOTYPE_KEY = 'genotype'
    FILTER_PASSING_READS_COUNT_KEY = 'filter_passing_reads_count'
    GENOTYPE_LIKELIHOODS_KEY = 'genotype_likelihoods'
    SAMPLES_KEY = 'samples'
    SAMPLE_ID_KEY = 'sample_id'
    GENOTYPE_SUBCLASS_BY_CLASS_KEY = 'genotype_subclass_by_class'
    HGVS_ID_KEY = 'hgvs_id'
    ALLELE_DEPTH_KEY = 'AD'

    @staticmethod
    def _list_has_valid_content(a_list):
        """Return true if list has at least one non-None entry, false otherwise.

        Args:
            a_list (List[Any]): The list to evaluate for validity.

        Returns:
            bool: True if list has at least one non-None entry, false otherwise.
        """

        is_valid = len(a_list) > 0 and not all(curr_item is None for curr_item in a_list)
        return is_valid

    @classmethod
    def make_per_variant_annotation_dict(cls, fields_by_annovar_header, hgvs_id, format_string,
                                         genotype_field_strings_by_sample_name):
        result = fields_by_annovar_header
        result[cls.HGVS_ID_KEY] = hgvs_id

        # parse sample-level info into a list of per-sample dicts and add that list to the variant-level dict
        sample_specific_dicts_list = []
        for curr_sample_name, genotype_fields_string_for_curr_sample in genotype_field_strings_by_sample_name.items():
            sample_specific_dict = cls._make_per_sample_annotation_dict(curr_sample_name, format_string,
                                                                        genotype_fields_string_for_curr_sample)

            if sample_specific_dict is not None:
                sample_specific_dicts_list.append(sample_specific_dict)

        result[cls.SAMPLES_KEY] = sample_specific_dicts_list
        return result

    @classmethod
    def _make_per_sample_annotation_dict(cls, sample_name, format_string, genotype_fields_string):

        if not VCFGenotypeParser.is_valid_genotype_fields_string(genotype_fields_string):
            return None

        genotype_info = VCFGenotypeParser.parse(format_string, genotype_fields_string)
        if genotype_info is None:
            return None

        # Always include a genotype key, EVEN IF the value for that key is None
        result = {cls.SAMPLE_ID_KEY: sample_name, cls.GENOTYPE_KEY: genotype_info.genotype}

        # for all other keys, include them only if they have meaningful values
        if genotype_info.genotype is not None:
            result[cls.GENOTYPE_SUBCLASS_BY_CLASS_KEY] = genotype_info.genotype_subclass_by_class
        if genotype_info.filter_passing_reads_count is not None:
            result[cls.FILTER_PASSING_READS_COUNT_KEY] = genotype_info.filter_passing_reads_count

        genotype_likelihoods_list = [i.likelihood_neg_exponent for i in genotype_info.genotype_likelihoods]
        if cls._list_has_valid_content(genotype_likelihoods_list):
            result[cls.GENOTYPE_LIKELIHOODS_KEY] = genotype_likelihoods_list

        unfiltered_read_depths_list = [i.unfiltered_read_counts for i in genotype_info.alleles]
        if cls._list_has_valid_content(unfiltered_read_depths_list):
            result[cls.ALLELE_DEPTH_KEY] = unfiltered_read_depths_list

        return result

    # # Gets its own class because it is particularly pesky to parse
    # class CytoBand(object):
    #     NAME_KEY = "Name"
    #     CHROMOSOME_KEY = "Chromosome"
    #     BAND_KEY = "Band"
    #     REGION_KEY = "Region"
    #     SUBBAND_KEY = "Sub_Band"
    #
    #     def __init__(self, cyto_band_name):
    #
    #         self.letters = set('XY')
    #         self.name = cyto_band_name
    #         self.processed = self.fill()
    #
    #     def fill(self):
    #
    #         processed = {self.NAME_KEY: self.name}
    #         spliced = re.split('(\D+)', self.name)  # \D means "any non-digit", so \D+ means 1 or more non-digit
    #
    #         if any((c in self.letters) for c in self.name):
    #             processed[self.CHROMOSOME_KEY] = spliced[1][0]
    #             processed[self.BAND_KEY] = spliced[1][1]
    #         else:
    #             processed[self.CHROMOSOME_KEY] = int(spliced[0])
    #             processed[self.BAND_KEY] = spliced[1]
    #
    #         processed[self.REGION_KEY] = spliced[2]
    #
    #         if '.' in spliced:
    #             processed[self.SUBBAND_KEY] = spliced[-1]
    #         return processed
