# built-in libraries
import csv
import itertools
import logging
import sys

# third-party libraries
import myvariant

# project libraries
import VAPr.vcf_genotype_fields_parsing

# TODO: Understand, vet this logging set-up
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


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

    ANNOVAR_OUTPUT_COLS = [CHR_HEADER,
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
                           THOU_G_2015_ALL_HEADER,
                           # ESP6500_ALL_HEADER,
                           # 'cosmic70',
                           # NCI60_HEADER,
                           OTHERINFO_HEADER]

    @staticmethod
    def _normalize_header(raw_headers_list):
        """Lower-case all header strings and replace periods with underscores."""

        normalized_headers_list = []

        for curr_header in raw_headers_list:
            normalized_headers_list.append(curr_header.lower().replace('.', '_'))

        return normalized_headers_list

    @classmethod
    def mine_chunk_of_annovar_annotations(cls, annovar_txt_file_path, sample_names_list, chunk_number, chunk_size):
        annotations_dict_per_variant_list = []
        hgvsid_list = []

        with open(annovar_txt_file_path, 'r') as txt:
            # read in the header line and normalize the header names
            reader = csv.reader(txt, delimiter='\t')
            header = cls._normalize_header(next(reader))

            # for each row in this chunk--which is to say, each variant
            for curr_line_fields in itertools.islice(reader, (chunk_number * chunk_size),
                                                     ((chunk_number + 1) * chunk_size)):
                # make a dictionary of all the fields in this line keyed by all the headers
                last_field_index = len(header) - 1
                all_fields_dict = dict(zip(header[0:last_field_index], curr_line_fields[0:last_field_index]))

                # TODO: figure out this slice again and make in-line comment to explain
                all_fields_dict[cls.OTHERINFO_HEADER] = curr_line_fields[-1 - len(sample_names_list)::]

                # drop out any key/value pairs for fields that are empty
                filled_fields_dict = {k: all_fields_dict[k] for k in cls.ANNOVAR_OUTPUT_COLS
                                      if all_fields_dict[k] != '.'}

                # generate the hgvs id for this variant
                hgvs_id = myvariant.format_hgvs(filled_fields_dict[cls.CHR_HEADER],
                                                filled_fields_dict[cls.START_HEADER],
                                                filled_fields_dict[cls.REF_HEADER],
                                                filled_fields_dict[cls.ALT_HEADER])
                hgvsid_list.append(hgvs_id)

                # turn the dictionary of annovar fields into a dictionary of annotations for the variant, including
                # nested structures containing sample-specific genotype-related info
                annotations_dict_for_curr_variant = AnnovarAnnotatedVariant.make_annotations_dict(
                    hgvs_id, filled_fields_dict, sample_names_list)
                annotations_dict_per_variant_list.append(annotations_dict_for_curr_variant)

        return annotations_dict_per_variant_list, hgvsid_list


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
    def to_dict(delimited_str):
        as_dict = dict(item.split("=") for item in delimited_str.split(";"))
        as_dict["Score"] = float(as_dict["Score"])
        return as_dict

    @staticmethod
    def _list_has_valid_content(a_list):
        # if there are entries in the list and they are not ALL none (ok if some are none):
        is_valid = len(a_list) > 0 and not all(curr_item is None for curr_item in a_list)
        return is_valid

    @classmethod
    def make_annotations_dict(cls, hgvs_id, fields_by_annovar_header, sample_names_list):
        # result = {} if extra_data is None else dict(extra_data.items())
        result = {cls.HGVS_ID_KEY: hgvs_id}

        # first, loop over the keys, fix up any entries whose keys need special coddling, and create new clean dict
        for curr_header, curr_value in fields_by_annovar_header.items():
            if curr_header == AnnovarTxtParser.CHR_HEADER:
                if curr_value == 'chrM':
                    curr_value = 'chrMT'
            elif curr_header in [AnnovarTxtParser.THOU_G_2015_ALL_HEADER, AnnovarTxtParser.ESP6500_ALL_HEADER,
                                 AnnovarTxtParser.NCI60_HEADER]:
                curr_value = float(curr_value)
            elif curr_header in [AnnovarTxtParser.START_HEADER, AnnovarTxtParser.END_HEADER]:
                curr_value = int(curr_value)
            # elif curr_header == AnnovarTxtParser.CYTOBAND_HEADER:
            #    cytoband_data = CytoBand(curr_value)
            #    curr_value = cytoband_data.fill()
            elif curr_header in [AnnovarTxtParser.GENOMIC_SUPERDUPS_HEADER, AnnovarTxtParser.TFBS_CONS_SITES_HEADER]:
                curr_value = cls.to_dict(curr_header)
            # end checking if this is a field that needs special coddling

            result[curr_header] = curr_value
        # end loop over input dict to create clean dict

        # then parse the genotype info out of the cleaned-up dictionary
        sample_specific_dicts_list = cls._generate_sample_specific_dicts_list(result, sample_names_list)
        result[cls.SAMPLES_KEY] = sample_specific_dicts_list

        # remove the otherinfo field--we have parsed the info out of it, so it is now redundant
        result.pop(AnnovarTxtParser.OTHERINFO_HEADER, None)

        return result

    @classmethod
    def _generate_sample_specific_dicts_list(cls, fields_by_annovar_header, sample_names_list):
        format_string = fields_by_annovar_header[AnnovarTxtParser.OTHERINFO_HEADER][0]
        sample_specific_dicts_list = []

        for index, curr_sample_name in enumerate(sample_names_list):
            sample_specific_dict = {cls.SAMPLE_ID_KEY: curr_sample_name}

            genotype_fields_string_for_curr_sample = fields_by_annovar_header[AnnovarTxtParser.OTHERINFO_HEADER][
                index + 1]
            # TODO: need to handle more complexity here: '.' OR a string like './.:.<etc>'
            if genotype_fields_string_for_curr_sample == '.':
                continue

            genotype_info = VAPr.vcf_genotype_fields_parsing.VCFGenotypeParser.parse(
                format_string, genotype_fields_string_for_curr_sample)

            # Always include a genotype key, EVEN IF the value for that key is None
            sample_specific_dict[cls.GENOTYPE_KEY] = genotype_info.genotype

            # for all other keys, include them only if they have meaningful values
            if genotype_info.genotype is not None:
                sample_specific_dict[cls.GENOTYPE_SUBCLASS_BY_CLASS_KEY] = genotype_info.genotype_subclass_by_class
            if genotype_info.filter_passing_reads_count is not None:
                sample_specific_dict[cls.FILTER_PASSING_READS_COUNT_KEY] = genotype_info.filter_passing_reads_count

            genotype_likelihoods_list = [i.likelihood_neg_exponent for i in genotype_info.genotype_likelihoods]
            if cls._list_has_valid_content(genotype_likelihoods_list):
                sample_specific_dict[cls.GENOTYPE_LIKELIHOODS_KEY] = genotype_likelihoods_list

            unfiltered_read_depths_list = [i.unfiltered_read_counts for i in genotype_info.alleles]
            if cls._list_has_valid_content(unfiltered_read_depths_list):
                sample_specific_dict[cls.ALLELE_DEPTH_KEY] = unfiltered_read_depths_list

            sample_specific_dicts_list.append(sample_specific_dict)

        return sample_specific_dicts_list

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

    # class HgvsParser(object):
    #     """ Process a vcf file and extract its hgvs ids."""
    #
    #     @staticmethod
    #     def complete_chromosome(expanded_list):
    #         """ Ensuring syntax consistency """
    #
    #         for i in range(0, len(expanded_list)):
    #             if 'M' in expanded_list[i]:
    #                 one = expanded_list[i].split(':')[0]
    #                 two = expanded_list[i].split(':')[1]
    #                 if 'MT' not in one:
    #                     one = 'chrMT'
    #                 expanded_list[i] = "".join([one, ':', two])
    #         return expanded_list
    #
    # def __init__(self, vcf_file):
    #     self.vcf = vcf_file
    #     self.chunksize = definitions.chunk_size
    #     self.samples = vcf.Reader(open(self.vcf, 'r')).samples
    #     self.num_samples = len(self.samples)
    #
    #
    # # For debugging mostly
    # def get_all_variants_from_vcf(self):
    #     list_ids = []
    #
    #     reader = vcf.Reader(open(self.vcf, 'r'))
    #     for record in reader:
    #         list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS,
    #                                               record.REF, str(record.ALT[0])))
    #
    #     return self.complete_chromosome(list_ids)
    #
    # def get_variants_from_vcf(self, step):
    #     """ Retrieve variant names from a LARGE vcf file.
    #
    #     :param step: tells the parallel processing where to start and end the hgvs id creation
    #     :return: a list of variants formatted according to HGVS standards
    #     """
    #     reader = vcf.Reader(open(self.vcf, 'r'))
    #     list_ids = []
    #
    #     for record in itertools.islice(reader, step * self.chunksize, (step + 1) * self.chunksize):
    #         list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS,
    #                                               record.REF, str(record.ALT[0])))
    #
    #     return self.complete_chromosome(list_ids)
