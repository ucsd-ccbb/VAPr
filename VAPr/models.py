import myvariant
import vcf
import csv
import re
import itertools
import logging
import sys
import VAPr.vcf_parsing as vvp
import VAPr.definitions as definitions
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class HgvsParser(object):

    """ Class that process vcf files and extracts their hgvs ids"""

    def __init__(self, vcf_file):

        self.vcf = vcf_file
        # self.num_lines = sum(1 for _ in open(self.vcf))
        self.chunksize = definitions.chunk_size
        self.samples = vcf.Reader(open(self.vcf, 'r')).samples
        self.num_samples = len(self.samples)

    def get_num_lines(self):
        return sum(1 for _ in open(self.vcf))

    def get_all_variants_from_vcf(self):
        """ For debugging mostly """

        list_ids = []

        reader = vcf.Reader(open(self.vcf, 'r'))
        for record in reader:
            list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS,
                                                  record.REF, str(record.ALT[0])))

        return self.complete_chromosome(list_ids)

    def get_variants_from_vcf(self, step):
        """
        Retrieves variant names from a LARGE vcf file.
        :param step: tells the parallel processing where to start and end the hgvs id creation
        :return: a list of variants formatted according to HGVS standards
        """
        reader = vcf.Reader(open(self.vcf, 'r'))
        list_ids = []

        for record in itertools.islice(reader, step * self.chunksize, (step + 1) * self.chunksize):
            list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS,
                                                  record.REF, str(record.ALT[0])))

        return self.complete_chromosome(list_ids)

    @staticmethod
    def complete_chromosome(expanded_list):
        """ Ensuring syntax consistency """

        for i in range(0, len(expanded_list)):
            if 'M' in expanded_list[i]:
                one = expanded_list[i].split(':')[0]
                two = expanded_list[i].split(':')[1]
                if 'MT' not in one:
                    one = 'chrMT'
                expanded_list[i] = "".join([one, ':', two])
        return expanded_list


class TxtParser(object):

    """ Class that process an Annovar created csv file """

    def __init__(self, txt_file, samples=None, extra_data=None):

        self.samples = samples
        self.txt_file = txt_file
        self.num_lines = sum(1 for _ in open(self.txt_file))
        self.chunksize = definitions.chunk_size
        self.offset = 0
        self.extra_data = extra_data
        self.hg_19_columns = ['chr',
                              'start',
                              'end',
                              'ref',
                              'alt',
                              'func_knowngene',
                              'gene_knowngene',
                              'genedetail_knowngene',
                              'exonicfunc_knowngene',
                              'tfbsconssites',
                              'cytoband',
                              'genomicsuperdups',
                              '1000g2015aug_all',
                              'esp6500siv2_all',
                              'cosmic70',
                              'nci60',
                              'otherinfo']

        self.hg_18_columns = self.hg_19_columns

        self.hg_38_columns = ['chr',
                              'start',
                              'end',
                              'ref',
                              'alt',
                              'func_knowngene',
                              'gene_knowngene',
                              'genedetail_knowngene',
                              'exonicfunc_knowngene',
                              'cytoband',
                              'genomicsuperdups',
                              '1000g2015aug_all',
                              'esp6500siv2_all',
                              'cosmic70',
                              'nci60',
                              'otherinfo']

    def open_and_parse_chunks(self, step, build_ver=None, offset=0):

        """
        Parsing function that retrieves data from a specific location from a csv file.
        It will process every (significant) data field in an annovar-annotated csv and
        conveniently return it into a dictionary

        :param step: tells the parallel processing where to start and end the hgvs id creation
        :param build_ver: genome build version
        :param offset: deprecated, optional argument may be used later. Used to spot discrepancies between
                       csv and vcf file
        :return: list of dictionaries, where each dictionary is a variant document
        """

        listofdicts = []

        with open(self.txt_file, 'r') as txt:

            reader = csv.reader(txt, delimiter='\t')
            header = self._normalize_header(next(reader))

            for i in itertools.islice(reader, (step*self.chunksize) + self.offset,
                                      ((step+1)*self.chunksize) + offset + self.offset):

                sparse_dict = dict(zip(header[0:len(header)-1], i[0:len(header)-1]))
                sparse_dict['otherinfo'] = i[-1-len(self.samples)::]

                if build_ver == 'hg19':
                    dict_filled = {k: sparse_dict[k] for k in self.hg_19_columns if sparse_dict[k] != '.'}
                elif build_ver == 'hg18':
                    dict_filled = {k: sparse_dict[k] for k in self.hg_18_columns if sparse_dict[k] != '.'}
                else:
                    dict_filled = {k: sparse_dict[k] for k in self.hg_38_columns if sparse_dict[k] != '.'}

                modeled = AnnovarModels(dict_filled, self.samples, extra_data=self.extra_data)
                listofdicts.append(modeled.final_list_dict)

            self.offset += offset

        return listofdicts

    @staticmethod
    def _normalize_header(header):
        """ Ensuring syntax consistency """

        normalized = []

        for item in header:
            normalized.append(item.lower().replace('.', '_'))

        return normalized


class AnnovarModels(object):

    """
    The actuall class that process a single annovar-annotated variant.
    It performs specific and necessary manipulations to ensure that the data will be parsed appropriately
    and can be queried easily from MongoDB. Basically our ODM (Object Data Model) for the Annovar data.
    """

    def __init__(self, dictionary, samples, extra_data=None):

        if extra_data:
            self.dictionary = dict(extra_data.items() + dictionary.items())
        else:
            self.dictionary = dictionary
        self.samples = samples
        self.existing_keys = self.dictionary.keys()
        self.errors = None
        self.extra_data = extra_data
        self.final_list_dict = self.process()

    def process(self):
        for key in self.dictionary.keys():

            if self.dictionary['chr'] == 'chrM':
                self.dictionary['chr'] = 'chrMT'

            if key in ['1000g2015aug_all', 'esp6500siv2_all', 'nci60']:
                self.dictionary[key] = float(self.dictionary[key])

            if key in ['start', 'end']:
                self.dictionary[key] = int(self.dictionary[key])

            if key == 'cytoband':
                cytoband_data = CytoBand(self.dictionary['cytoband'])
                self.dictionary['cytoband'] = cytoband_data.fill()

            if key in ['genomicsuperdups', 'tfbsconssites']:
                self.dictionary[key] = self.to_dict(key)

            if key == 'otherinfo':
                self.dictionary[key] = [i for i in self.dictionary[key] if i != '.']

        final_annovar_list_of_dicts, self.errors = self.parse_genotype(self.dictionary)

        return final_annovar_list_of_dicts

    def parse_genotype(self, dictionary):
        """ Implements the genotype parsing scheme. Many thanks to Amanda Birmingham """

        read_depth_error = genotype_lik_error = allele_error = 0
        parser = vvp.VCFGenotypeStrings()

        list_dictionaries = []

        for index, sample in enumerate(self.samples):
            sample_specific_dict = {k: v for k, v in dictionary.items()}  # make copy, propagate genotype info over alleles

            genotype_to_fill = parser.parse(dictionary['otherinfo'][0],
                                            dictionary['otherinfo'][index + 1])

            sample_specific_dict['sample_id'] = sample
            sample_specific_dict['genotype'] = genotype_to_fill.genotype

            try:
                sample_specific_dict['filter_passing_reads_count'] = [float(genotype_to_fill.filter_passing_reads_count)]
            except ValueError:
                read_depth_error += 1
            except TypeError:
                read_depth_error += 1
            try:
                sample_specific_dict['genotype_likelihoods'] = [float(i.likelihood_neg_exponent) for i in
                                                                genotype_to_fill.genotype_likelihoods]
            except IndexError:
                genotype_lik_error += 1

            try:
                sample_specific_dict['alleles'] = [float(i.read_counts) for i in genotype_to_fill.alleles]
            except IndexError:
                allele_error += 1
            except ValueError:
                allele_error += 1

            if sample_specific_dict['genotype'] is not None:
                genotype_class_to_fill = vvp.VCFGenotypeInfo('')
                vvp.fill_genotype_class(sample_specific_dict['genotype'], genotype_class_to_fill)
                sample_specific_dict['genotype_subclass_by_class'] = genotype_class_to_fill.genotype_subclass_by_class

            list_dictionaries.append(sample_specific_dict)

        errors = (read_depth_error, genotype_lik_error, allele_error)

        return list_dictionaries, errors

    @staticmethod
    def _to_int(val):
        """ This is.... Bad """
        try:
            val = int(val)
        except:
            pass
        return val

    def to_dict(self, key):
        as_dict = dict(item.split("=") for item in self.dictionary[key].split(";"))
        as_dict["Score"] = float(as_dict["Score"])
        return as_dict


class CytoBand(object):

    """ Gets its own class cause it is particularly pesky to parse"""

    def __init__(self, cyto_band_name):

        self.letters = set('XY')
        self.name = cyto_band_name
        self.processed = self.fill()

    def fill(self):

        processed = {'Name': self.name}
        spliced = re.split('(\D+)', self.name)

        if any((c in self.letters) for c in self.name):
            processed['Chromosome'] = spliced[1][0]
            processed['Band'] = spliced[1][1]
        else:
            processed['Chromosome'] = int(spliced[0])
            processed['Band'] = spliced[1]

        processed['Region'] = spliced[2]

        if '.' in spliced:
            processed['Sub_Band'] = spliced[-1]
        return processed

