import myvariant
import vcf
import csv
import re
import itertools
import logging
import sys
import vcf_parsing as vvp
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class HgvsParser(object):

    def __init__(self, vcf_file):

        self.vcf = vcf_file
        self.num_lines = sum(1 for _ in open(self.vcf))
        self.chunksize = 950

    def get_variants_from_vcf(self, step):
        """
        Retrieves variant names from a LARGE vcf file.
        :param step: ...
        :return: a list of variants formatted according to HGVS standards
        """
        list_ids = []
        reader = vcf.Reader(open(self.vcf, 'r'))

        for record in itertools.islice(reader, step * self.chunksize, (step + 1) * self.chunksize):
            if len(record.ALT) > 1:
                for alt in record.ALT:
                    list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS,
                                                          record.REF, str(alt)))
            else:
                list_ids.append(myvariant.format_hgvs(record.CHROM, record.POS,
                                                      record.REF, str(record.ALT[0])))

        return self.complete_chromosome(list_ids)

    @staticmethod
    def complete_chromosome(expanded_list):
        for i in range(0, len(expanded_list)):
            if 'M' in expanded_list[i]:
                one = expanded_list[i].split(':')[0]
                two = expanded_list[i].split(':')[1]
                if 'MT' not in one:
                    one = 'chrMT'
                expanded_list[i] = one + ':' + two
        return expanded_list


class TxtParser(object):

    def __init__(self, txt_file):

        self.txt_file = txt_file
        self.num_lines = sum(1 for _ in open(self.txt_file))
        self.chunksize = 950
        self.offset = 0
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

        listofdicts = []
        with open(self.txt_file, 'r') as txt:

            reader = csv.reader(txt, delimiter='\t')
            header = self._normalize_header(next(reader))

            for i in itertools.islice(reader, (step*self.chunksize) + self.offset,
                                      ((step+1)*self.chunksize) + offset + self.offset):

                sparse_dict = dict(zip(header[0:len(header)-1], i[0:len(header)-1]))
                sparse_dict['otherinfo'] = i[-2::]

                if build_ver == 'hg19':
                    dict_filled = {k: sparse_dict[k] for k in self.hg_19_columns if sparse_dict[k] != '.'}
                elif build_ver == 'hg18':
                    dict_filled = {k: sparse_dict[k] for k in self.hg_18_columns if sparse_dict[k] != '.'}
                else:
                    dict_filled = {k: sparse_dict[k] for k in self.hg_38_columns if sparse_dict[k] != '.'}

                modeled = AnnovarModels(dict_filled)
                listofdicts.append(modeled.final_dict)

            self.offset += offset

        return listofdicts

    @staticmethod
    def _normalize_header(header):
        normalized = []

        for item in header:
            normalized.append(item.lower().replace('.', '_'))

        return normalized


class AnnovarModels(object):

    def __init__(self, dictionary):

        self.dictionary = dictionary
        self.existing_keys = self.dictionary.keys()
        self.errors = None
        self.final_dict = self.process()

    def process(self):
        for key in self.dictionary.keys():

            if self.dictionary['chr'] == 'chrM':
                self.dictionary['chr'] = 'chrMT'

            if key in ['1000g2015aug_all', 'esp6500si_all', 'nci60']:
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

        self.dictionary['genotype'], self.errors = self.parse_genotype()

        return self.dictionary

    def parse_genotype(self):

        read_depth_error = genotype_lik_error = allele_error = 0
        parser = vvp.VCFGenotypeStrings()

        genotype_to_fill = parser.parse(self.dictionary['otherinfo'][0], self.dictionary['otherinfo'][1])

        if not genotype_to_fill:
            print('skipped')
            return None
        gen_dic = {'genotype': genotype_to_fill.genotype}

        try:
            gen_dic['filter_passing_reads_count'] = [genotype_to_fill.filter_passing_reads_count]
        except ValueError:
            read_depth_error += 1
            # print('Read depth returned null value')

        try:
            gen_dic['genotype_likelihoods'] = [float(genotype_to_fill.genotype_likelihoods[0].likelihood_neg_exponent),
                                               float(genotype_to_fill.genotype_likelihoods[1].likelihood_neg_exponent),
                                               float(genotype_to_fill.genotype_likelihoods[2].likelihood_neg_exponent)]
        except IndexError:
            genotype_lik_error += 1

        try:
            gen_dic['alleles'] = [genotype_to_fill.alleles[0].read_counts,
                                  genotype_to_fill.alleles[1].read_counts]
        except IndexError or ValueError:
            allele_error += 1
        errors = (read_depth_error, genotype_lik_error, allele_error)
        return gen_dic, errors

    @staticmethod
    def _to_int(val):
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
