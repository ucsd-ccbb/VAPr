from __future__ import division, print_function
import vcf
import myvariant
import csv
import re
import os
import itertools
import sys
from multiprocessing import Pool
from pymongo import MongoClient
from VAPr.base import AnnotationProject
import VAPr.vcf_parsing as vvp
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
logger.handlers[0].stream = sys.stdout


class VariantParsing(AnnotationProject):

    def __init__(self,
                 input_dir,
                 output_csv_path,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None):

        super(AnnotationProject, self).__init__(input_dir,
                                                output_csv_path,
                                                annovar_path,
                                                project_data,
                                                design_file=design_file,
                                                build_ver=build_ver)

        self.chunksize = 950
        self.step = 0
        self.csvs, self.vcfs = self.get_file_names()
        self.collection = project_data['project_name']
        self.db = project_data['db_name']
        self._buffer_len = 50000
        self._last_round = False

    def annotate_and_save(self, buffer=False):

        for csv, vcf in list(zip(self.csvs, self.vcfs)):

            csv_parsing = TxtParser(csv)
            hgvs = HgvsParser(vcf)
            sample_id = os.path.splitext(os.path.basename(vcf))[0]

            if not csv:

                while hgvs.num_lines > self.step * self.chunksize:

                    list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
                    myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, sample_id)

                    if len(myvariants_variants) < self.chunksize:
                        self._last_round = True

                    if self._last_round:
                        return 'Done'
                    else:
                        self.export(myvariants_variants)
                        self.step += 1

                return 'Done'

            else:

                if buffer:
                    variant_buffer = []
                    while csv_parsing.num_lines > self.step * self.chunksize:

                        list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
                        myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, sample_id)
                        offset = len(list_hgvs_ids) - self.chunksize
                        csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver, offset=offset)

                        merged_list = []
                        for i, _ in enumerate(myvariants_variants):
                            merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

                        variant_buffer.extend(merged_list)
                        self.step += 1

                        if len(merged_list) < self.chunksize:
                            self._last_round = True

                        if (len(variant_buffer) > self._buffer_len) or self._last_round:
                            logging.info('Parsing Buffer...')
                            self.export(variant_buffer)
                            variant_buffer = []

                            if self._last_round:
                                return 'Done2'

                    return 'Done1'

                else:

                    while csv_parsing.num_lines > self.step*self.chunksize:

                        list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
                        myvariants_variants = self.get_dict_myvariant(list_hgvs_ids, sample_id)
                        offset = len(list_hgvs_ids) - self.chunksize
                        csv_variants = csv_parsing.open_and_parse_chunks(self.step, offset=offset)

                        merged_list = []
                        for i, _ in enumerate(myvariants_variants):
                            merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

                        if len(merged_list) < self.chunksize:
                            self._last_round = True

                        if self._last_round:
                            return 'Done'
                        else:
                            self.export(merged_list)
                            self.step += 1

                return 'Done'

    def get_file_names(self):
        vcfs = [i[0] for i in self.mapping]
        csvs = []
        all_csvs = os.listdir(self.output_csv_path)

        for csv in all_csvs:
            for csv_init in [i[1] for i in self.mapping]:
                if csv.startswith(csv_init):
                    csvs.append(csv)

        return vcfs, csvs

    def check_ver(self, build_ver):
        """ Checking user input validity for genome build version """

        if not build_ver:
            return 'hg19'  # Default genome build vesion

        if build_ver not in self.supported_build_vers:
            raise ValueError('Build version must not recognized. Supported builds are'
                             ' %s, %s, %s' % (self.supported_build_vers[0],
                                              self.supported_build_vers[1],
                                              self.supported_build_vers[2]))
        else:
            return build_ver

    def export(self, list_docs):
        """
        Export data do a MongoDB server
        :param list_docs: list of dictionaries containing variant information
        :return: null
        """
        client = MongoClient()
        db = getattr(client, self.db)
        collection = getattr(db, self.collection)
        collection.insert_many(list_docs, ordered=False)

    @staticmethod
    def merge_dict_lists(*dict_args):
        """
        Given any number of dicts, shallow copy and merge into a new dict,
        precedence goes to key value pairs in latter dicts.
        """
        result = {}
        for dictionary in dict_args:
            result.update(dictionary)
        return result

    def get_dict_myvariant(self, variant_list, sample_id):
        """
        Function designated to place the queries on myvariant.info servers.

        :param variant_list: list of HGVS variant ID's. Usually retrived beforehand using the method
        get_variants_from_vcf
        from the class VariantParsing.
        :return: list of dictionaries. Each dictionary contains data about a single variant.
        """

        mv = myvariant.MyVariantInfo()
        # This will retrieve a list of dictionaries
        variant_data = mv.getvariants(variant_list, as_dataframe=False)
        variant_data = self.remove_id_key(variant_data, sample_id)

        return variant_data

    def remove_id_key(self, variant_data, sample_id):

        for dic in variant_data:
            dic['hgvs_id'] = dic.pop("_id", None)
            dic['hgvs_id'] = dic.pop("query", None)
            dic['sample_id'] = sample_id

        return variant_data

    def parallel_annotation(self, n_processes):

        mapping_tpls = [(k, v) for k, v in self.mapping.items()]
        self.pooling(n_processes, self._variant_parsing, mapping_tpls)

    def _variant_parsing(self, mapping):

        hgvs = HgvsParser(mapping[0])
        csv_parsing = TxtParser(mapping[1])

        variant_buffer = []
        while csv_parsing.num_lines > self.step * self.chunksize:

            list_hgvs_ids = hgvs.get_variants_from_vcf(self.step)
            myvariants_variants = self.get_dict_myvariant(list_hgvs_ids)
            offset = len(list_hgvs_ids) - self.chunksize
            csv_variants = csv_parsing.open_and_parse_chunks(self.step, build_ver=self.buildver, offset=offset)

            merged_list = []
            for i, _ in enumerate(myvariants_variants):
                merged_list.append(self.merge_dict_lists(myvariants_variants[i], csv_variants[i]))

            variant_buffer.extend(merged_list)
            self.step += 1

            if len(merged_list) < self.chunksize:
                self._last_round = True

            if (len(variant_buffer) > self._buffer_len) or self._last_round:
                logging.info('Parsing Buffer...')
                self.export(variant_buffer)
                variant_buffer = []

                if self._last_round:
                    return 'Done2'

        return 'Done1'

    @staticmethod
    def pooling(n_processes, input_1, input_2):

        pool = Pool(processes=n_processes)
        pool.map(input_1, input_2)
        pool.close()
        pool.join()


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
                if build_ver == 'hg19' or 'hg18':
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

        self.dictionary['genotype'] = self.parse_genotype()

        return self.dictionary

    def parse_genotype(self):

        parser = vvp.VCFGenotypeStrings()
        genotype_to_fill = parser.parse(self.dictionary['otherinfo'][0], self.dictionary['otherinfo'][1])
        gen_dic = {'genotype': genotype_to_fill.genotype,
                   'filter_passing_reads_count': [self._to_int(genotype_to_fill.filter_passing_reads_count)],
                   'genotype_likelihoods': [float(genotype_to_fill.genotype_likelihoods[0].likelihood_neg_exponent),
                                            float(genotype_to_fill.genotype_likelihoods[1].likelihood_neg_exponent),
                                            float(genotype_to_fill.genotype_likelihoods[2].likelihood_neg_exponent)],
                   'alleles': [self._to_int(genotype_to_fill.alleles[0].read_counts),
                               self._to_int(genotype_to_fill.alleles[1].read_counts)]
                   }

        return gen_dic

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
