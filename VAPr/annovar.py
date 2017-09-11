import os
import sys
import shlex
import glob
import time
import subprocess
import logging
from collections import OrderedDict

# TODO: Is this really necessary?
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout     # Enables logging on jupyter notebooks
except:
    pass

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'


def listen(out_path, num_batch_jobs, num_files):
    """ Check for newly created annotated files """

    txt_extension = "txt"
    num_added_text_files = 0

    # TODO: Refactor to do one initial read before loop (to get initial files), then one read within each loop,
    # with state saved across loops (unlike here) instead of two reads every loop?
    while True:
        initial_text_file_names = []
        walker = os.walk(out_path)
        for folder, subfolders, files in walker:
            for _file in files:
                if _file.endswith(txt_extension):
                    initial_text_file_names.append(_file)

        time.sleep(5)
        # TODO: refactor to remove duplication of code
        later_text_file_names = []
        walker = os.walk(out_path)
        for folder, subfolders, files in walker:
            for _file in files:
                if _file.endswith(txt_extension):
                    later_text_file_names.append(_file)

        new_text_file_names = [x for x in later_text_file_names if x not in initial_text_file_names]

        for curr_new_text_file in new_text_file_names:
            num_added_text_files += 1
            logging.info('File %i/%i: Annovar finished working on file : ' % (num_added_text_files, num_files) +
                         os.path.basename(curr_new_text_file) +
                         '.\n A text file has been created in the %s directory\n' % out_path)
        # end for

        # TODO: Possible bug:
        # why check whether the number of text files *added this time around* matches the number of batch
        # jobs??  Why check if the number of *initial* text files this time around matches the expected total
        # number?Why not just check if the total number of text files ever seen (later_text_file_names)
        # matches
        if num_added_text_files == num_batch_jobs:
            break
        if len(initial_text_file_names) >= num_files:
            break


class AnnovarWrapper:
    """ Wrapper around ANNOVAR download and annotation functions """

    # TODO: Pull out string keys into symbolic constants
    hg_18_databases = OrderedDict({'knownGene': 'g',
                                   'tfbsConsSites': 'r',
                                   'cytoBand': 'r',
                                   'targetScanS': 'r',
                                   'genomicSuperDups': 'r',
                                   'esp6500siv2_all': 'f',
                                   '1000g2015aug': 'f',
                                   'cosmic70': 'f',
                                   'nci60': 'f'})

    hg_19_databases = OrderedDict({'knownGene': 'g',
                                   'tfbsConsSites': 'r',
                                   'cytoBand': 'r',
                                   'targetScanS': 'r',
                                   'genomicSuperDups': 'r',
                                   'esp6500siv2_all': 'f',
                                   '1000g2015aug': 'f',
                                   'popfreq_all_20150413': 'f',
                                   'clinvar_20161128': 'f',
                                   'cosmic70': 'f',
                                   'nci60': 'f'})

    hg_38_databases = OrderedDict({'knownGene': 'g',
                                   'cytoBand': 'r',
                                   'genomicSuperDups': 'r',
                                   'esp6500siv2_all': 'f',
                                   '1000g2015aug': 'f',
                                   'clinvar_20161128': 'f',
                                   'cosmic70': 'f',
                                   'nci60': 'f'})

    def __init__(self, input_dir, output_csv_path, annovar_path, mongo_db_and_collection_names_dict,
                 list_of_vcf_mapping_dicts, design_file=None, genome_build_version=None):

        self.HUMANDB_FOLDER_NAME = "/humandb/"

        # TODO: Figure out what this string actually *is*
        self.DOWN_DD = '-webfrom annovar'

        # TODO: It appears these properties are no longer used
        # self.dl_list_command = 'avdblist'
        # self.manual_update = {'clinvar_20161128': [datetime.datetime(2016, 11, 28)],
        #          '1000g2015aug': [datetime.datetime(2016, 8, 30)],
        #          'popfreq_all_20150413': [datetime.datetime(2015, 4, 13)]}

        # TODO: Pull out string keys into symbolic constants
        self.ANNOVAR_DB_IS_HOSTED_BY_ANNOVAR = OrderedDict({'knownGene': True,
                                                            'tfbsConsSites': False,
                                                            'cytoBand': False,
                                                            'targetScanS': False,
                                                            'genomicSuperDups': False,
                                                            'esp6500siv2_all': True,
                                                            '1000g2015aug': True,
                                                            'popfreq_all_20150413': True,
                                                            'clinvar_20161128': True,
                                                            'cosmic70': True,
                                                            'nci60': True,
                                                            'avdblist': True})

        # Project data
        self.input_dir = input_dir
        self.output_csv_path = output_csv_path
        self.annovar_path = annovar_path
        self.mongo_db_and_collection_names_dict = mongo_db_and_collection_names_dict
        self.design_file = design_file
        self.genome_build_version = genome_build_version
        self.list_of_vcf_mapping_dicts = list_of_vcf_mapping_dicts

        # Databases data
        self.annovar_dbs_to_use = self._get_annovar_dbs_to_use_for_build_version()

    # TODO: having both a "get_all_dbs" parameter and an "annovar_dbs_to_get" parameter introduces the potential for
    # unexpected behavior--if the user sets get_all_dbs to True (not understanding precisely what it implies) then any
    # value in annovar_dbs_to_get will be ignored without comment or error.  It would be preferable to have only
    # annovar_dbs_to_get and to set it to None if you wanted to get all annovar dbs.  However, I have not made that
    # change yet in order to maintain the public interface.
    def download_dbs(self, annovar_dbs_to_get=None):
        """
        Wrapper around annotate_variation.pl with -downdb as optional arg
        First, it cleans up the humandb/ directory to avoid conflicts, then gets newest versions of databases
        by spawning the jobs using subprocesses
        """

        # TODO: Do we really want to delete all the user's past annovar dbs without any warning??
        # Could make it very hard for them to reproduce any past results!
        if len(os.listdir(os.path.join(self.annovar_path, self.HUMANDB_FOLDER_NAME))) > 0:
            files = glob.glob(os.path.join(self.annovar_path, self.HUMANDB_FOLDER_NAME, '*'))
            for f in files:
                os.remove(f)

        list_commands = self._build_annovar_database_download_command_str(annovar_dbs_to_get)
        for command in list_commands:
            args = shlex.split(command)
            subprocess.Popen(args, stdout=subprocess.PIPE)

            logging.info('Finished downloading databases to {}'.format(
                os.path.join(self.annovar_path, self.HUMANDB_FOLDER_NAME)))

    def run_annovar(self, num_batch_jobs=10, vcf_is_multisample=False):
        """ Spawn ANNOVAR VCF annotation jobs in batches of five files at a time to prevent memory overflow """

        handler = AnnovarJobHandler(num_batch_jobs, vcf_is_multisample, self.list_of_vcf_mapping_dicts)
        num_files_created = 0
        for index, job in enumerate(handler.chunkenize):
            logging.info('Job %i/%i sent for processing' %
                         (index + 1, len(self.list_of_vcf_mapping_dicts) / num_batch_jobs + 1))
            num_files_created += len(job)

# break into separate function(s) and test for functionality
            for idx, _map in enumerate(job):
                annotation_dir = _map['csv_file_full_path']
                if os.path.isdir(annotation_dir):
                    logging.info('Directory already exists for %s. '
                                 'Writing output files there for file %s.' % (annotation_dir,
                                                                              _map['raw_vcf_file_full_path']))
                else:
                    os.makedirs(annotation_dir)

                vcf_path = _map['raw_vcf_file_full_path']
                csv_path = os.path.join(_map['csv_file_full_path'], _map['csv_file_basename'])
                cmd_string = self._build_table_annovar_command_str(vcf_path, csv_path,
                                                                   vcf_is_multisample=vcf_is_multisample)
                args = shlex.split(cmd_string)

                subprocess.Popen(args, stdout=subprocess.PIPE)

            logging.info('Annovar jobs submitted for %i files: %s' % (len(job),
                                                                      ', '.join([os.path.basename(
                                                                          i['raw_vcf_file_full_path']) for i in
                                                                                 job])))

            listen(self.output_csv_path, len(job), num_files_created)
            logging.info('Finished running Annovar on this batch')
        logging.info('Finished running Annovar on all files')

    def _build_table_annovar_command_str(self, vcf_path, csv_path, vcf_is_multisample=False):
        """Generate command string to run table_annovar.pl, which annotates a VCF file."""

        # TODO: check for newer version of databases
        dbs = ",".join(list(self.annovar_dbs_to_use.keys()))
        dbs_args = ",".join(list(self.annovar_dbs_to_use.values()))

        # TODO: Explain why this replacement is being done?
        if '1000g2015aug' in dbs:
            dbs = dbs.replace('1000g2015aug', '1000g2015aug_all')
        command = " ".join([
            'perl', os.path.join(self.annovar_path, 'table_annovar.pl'),
            vcf_path, "".join([self.annovar_path, self.HUMANDB_FOLDER_NAME]),
            '-genome_build_version', self.genome_build_version,
            '-out', csv_path, '-remove -protocol', dbs,
            '-operation', dbs_args, '-nastring .', '-otherinfo -vcfinput'
        ])

        if vcf_is_multisample:
            command += ' -format vcf4 -allsample -withfreq'

        return command

    def _build_annovar_database_download_command_str(self, annovar_dbs_to_get):
        """Concatenate command string arguments for Annovar download database jobs"""

        if annovar_dbs_to_get is not None:
            for annovar_db_name in annovar_dbs_to_get:
                if annovar_db_name not in self.annovar_dbs_to_use:
                    raise ValueError('Database %s not supported for build version %s' % (annovar_db_name,
                                                                                         self.genome_build_version))
                self.annovar_dbs_to_use = {db: self.annovar_dbs_to_use[db] for db in annovar_dbs_to_get}

        command_list = []
        # TODO: refactor to remove duplicated command components :(
        for annovar_db_name in self.annovar_dbs_to_use:
            if self.ANNOVAR_DB_IS_HOSTED_BY_ANNOVAR[annovar_db_name]:
                command_list.append(" ".join(
                    ['perl', self.annovar_path + '/annotate_variation.pl',
                     '-build', self.genome_build_version,
                     '-downdb', self.DOWN_DD, annovar_db_name,
                     self.annovar_path + self.HUMANDB_FOLDER_NAME]))
            else:
                command_list.append(" ".join(
                    ['perl', self.annovar_path + '/annotate_variation.pl',
                     '-build', self.genome_build_version,
                     '-downdb', annovar_db_name,
                     self.annovar_path + self.HUMANDB_FOLDER_NAME]))

        return command_list

    # TODO: This method is deprecated, and it appears it is unused
    # def check_for_database_updates(self):
    #     """ Deprecated """
    #
    #     self.download_dbs(all_dbs=False, dbs=['avdblist'])
    #
    #     with open(os.path.join(self.annovar, '/humandb' + self.genome_build_version + '_avdblist.txt'), 'r') as db_list:
    #         reader = csv.reader(db_list, delimiter='\t')
    #         db_dict = {}
    #         for i in reader:
    #             if i[0][-6:] != 'idx.gz':
    #                 db_dict[i[0][5:]] = [datetime.datetime(int(i[1][0:4]), int(i[1][4:6]), int(i[1][6:8])), i[2]]
    #
    #     for db in self.manual_update.keys():
    #         for db_ in db_dict.keys():
    #             if db_.startswith(db):
    #                 if db_dict[db_][0] > self.manual_update[db][0]:
    #                     logging.info('Database %s outdated, will download newer version' % db_)
    #                     self.download_dbs(all_dbs=False, dbs=[os.path.splitext(os.path.splitext(db_)[0])[0]])

    def _get_annovar_dbs_to_use_for_build_version(self):
        if self.genome_build_version == 'hg18':
            databases = AnnovarWrapper.hg_18_databases
        elif self.genome_build_version == 'hg19':
            databases = AnnovarWrapper.hg_19_databases
        elif self.genome_build_version == 'hg38':
            databases = AnnovarWrapper.hg_38_databases
        else:
            raise NameError('Genome {0} not supported by VAPr, genome must be one of the following: hg18, hg19, or hg38'.format(self.genome_build_version))

        return databases


class AnnovarJobHandler:
    def __init__(self, num_batch_jobs, multisample, list_of_vcf_mapping_dicts):

        self.list_of_vcf_mapping_dicts = list_of_vcf_mapping_dicts

        # TODO: It does not appear that this property is ever used; if not, take out of init params
        # self.multisample = multisample

        self.num_batch_jobs = num_batch_jobs
        # TODO: This should probably be changed to integer division, as division to produce a float
        # causes linter to be concerned when value is used as step increment in for loop below
        if self.num_batch_jobs > len(self.list_of_vcf_mapping_dicts):
            self.num_batch_jobs = len(self.list_of_vcf_mapping_dicts) // 2

        # TODO: Do we really need to assign _get_chunk_of_mappings_to_process to a new property?
        # After all, it is already a method.  Couldn't we use it directly in _next?
        self.chunkenize = self._get_chunk_of_mappings_to_process()

    def _get_chunk_of_mappings_to_process(self):
        """Yield successive n-sized chunks of vcf info mappings to process."""
        for i in range(0, len(self.list_of_vcf_mapping_dicts), self.num_batch_jobs):
            yield self.list_of_vcf_mapping_dicts[i:i + self.num_batch_jobs]

    def _next(self):
        return next(self.chunkenize)


# if __name__ == '__main__':
#
#     def download_dbs(all_dbs=True, dbs=None):
#         """
#         Implementation of the wrapper around annotate_variation.pl with -downdb as optional arg
#         First, it cleans up the humandb/ directory to avoid conflicts, then gets newest versions of databases
#         by spawning the jobs using subprocesses
#
#         """
#
#         if len(os.listdir(os.path.join(annovar, 'humandb/'))) > 0:
#             files = glob.glob(os.path.join(annovar, 'humandb/*'))
#             for f in files:
#                 os.remove(f)
#
#         list_commands = self._build_annovar_database_download_command_str(all_dbs, dbs)
#         for command in list_commands:
#             args = shlex.split(command)