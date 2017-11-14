import os
import sys
import shlex
import subprocess
import logging
from collections import OrderedDict

# TODO: Understand, vet this logging set-up
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout  # Enables logging on jupyter notebooks
except:
    pass


class AnnovarWrapper(object):
    """ Wrapper around ANNOVAR download and annotation functions """

    # TODO: someday: Refactor to dict, so there is guaranteed to be a list of dbs for each supported build version?
    hg_19_databases = OrderedDict({'knownGene': 'g',
                                   # 'tfbsConsSites': 'r',
                                   # 'cytoBand': 'r',
                                   # 'targetScanS': 'r',
                                   # 'genomicSuperDups': 'r',
                                   # 'esp6500siv2_all': 'f',
                                   '1000g2015aug': 'f'})
    # 'popfreq_all_20150413': 'f',
    # 'clinvar_20161128': 'f',
    # 'cosmic70': 'f',
    # 'nci60': 'f'

    hg_38_databases = OrderedDict({'knownGene': 'g',
                                   # 'cytoBand': 'r',
                                   # 'genomicSuperDups': 'r',
                                   # 'esp6500siv2_all': 'f',
                                   '1000g2015aug': 'f'})

    # 'clinvar_20161128': 'f',
    # 'cosmic70': 'f',
    # 'nci60': 'f'})

    def __init__(self, input_dir, output_dir, annovar_path, single_vcf_path, annovar_output_basename,
                 design_file=None, genome_build_version=None, custom_annovar_dbs_to_use=None):

        self.HUMANDB_FOLDER_NAME = "/humandb/"

        # TODO: Figure out what this string actually *is*
        self.DOWN_DD = '-webfrom annovar'

        # TODO: Pull out string keys into symbolic constants
        self.ANNOVAR_DB_IS_HOSTED_BY_ANNOVAR = OrderedDict({'knownGene': True,
                                                            # 'tfbsConsSites': False,
                                                            # 'cytoBand': False,
                                                            # 'targetScanS': False,
                                                            # 'genomicSuperDups': False,
                                                            # 'esp6500siv2_all': True,
                                                            '1000g2015aug': True
                                                            # 'popfreq_all_20150413': True,
                                                            # 'clinvar_20161128': True,
                                                            # 'cosmic70': True,
                                                            # 'nci60': True,
                                                            # 'avdblist': True
                                                            })

        # Project data
        self.input_dir = input_dir
        self.annovar_path = annovar_path
        self.design_file = design_file
        self._single_vcf_path = single_vcf_path
        self._output_dir = output_dir
        self._annovar_output_basename = annovar_output_basename
        # self.vcf_mapping_dict = vcf_mapping_dict

        # Databases data
        self._set_annovar_dbs_to_use(genome_build_version, custom_annovar_dbs_to_use)

    @property
    def genome_build_version(self):
        return self._genome_build_version

    @genome_build_version.setter
    def genome_build_version(self, value):
        raise NotImplementedError("genome_build_version setter not implemented")

    @property
    def annovar_dbs_to_use(self):
        return self._annovar_dbs_to_use

    @annovar_dbs_to_use.setter
    def annovar_dbs_to_use(self, value):
        raise NotImplementedError("annovar_dbs_to_use setter not implemented")

    def download_dbs(self):
        """
        Wrapper around annotate_variation.pl with -downdb as optional arg
        First, it cleans up the humandb/ directory to avoid conflicts, then gets newest versions of databases
        by spawning the jobs using subprocesses
        """
        list_commands = self._build_annovar_database_download_command_str()
        for command in list_commands:
            args = shlex.split(command)
            subprocess.call(args)

            logging.info('Finished downloading databases to {}'.format(
                os.path.join(self.annovar_path, self.HUMANDB_FOLDER_NAME)))

    def run_annovar(self, vcf_is_multisample=False):
        """ Spawn ANNOVAR VCF annotation jobs in batches of five/ten? files at a time to prevent memory overflow """

        annovar_txt_output_fp = os.path.join(self._output_dir, self._annovar_output_basename)
        cmd_string = self._build_table_annovar_command_str(self._single_vcf_path, annovar_txt_output_fp,
                                                           vcf_is_multisample=vcf_is_multisample)
        args = shlex.split(cmd_string)
        logging.info('Running Annovar')
        subprocess.call(args)
        logging.info('Finished running Annovar')

        return annovar_txt_output_fp

    def _build_table_annovar_command_str(self, vcf_path, csv_path, vcf_is_multisample=False):
        """Generate command string to run table_annovar.pl, which annotates a VCF file."""
        dbs = ",".join(list(self.annovar_dbs_to_use.keys()))
        dbs_args = ",".join(list(self.annovar_dbs_to_use.values()))
        if '1000g2015aug' in dbs:
            dbs = dbs.replace('1000g2015aug', '1000g2015aug_all')
        command = " ".join([
            'perl', os.path.join(self.annovar_path, 'table_annovar.pl'),
            vcf_path, "".join([self.annovar_path, self.HUMANDB_FOLDER_NAME]),
            '--buildver', self.genome_build_version,
            '-out', csv_path, '-remove -protocol', dbs,
            '-operation', dbs_args, '-nastring .', '-otherinfo -vcfinput'
        ])
        return command

    def _build_annovar_database_download_command_str(self):
        """Concatenate command string arguments for Annovar download database jobs"""

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

    def _set_annovar_dbs_to_use(self, genome_build_version, custom_annovar_dbs_to_use=None):
        self._genome_build_version = genome_build_version
        annovar_dbs_for_build_version_dict = self._get_annovar_dbs_to_use_for_build_version()

        annovar_dbs_to_get_dict = annovar_dbs_for_build_version_dict  # by default, assume we use all for build version
        if custom_annovar_dbs_to_use is not None:
            for annovar_db_name in custom_annovar_dbs_to_use:
                if annovar_db_name not in annovar_dbs_for_build_version_dict:
                    raise ValueError('Database %s not supported for build version %s' % (annovar_db_name,
                                                                                         self.genome_build_version))
                annovar_dbs_to_get_dict = {db: annovar_dbs_for_build_version_dict[db] for db in
                                           custom_annovar_dbs_to_use}

        self._annovar_dbs_to_use = annovar_dbs_to_get_dict

    def _get_annovar_dbs_to_use_for_build_version(self):
        if self.genome_build_version == 'hg19':
            databases = AnnovarWrapper.hg_19_databases
        elif self.genome_build_version == 'hg38':
            databases = AnnovarWrapper.hg_38_databases
        else:
            raise NameError(
                'Genome {0} not supported by VAPr, genome must be one of the following: hg19, or hg38'.format(
                    self.genome_build_version))

        return databases
