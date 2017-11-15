import logging
import os
import shlex
import subprocess
from collections import OrderedDict


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

    hg_38_databases = OrderedDict({'knownGene': 'g',
                                   # 'cytoBand': 'r',
                                   # 'genomicSuperDups': 'r',
                                   # 'esp6500siv2_all': 'f',
                                   '1000g2015aug': 'f'})

    @classmethod
    def _get_annovar_dbs_to_use_for_build_version(cls, genome_build_version):
        # TODO: I'm not happy about having these strings here; they should be symbolic constants
        if genome_build_version == 'hg19':
            databases = cls.hg_19_databases
        elif genome_build_version == 'hg38':
            databases = cls.hg_38_databases
        else:
            raise ValueError("Genome build version '{0}' not supported.".format(genome_build_version))

        return databases

    def __init__(self, annovar_install_path, genome_build_version, custom_annovar_dbs_to_use=None):
        self._HUMANDB_FOLDER_NAME = "/humandb/"
        self._DOWN_DD = '-webfrom annovar'

        # TODO: Pull out string keys into symbolic constants
        self._ANNOVAR_DB_IS_HOSTED_BY_ANNOVAR = OrderedDict({'knownGene': True,
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

        self._annovar_install_path = annovar_install_path
        self._genome_build_version = genome_build_version
        self._annovar_dbs_to_use = self._get_annovar_dbs_to_use(custom_annovar_dbs_to_use)

    def download_databases(self):
        list_commands = self._build_annovar_database_download_command_str()
        for command in list_commands:
            args = shlex.split(command)
            subprocess.call(args)

        logging.info('Finished downloading databases to {0}'.format(
            os.path.join(self._annovar_install_path, self._HUMANDB_FOLDER_NAME)))

    def run_annotation(self, single_vcf_path, output_basename, output_dir):
        annovar_output_basename = output_basename + '_annotated'
        annovar_output_base = os.path.join(output_dir, annovar_output_basename)
        annovar_txt_output_fp = annovar_output_base + '.' + self._genome_build_version + '_multianno.txt'

        cmd_string = self._build_table_annovar_command_str(single_vcf_path, annovar_output_base)
        args = shlex.split(cmd_string)
        logging.info('Running Annovar')
        subprocess.call(args)
        logging.info('Finished running Annovar')
        return annovar_txt_output_fp

    def _get_annovar_dbs_to_use(self, custom_annovar_dbs_to_use=None):
        annovar_dbs_for_build_version_dict = self._get_annovar_dbs_to_use_for_build_version(
            self._genome_build_version)

        annovar_dbs_to_get_dict = annovar_dbs_for_build_version_dict  # by default, assume use all for build version
        if custom_annovar_dbs_to_use is not None:
            for annovar_db_name in custom_annovar_dbs_to_use:
                if annovar_db_name not in annovar_dbs_for_build_version_dict:
                    raise ValueError('Database %s not supported for build version %s' % (annovar_db_name,
                                                                                         self._genome_build_version))
                annovar_dbs_to_get_dict = {db: annovar_dbs_for_build_version_dict[db] for db in
                                           custom_annovar_dbs_to_use}

        return annovar_dbs_to_get_dict

    def _build_table_annovar_command_str(self, vcf_path, csv_path):
        """Generate command string to run table_annovar.pl, which annotates a VCF file."""

        dbs = ",".join(list(self._annovar_dbs_to_use.keys()))
        dbs_args = ",".join(list(self._annovar_dbs_to_use.values()))
        if '1000g2015aug' in dbs:
            dbs = dbs.replace('1000g2015aug', '1000g2015aug_all')

        command = " ".join([
            'perl', os.path.join(self._annovar_install_path, 'table_annovar.pl'),
            vcf_path, "".join([self._annovar_install_path, self._HUMANDB_FOLDER_NAME]),
            '--buildver', self._genome_build_version,
            '-out', csv_path, '-remove -protocol', dbs,
            '-operation', dbs_args, '-nastring .', '-otherinfo -vcfinput'
        ])
        return command

    def _build_annovar_database_download_command_str(self):
        """Concatenate command string arguments for Annovar download database jobs"""

        command_list = []
        # TODO: someday: refactor to remove duplicated command components :(
        for annovar_db_name in self._annovar_dbs_to_use:
            if self._ANNOVAR_DB_IS_HOSTED_BY_ANNOVAR[annovar_db_name]:
                command_list.append(" ".join(
                    ['perl', self._annovar_install_path + '/annotate_variation.pl',
                     '-build', self._genome_build_version,
                     '-downdb', self._DOWN_DD, annovar_db_name,
                     self._annovar_install_path + self._HUMANDB_FOLDER_NAME]))
            else:
                command_list.append(" ".join(
                    ['perl', self._annovar_install_path + '/annotate_variation.pl',
                     '-build', self._genome_build_version,
                     '-downdb', annovar_db_name,
                     self._annovar_install_path + self._HUMANDB_FOLDER_NAME]))

        return command_list
