from __future__ import division, print_function
import subprocess
import shlex
import os
import csv
import glob
import pandas
import datetime
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from collections import OrderedDict
from VAPr import parser_models
import logging
logger = logging.getLogger(__name__)


class ProjectData(object):

    supported_build_vers = ['hg19', 'hg18', 'hg38']

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None):

        self.input_dir = input_dir
        self.output_csv_path = output_dir
        self.vcf_files = self.find_vcfs()
        self.design_file = design_file
        self.input_output_mapping = self.get_mapping()
        self.annovar = annovar_path
        self.project_data = project_data
        self.buildver = build_ver
        self.mapping = [(k, v) for k, v in self.get_mapping().items()]

    def get_mapping(self):

        if self.design_file:
            design_df = pandas.read_csv(self.design_file)
            samples = design_df.Sample.values
            sample_to_vcf_map = self.check_file_existance_return_mapping(samples)

        else:
            sample_to_vcf_map = {os.path.splitext(os.path.basename(vcf))[0]: vcf for vcf in self.vcf_files}

        return {os.path.join(self.input_dir, vcf): self.output_csv_path + os.path.splitext(os.path.basename(vcf))[0] +
                '_annotated' for vcf in sample_to_vcf_map.values()}

    def find_vcfs(self):
        vcf_files = os.listdir(self.input_dir)
        return [vcf for vcf in vcf_files if vcf.endswith('vcf')]

    def check_file_existance_return_mapping(self, sample_names):

        mappings = {}
        combined = '\t'.join(self.vcf_files)

        for sample in sample_names:
            if sample not in combined:
                raise FileNotFoundError('Sample %s does not have an associated vcf file' % sample)

            for vcf in self.vcf_files:
                if sample in vcf:
                    mappings[sample] = vcf

        return mappings

    def check_ver(self, build_ver):

        if not build_ver:
            self.buildver = 'hg19'  # Default genome build vesion

        if build_ver not in self.supported_build_vers:
            raise ValueError('Build version must not recognized. Supported builds are'
                             ' %s, %s, %s' % (self.supported_build_vers[0],
                                              self.supported_build_vers[1],
                                              self.supported_build_vers[2]))
        else:
            self.buildver = build_ver

        return self.buildver


class AnnotationProject(ProjectData):

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None):

        super().__init__(input_dir,
                         output_dir,
                         annovar_path,
                         project_data,
                         design_file=design_file,
                         build_ver=build_ver)

        self.annovar_wrapper = AnnovarWrapper(self.input_dir,
                                              self.output_csv_path,
                                              self.annovar,
                                              self.project_data,
                                              design_file=self.design_file,
                                              build_ver=self.buildver)

        self.annotator_wrapper = parser_models.VariantParsing(self.input_dir,
                                                              self.output_csv_path,
                                                              self.annovar,
                                                              self.project_data,
                                                              design_file=self.design_file,
                                                              build_ver=self.buildver)

    def download_dbs(self, all_dbs=True, dbs=None):
        self.annovar_wrapper.download_dbs(all_dbs=all_dbs, dbs=dbs)

    def run_annovar(self):
        self.annovar_wrapper.run_annovar()

    def annotate_and_save(self):
        self.annotator_wrapper.annotate_and_save(buffer=False)

    def parallel_annotation_and_saving(self):
        self.annotator_wrapper.parallel_annotation(n_processes=4)


class AnnovarWrapper(AnnotationProject):

    """Run Annovar as subprocess"""

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

        self.down_dd = '-webfrom annovar'
        self.annovar_hosted = OrderedDict({'knownGene': True,
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

        self.dl_list_command = 'avdblist'
        self.manual_update = {'clinvar_20161128': [datetime.datetime(2016, 11, 28)],
                              '1000g2015aug':  [datetime.datetime(2016, 8, 30)],
                              'popfreq_all_20150413': [datetime.datetime(2015, 4, 13)]
                              }

        self.hg_18_databases = OrderedDict({'knownGene': 'g',
                                            'tfbsConsSites': 'r',
                                            'cytoBand': 'r',
                                            'targetScanS': 'r',
                                            'genomicSuperDups': 'r',
                                            'esp6500siv2_all': 'f',
                                            '1000g2015aug': 'f',
                                            'cosmic70': 'f',
                                            'nci60': 'f',
                                            })

        self.hg_19_databases = OrderedDict({'knownGene': 'g',
                                            'tfbsConsSites': 'r',
                                            'cytoBand': 'r',
                                            'targetScanS': 'r',
                                            'genomicSuperDups': 'r',
                                            'esp6500siv2_all': 'f',
                                            '1000g2015aug': 'f',
                                            'popfreq_all_20150413': 'f',
                                            'clinvar_20161128': 'f',
                                            'cosmic70': 'f',
                                            'nci60': 'f',
                                            })

        self.hg_38_databases = OrderedDict({'knownGene': 'g',
                                            'cytoBand': 'r',
                                            'genomicSuperDups': 'r',
                                            'esp6500siv2_all': 'f',
                                            '1000g2015aug': 'f',
                                            'clinvar_20161128': 'f',
                                            'cosmic70': 'f',
                                            'nci60': 'f',
                                            })

        self.databases = self.get_databases()

    def download_dbs(self, all_dbs=True, dbs=None):

        if len(os.listdir(os.path.join(self.annovar, 'humandb/'))) > 0:
            files = glob.glob(os.path.join(self.annovar, 'humandb/*'))
            for f in files:
                os.remove(f)

        list_commands = self.build_db_dl_command_str(all_dbs, dbs)
        for command in list_commands:
            args = shlex.split(command)
            subprocess.Popen(args, stdout=subprocess.PIPE)
            run_handler(os.path.join(self.annovar, 'humandb/'), cmds=list_commands, download=True,
                        annovar_path=self.annovar)

        return 'Finished downloading databases to {}'.format(os.path.join(self.annovar, 'humandb/'))

    def run_annovar(self):
        """
        files = glob.glob(self.output_csv_path + '*')
        for f in files:
            if os.path.basename(f).startswith(os.path.basename(self.input).split('.')[0]):
                os.remove(f)
        """
        for vcf, csv in self.input_output_mapping.items():
            cmd_string = self.build_annovar_command_str(vcf, csv)
            args = shlex.split(cmd_string)
            subprocess.Popen(args, stdout=subprocess.PIPE)

        logging.INFO('Annovar jobs submitted')
        run_handler(self.output_csv_path, annovar_path=self.annovar, num_commands=len(self.input_output_mapping.keys()))

        return 'Finished running ANNOVAR on {}'.format(self.input_dir)

    def build_annovar_command_str(self, vcf, csv):

        dbs = ",".join(list(self.databases.keys()))
        dbs_args = ",".join(list(self.databases.values()))

        if '1000g2015aug' in dbs:
            dbs = dbs.replace('1000g2015aug', '1000g2015aug_all')
        command = " ".join(['perl', os.path.join(self.annovar, 'table_annovar.pl'), vcf,
                            os.path.join(self.annovar, 'humandb/'), '-buildver', self.buildver, '-out',
                            csv, '-remove -protocol', dbs,  '-operation',
                            dbs_args, '-nastring .', '-otherinfo -vcfinput'])

        return command

    def build_db_dl_command_str(self, all_dbs, dbs):

        if not all_dbs:
            for db in dbs:
                if db not in self.databases:
                    raise ValueError('Database %s not supported for buil version %s' % (db, self.buildver))
            self.databases = {db: self.databases[db] for db in dbs}

        command_list = []

        for db in self.databases:

            if self.annovar_hosted[db]:
                command_list.append(" ".join(['perl', self.annovar + 'annotate_variation.pl', '-build', self.buildver,
                                              '-downdb', self.down_dd, db, self.annovar + 'humandb/']))
            else:
                command_list.append(" ".join(['perl', self.annovar + 'annotate_variation.pl', '-build', self.buildver,
                                              '-downdb', db, self.annovar + 'humandb/']))
        return command_list

    def check_for_database_updates(self):

        self.download_dbs(all_dbs=False, dbs=['avdblist'])

        with open(os.path.join(self.annovar, '/humandb' + self.buildver + '_avdblist.txt'), 'r') as db_list:
            reader = csv.reader(db_list, delimiter='\t')
            db_dict = {}
            for i in reader:

                if i[0][-6:] != 'idx.gz':
                    db_dict[i[0][5:]] = [datetime.datetime(int(i[1][0:4]), int(i[1][4:6]), int(i[1][6:8])), i[2]]

        for db in self.manual_update.keys():
            for db_ in db_dict.keys():
                if db_.startswith(db):
                    if db_dict[db_][0] > self.manual_update[db][0]:
                        logging.INFO('Database %s outdated, will download newer version' % db_)
                        self.download_dbs(all_dbs=False, dbs=[os.path.splitext(os.path.splitext(db_)[0])[0]])

    def get_databases(self):

        if self.buildver == 'hg18':
            databases = self.hg_18_databases
        elif self.buildver == 'hg19':
            databases = self.hg_19_databases
        else:
            databases = self.hg_38_databases

        return databases


class MyHandler(FileSystemEventHandler):
    """
    Overwrite the methods for creation of files as annovar runs. Once the .csv file is detected, exit process
    and proceed to next file.
    """

    def __init__(self, observer, cmds=None, download=False, annovar_path='ANNOVAR_PATH', num_commands=0):

        self.observer = observer
        self.dl = download
        self.annovar = annovar_path
        self.cmds = cmds
        self.num_commands = num_commands

    def on_created(self, event):

        files_created = 0
        if self.dl:

            logging.INFO('Downloading: ' + event.src_path)

            if self.cmds[-1][-2] in event.src_path:
                self.observer.stop()

        else:
            """

            if event.src_path.split('.')[-1] not in ["invalid_input", "log", "avinput"]:
                logging.INFO("Currently working on VCF file: " + event.src_path.split('/')[-1].split('.')[0] + ", field " +
                      event.src_path.split('/')[-1].split('.')[-1])
            """
            final = event.src_path[-3:]
            if final == 'txt':
                files_created += 1
                logging.INFO('\nAnnovar finished working on file : ' + event.src_path.split('/')[-1].split('.')[0] +
                             '. A .txt file has been created in the OUT_PATH directory\n')

                if files_created == self.num_commands:
                    self.observer.stop()


def run_handler(output_csv_path, cmds=None, download=False, annovar_path='ANNOVAR_PATH', num_commands=0):
    observer = Observer()
    event_handler = MyHandler(observer,
                              cmds=cmds,
                              download=download,
                              annovar_path=annovar_path,
                              num_commands=num_commands)

    observer.schedule(event_handler, output_csv_path, recursive=True)
    observer.start()
    observer.join()

    return None
