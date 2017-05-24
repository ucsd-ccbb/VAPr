import os
import sys
import shlex
import glob
import csv
import time
import datetime
import subprocess
from collections import OrderedDict
import logging
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from base import AnnotationProject
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout     # Enables logging on jupyter notebooks
except:
    pass


class AnnovarWrapper(AnnotationProject):

    """ Wrapper around Annovar download and annotation functions """

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
            run_handler(os.path.join(self.annovar, 'humandb/'), cmds=list_commands, annovar_path=self.annovar)

        return 'Finished downloading databases to {}'.format(os.path.join(self.annovar, 'humandb/'))

    def run_annovar(self):
        """ Spawning Annovar jobs """

        print(self.mapping)
        for sample in self.mapping.keys():
            annotation_dir = os.path.join(self.output_csv_path, sample)
            if os.path.isdir(annotation_dir):
                raise TypeError('Directory already exists for %s. Aborting.' % annotation_dir)
            else:
                os.makedirs(annotation_dir)

            for vcf, csv in self.mapping[sample]['vcf_csv']:
                vcf_path = os.path.join(self.input_dir, os.path.join(sample, vcf))
                csv_path = os.path.join(self.output_csv_path, os.path.join(sample, csv))
                cmd_string = self.build_annovar_command_str(vcf_path, csv_path)
                args = shlex.split(cmd_string)

                subprocess.Popen(args, stdout=subprocess.PIPE)

            n_commands = len(self.mapping[sample]['vcf_csv'])
            logging.info('Annovar jobs submitted for sample %s' % sample)
            listen(os.path.join(self.output_csv_path, sample), n_commands)
            logging.info('Finished running Annovar on sample %s' % sample)

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
                        logging.info('Database %s outdated, will download newer version' % db_)
                        self.download_dbs(all_dbs=False, dbs=[os.path.splitext(os.path.splitext(db_)[0])[0]])

    def get_databases(self):

        if self.buildver == 'hg18':
            databases = self.hg_18_databases
        elif self.buildver == 'hg19':
            databases = self.hg_19_databases
        else:
            databases = self.hg_38_databases

        return databases


def listen(out_path, n_files):

    added = 0

    while True:

        files = os.listdir(out_path)
        txts = sorted([i for i in files if i.endswith('txt')])

        time.sleep(5)

        new_files = os.listdir(out_path)
        new_txts = sorted([i for i in new_files if i.endswith('txt')])

        newly_created = [x for x in new_txts if x not in txts]

        if len(newly_created) > 0:
            added += 1
            logging.info('File %i/%i: Annovar finished working on file : ' % (added, n_files) +
                         os.path.basename(newly_created[0]) +
                         '.\n A text file has been created in the %s directory\n' % out_path)

        if added == n_files:
            break
        if len(txts) >= n_files:
            break


class MyHandler(FileSystemEventHandler):
    """
    Overwrite the methods for creation of files as annovar runs. Once the .csv file is detected, exit process
    and proceed to next file.
    """

    def __init__(self, observer, cmds=None, annovar_path='ANNOVAR_PATH'):

        self.observer = observer
        self.annovar = annovar_path
        self.cmds = cmds

    def on_created(self, event):

        logging.info('Downloading: ' + event.src_path)

        if self.cmds[-1][-2] in event.src_path:
            self.observer.stop()


def run_handler(output_csv_path, cmds=None, annovar_path='ANNOVAR_PATH'):
    observer = Observer()
    event_handler = MyHandler(observer,
                              cmds=cmds,
                              annovar_path=annovar_path)

    observer.schedule(event_handler, output_csv_path, recursive=True)
    observer.start()
    observer.join()

    return None
