from __future__ import division, print_function
import subprocess
import shlex
import os
import csv
import glob
import datetime
import warnings
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from collections import OrderedDict


class AnnovarWrapper(object):

    """Run Annovar as subprocess"""

    def __init__(self, input_vcf, output_csv, annovar_path):

        self.input = input_vcf
        self.output = output_csv + os.path.splitext(os.path.basename(self.input))[0] + '_annotated'
        self.output_csv_path = output_csv
        self.path = annovar_path
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

        self.supported_databases = OrderedDict({'knownGene': 'g',
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

        self.annovar_command_str = self.build_annovar_command_str()

    def build_annovar_command_str(self):
        dbs = ",".join(list(self.supported_databases.keys()))
        dbs = dbs.replace('1000g2015aug', '1000g2015aug_all')
        dbs_args = ",".join(list(self.supported_databases.values()))
        command = " ".join(['perl', self.path + 'table_annovar.pl', self.input,
                            self.path + 'humandb/', '-buildver', 'hg19', '-out',
                            self.output, '-remove -protocol', dbs,  '-operation',
                            dbs_args, '-nastring .', '-otherinfo -vcfinput'])

        return command

    def run_annovar(self):

        files = glob.glob(self.output_csv_path + '*')
        for f in files:
            if os.path.basename(f).startswith(os.path.basename(self.input).split('.')[0]):
                os.remove(f)

        args = shlex.split(self.annovar_command_str)
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        run_handler(self.output_csv_path)

        return 'Finished running ANNOVAR on {}'.format(self.input)

    def build_db_dl_command_str(self, all, dbs):

        if not all:
            databases = dbs
        else:
            databases = self.supported_databases.keys()

        command_list = []

        for db in databases:
            if self.annovar_hosted[db]:
                command_list.append(" ".join(['perl', self.path + 'annotate_variation.pl', '-buildver', 'hg19', '-downdb',
                                    self.down_dd, db, self.path + 'humandb/']))
            else:
                command_list.append(" ".join(['perl', self.path + 'annotate_variation.pl', '-buildver', 'hg19', '-downdb',
                                              db, self.path + 'humandb/']))
        return command_list

    def check_for_database_updates(self):
        self.download_dbs(all=False, dbs=['avdblist'])

        with open(self.path + '/humandb/hg19_avdblist.txt', 'r') as db_list:
            reader = csv.reader(db_list, delimiter='\t')
            db_dict = {}
            for i in reader:

                if i[0][-6:] != 'idx.gz':
                    db_dict[i[0][5:]] = [datetime.datetime(int(i[1][0:4]), int(i[1][4:6]), int(i[1][6:8])), i[2]]

        for db in self.manual_update.keys():
            for db_ in db_dict.keys():
                if db_.startswith(db):
                    if db_dict[db_][0] > self.manual_update[db][0]:
                        print('Database %s outdated, will download newer version' % db_)
                        self.download_dbs(all=False, dbs=[os.path.splitext(os.path.splitext(db_)[0])[0]])

    def download_dbs(self, all=True, dbs=None):
        if len(os.listdir(self.path + 'humandb/')) > 0:
            files = glob.glob(self.path + 'humandb/*')
            for f in files:
                os.remove(f)

        list_commands = self.build_db_dl_command_str(all, dbs)
        for command in list_commands:
            args = shlex.split(command)

            p = subprocess.Popen(args, stdout=subprocess.PIPE)
            run_handler(self.path + 'humandb/', download=True)

        return 'Finished downloading databases to {}'.format(self.path + 'humandb/')


class MyHandler(FileSystemEventHandler):
    """
    Overwrite the methods for creation of files as annovar runs. Once the .csv file is detected, exit process
    and proceed to next file.
    """

    def __init__(self, observer, download=False):

        self.observer = observer
        self.dl = download

    def on_created(self, event):

        if self.dl:
            dowloaded = ["annovar_date"]
            if event.src_path.split('.')[-1] not in dowloaded:
                print("Currently downloading database file: " + event.src_path.split('/')[-1].split('.')[0])
                dowloaded.append(event.src_path.split('/')[-1].split('.')[0])

            final = event.src_path[-3:]
            if final == 'txt':
                print("\nAnnovar finished dowloading on file : " + event.src_path.split('/')[-1].split('.')[0] + \
                      ". A .txt file has been created in the ANNOVAR_PATH directory\n")
                self.observer.stop()

        else:
            if event.src_path.split('.')[-1] not in ["invalid_input", "log", "avinput"]:
                print("Currently working on VCF file: " + event.src_path.split('/')[-1].split('.')[0] + ", field " +\
                      event.src_path.split('/')[-1].split('.')[-1])

                final = event.src_path[-3:]
                if final == 'txt':
                    print("\nAnnovar finished working on file : " + event.src_path.split('/')[-1].split('.')[0] +\
                          ". A .txt file has been created in the OUT_PATH directory\n")
                    self.observer.stop()


def run_handler(output_csv_path, download=False):
    observer = Observer()
    event_handler = MyHandler(observer, download=download)
    observer.schedule(event_handler, output_csv_path, recursive=True)
    observer.start()
    observer.join()

    return None
