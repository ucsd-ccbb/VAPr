import subprocess
import shlex
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer
from collections import OrderedDict

class AnnovarWrapper(object):

    """Run Annovar as subprocess"""

    def __init__(self, input_vcf, output_csv, annovar_path):

        self.input = input_vcf
        self.output = output_csv + 'annotated'
        self.output_csv_path = output_csv
        self.path = annovar_path
        self.down_dd = '-webfrom annovar'
        self.annovar_hosted = ['knownGene', 'esp6500siv2_all', '1000g2015aug', 'nci60',
                               'clinvar_20161128', 'popfreq_all_20150413', 'cosmic70']

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
        dbs_args = ",".join(list(self.supported_databases.values()))
        command = " ".join(['perl', self.path + 'table_annovar.pl', self.input,
                            self.path + 'humandb/', '-buildver', 'hg19', '-out',
                            self.output, '-remove -protocol', dbs,  '-operation',
                            dbs_args, '-nastring .', '-otherinfo -vcfinput'])

        return command

    def run_annovar(self):

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
            if db in self.annovar_hosted:
                command_list.append(" ".join(['perl', self.path + 'annotate_variation.pl', '-buildver', 'hg19', '-downdb',
                                    self.down_dd, db, self.path + 'humandb/']))
            else:
                command_list.append(" ".join(['perl', self.path + 'annotate_variation.pl', '-buildver', 'hg19', '-downdb',
                                              db, self.path + 'humandb/']))
        return command_list

    def download_dbs(self, all=True, dbs=None):

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
            if event.src_path.split('.')[-1] not in ["annovar_date"]:
                print("Currently downloading database file: " + event.src_path.split('/')[-1].split('.')[0])

                final = event.src_path[-3:]
                if final == 'txt':
                    print("\nAnnovar finished working on file : " + event.src_path.split('/')[-1].split('.')[0] + \
                          ". A .txt file has been created in the OUT_PATH directory\n")
                    self.observer.stop()

        else:
            if event.src_path.split('.')[-1] not in ["invalid_input", "log"]:
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
