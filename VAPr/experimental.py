import subprocess
import shlex
from watchdog.events import FileSystemEventHandler
from watchdog.observers import Observer


class AnnovarWrapper(object):

    """Run Annovar as subprocess"""

    def __init__(self, input_vcf, output_csv, annovar_path):

        self.input = input_vcf
        self.output = output_csv + 'annotated'
        self.output_csv_path = output_csv
        self.path = annovar_path
        self.supported_databases = ['knownGene',
                                    'tfbsConsSites',
                                    'cytoBand',
                                    'targetScanS',
                                    'genomicSuperDups',
                                    'esp6500siv2_all',
                                    '1000g2015aug_all',
                                    'snp138',
                                    'popfreq_all_20150413',
                                    'clinvar_20161128',
                                    'cosmic70',
                                    'nci60']

        self.annovar_command_str = self.build_arg_str()

    def build_arg_str(self):
        dbs = ",".join(self.supported_databases)
        command = " ".join(['sudo', 'perl', self.path, 'table_annovar.pl',
                            self.input, self.path, 'humandb/', '-buildver', 'hg19',
                            '-out', self.output, '-remove -protocol', dbs,  '-operation',
                            'g,r,r,r,r,f,f,f,f,f,f,f', '-nastring .', '-csvout'])

        return command

    def run_annovar(self):

        args = shlex.split(self.annovar_command_str)
        p = subprocess.Popen(args, stdout=subprocess.PIPE)
        #run_handler(self.output_csv_path)

        return 'Finished running ANNOVAR on {}'.format(self.input)



class MyHandler(FileSystemEventHandler):
    """
    Overwrite the methods for creation of files as annovar runs. Once the .csv file is detected, exit process
    and proceed to next file.
    """

    def __init__(self, observer):
        self.observer = observer

    def on_created(self, event):
        if event.src_path.split('/')[-1].split('.')[-1] != "invalid_input":
            if event.src_path.split('/')[-1].split('.')[-1] != "log":
                if event.src_path.split('/')[-1].split('.')[-1] != "log":
                    print("Currently working on VCF file: " + event.src_path.split('/')[-1].split('.')[0] + ", field " +\
                          event.src_path.split('/')[-1].split('.')[-1])
        if event.src_path[-3:] == 'csv':
            print("\nAnnovar finished working on file : " + event.src_path.split('/')[-1].split('.')[0] + \
                  " has finished. A .csv file has been created in the OUT_PATH directory\n")
            self.observer.stop()


def run_handler(output_csv_path):
    observer = Observer()
    event_handler = MyHandler(observer)
    observer.schedule(event_handler, output_csv_path, recursive=True)
    observer.start()
    observer.join()

    return None
