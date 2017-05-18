import os
import json
import logging
import sys

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'

# ------ # ------ # ------ # LOCAL FILES CONFIG # ------ # ------ # ------ #
ROOT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), os.pardir, os.pardir) # This is your Project Root
DATA_DIR = os.path.join(ROOT_DIR, 'data')
CONFIG_PATH = os.path.join(ROOT_DIR, 'config/config')
TOOL_SPECIFIC_CONFIG = os.path.join(ROOT_DIR, 'config/tool_specific_config')
MODES_PATH = os.path.join(DATA_DIR, 'modes')
AWS_DIR = os.path.join(ROOT_DIR, 'awsCluster')
TEST_DIR = os.path.join(AWS_DIR, 'tests')
ANALISYS_STEPS = ["fastqc", "bwa-alignment", "post-alignment", "gatk-haplotype"]

with open(MODES_PATH) as f:
    MODES = json.load(f)

# MODELS_DIR = os.path.join(ROOT_DIR, 'models')
# RESUTLS_DIR = os.path.join(ROOT_DIR, 'reports')

# ------ # ------ # ------ # SERVER FILES CONFIG # ------ # ------ # ------ #

SERVER_WORKSPACE = '/shared/workspace/'
SOFTWARE_DIR = os.path.join(SERVER_WORKSPACE, 'software')
SERVER_DATA_DIR = os.path.join(SERVER_WORKSPACE, 'data_archive/')

# ------ # ------ # ------ # SOFTWARE TOOLS PATHS # ------ # ------ # ------ #

KALLISTO = os.path.join(SOFTWARE_DIR, 'kallisto_linux-v0.42.1/kallisto')
KALLISTO_INDEX = os.path.join(SOFTWARE_DIR, 'kallisto_linux-v0.42.1/kallisto_index/kallisto_index')
BT2_HOME = os.path.join(SOFTWARE_DIR, 'bowtie2-2.3.0-legacy')
ADAPTER = 'TGGAATTCTCGGGTGCCAAGG'
INDEX = 'Homo sapiens'
FASTQC_PATH = os.path.join(SOFTWARE_DIR, 'FastQC/fastqc')
TRIMMO_PATH = os.path.join(SOFTWARE_DIR, 'Trimmomatic-0.36/trimmomatic-0.36.jar')
CUTADAPT_PATH = os.path.join(SOFTWARE_DIR, 'cutadapt-1.12/bin/cutadapt')
ALL_GENOMES_FA_FILEPATH = os.path.join(SERVER_WORKSPACE, 'data/hairpin.fa')  # TODO: add this fie, can't be found in server
INDICES_DIR = os.path.join(SERVER_WORKSPACE, 'data/indices')


logging.basicConfig(stream=sys.stdout, level=logging.INFO)


if __name__ == '__main__':
    import yaml

    with open(TOOL_SPECIFIC_CONFIG + "/BWA.yaml", 'r') as stream:
        try:
            print(yaml.load(stream))
        except yaml.YAMLError as exc:
            print(exc)