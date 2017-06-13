import logging
import sys
from collections import OrderedDict
import datetime

__author__ = 'Carlo Mazzaferro<cmazzafe@ucsd.edu>'

logging.basicConfig(stream=sys.stdout, level=logging.INFO)
chunk_size = 5000


# ------ #  ANNOVAR CONFIG # ------ #

down_dd = '-webfrom annovar'
annovar_hosted = OrderedDict({'knownGene': True,
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

dl_list_command = 'avdblist'
manual_update = {'clinvar_20161128': [datetime.datetime(2016, 11, 28)],
                 '1000g2015aug': [datetime.datetime(2016, 8, 30)],
                 'popfreq_all_20150413': [datetime.datetime(2015, 4, 13)]}

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