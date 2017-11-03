import time
x = time.time()
from VAPr.annotation_project import AnnotationProject
IN_PATH = "/Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_input_dir"
OUT_PATH = "/Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_out_csv_path/merged"
ANALYSIS_NAME = "merged_test"
ANNOVAR_PATH = '/Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir'
proj_data = {'db_name': 'VariantDBBenchmarkWES', 'collection_name': 'collect'}
Project = AnnotationProject(IN_PATH, OUT_PATH,ANALYSIS_NAME,ANNOVAR_PATH,proj_data,build_ver='hg19')
#Project.download_dbs()
# Project.run_annovar(multisample=True)
Project.parallel_annotation_and_saving(n_processes=8)
#Project.quick_annotate_and_save(n_processes=4)
y = time.time()
print((y-x)/60)

db_name="VariantDBBenchmarkWES"
collection_name="collect"
from pymongo import MongoClient
from VAPr.queries import Filters
client = MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)
db = getattr(client, db_name)
collection = getattr(db, collection_name)
filter_collection = Filters(db_name, collection_name)
rare_cancer_variants = filter_collection.rare_cancer_variant()


# import time
# x = time.time()
# rare_disease_variants = filter_collection.rare_disease_variants()
# y = time.time()
# (y-x)/60

# def process_format_string(field, anno, ):
#     processed_dict = {}
#     try:
#         processed_dict[field] = float(anno)
#     except:
#         try:
#             processed_dict[field] = anno.split(",")
#         except:
#             processed_dict[field] = anno
#     return processed_dict

# def process_info_string(info_string):
#     anno_list = info_string.split(";")
#     processed_dict = {}
#     for anno in anno_list:
#         split = anno.split("=")
#         k = split[0]
#         v = split[1]
#         try:
#             processed_dict[k] = float(v)
#         except:
#             processed_dict[k] = v
#     return processed_dict