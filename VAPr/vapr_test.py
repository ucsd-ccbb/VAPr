from VAPr.annotation_project import AnnotationProject
IN_PATH = "/Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_input_dir"
OUT_PATH = "/Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/test_out_csv_path/merged"
ANALYSIS_NAME = "merged_test"
ANNOVAR_PATH = '/Users/adammark/software/VAPr_code_review/VAPr/tests/test_files/annovar_dir'
proj_data = {'db_name': 'VariantDBBenchmarkWES', 'collection_name': 'collect'}
Project = AnnotationProject(IN_PATH, OUT_PATH,ANALYSIS_NAME,ANNOVAR_PATH,proj_data,build_ver='hg19')
#Project.run_annovar(multisample=True)
#Project.parallel_annotation_and_saving(n_processes=8)
Project.quick_annotate_and_save(n_processes=1)