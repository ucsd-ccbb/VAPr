import annotation_project

if __name__ == '__main__':
    # TODO: Figure out why this script is very similar to "run_annotations.py" script, get rid of one or other or both

    # TODO: All paths below must be changed to something else and contents put into github repository if
    # this module stays in project

    # Directory of input files to be annotated
    IN_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_benchmark/"

    # Output file directory
    OUT_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_benchmark/"
    # Location of your annovar dowload. The folder should contain the following files/directories:
    ANNOVAR_PATH = '/Volumes/Carlo_HD1/CCBB/annovar/'

    db_name = 'VariantDBMultiBenchmark'
    collection_name = 'test_collection_small_vcf'

    # Design file (optional)
    # design_file = '/Volumes/Carlo_HD1/CCBB/VAPr_files/guorong_single_sample.csv'

    # Database and collection names (optional)
    # TODO: I believe the comment above (that proj_data containing db_name and collection_name is optional) is incorrect
    proj_data = annotation_project.make_mongo_db_and_collection_names_dict(db_name, collection_name)

    Project = annotation_project.AnnotationProject(IN_PATH,
                                                   OUT_PATH,
                                                   ANNOVAR_PATH,
                                                   proj_data,
                                                   # design_file=design_file,
                                                   build_ver='hg19')

    Project.parallel_annotation(n_processes=8, verbose=2)
