from vapr import base

if __name__ == '__main__':
    # Directory of input files to be annotated
    IN_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_benchmark/"

    # Output file directory
    OUT_PATH = "/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_benchmark/"
    # Location of your annovar dowload. The folder should contain the following files/directories:
    ANNOVAR_PATH = '/Volumes/Carlo_HD1/CCBB/annovar/'

    db_name = 'VariantDBMultiBenchmark'
    collection_name = 'test_collection_small_vcf'

    # Design File (optional)
    # design_file = '/Volumes/Carlo_HD1/CCBB/VAPr_files/guorong_single_sample.csv'

    # Databse and Collection names (optional)
    proj_data = {'db_name': db_name,
                 'project_name': collection_name}

    Project = base.AnnotationProject(IN_PATH,
                                OUT_PATH,
                                ANNOVAR_PATH,
                                proj_data,
                                # design_file=design_file,
                                build_ver='hg19')

    Project.parallel_annotation_and_saving(n_processes=8, verbose=2)