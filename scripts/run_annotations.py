import annotation_project

if __name__ == '__main__':
    # TODO: Figure out why this script is very similar to "test.py" script, get rid of one or other or both

    # TODO: All paths below must be changed to something else and contents put into github repository if
    # this module stays in project

    # Directory of input files to be annotated
    IN_PATH = "/mnt/data1/carlom/samples/"

    # Output file directory
    OUT_PATH = "/mnt/data1/carlom/csv_results/"
    # Location of your annovar dowload. The folder should contain the following files/directories:
    ANNOVAR_PATH = '/mnt/data1/carlom/annovar/'

    # Design file (optional)
    design_file = '/mnt/data1/carlom/design_file_three_samples.csv'

    # Database and collection names (optional)
    # TODO: I believe the comment above (that proj_data containing db_name and collection_name is optional) is incorrect
    # TODO: Apparent bug:
    # This code sets "project_name", not "collection_name", which is what is dereferenced later (e.g., in
    # parsers.VariantParsing
    mongo_db_and_collection_names_dict = {'db_name': 'VariantDatabase',
                'project_name': 'collect'}

    Project = annotation_project.AnnotationProject(IN_PATH,
                                                   OUT_PATH,
                                                   ANNOVAR_PATH,
                                                   mongo_db_and_collection_names_dict,
                                                   design_file=design_file)