from vapr import base


if __name__ == '__main__':

        # Directory of input files to be annotated
        IN_PATH = "/mnt/data1/carlom/samples/"

        # file directory
        OUT_PATH = "/mnt/data1/carlom/csv_results/"
        # Location of your annovar dowload. The folder should contain the following files/directories:
        ANNOVAR_PATH = '/mnt/data1/carlom/annovar/'

        # Design File (optional)
        design_file = '/mnt/data1/carlom/design_file_three_samples.csv'

        # Databse and Collection names (optional)
        proj_data = {'db_name': 'VariantDatabase',
                    'project_name': 'collect'}

        Project = base.AnnotationProject(IN_PATH,
                                         OUT_PATH,
                                         ANNOVAR_PATH,
                                         proj_data,
                                         design_file=design_file)