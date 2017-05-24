from __future__ import division, print_function
import os
import sys
import pandas
import logging
logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class ProjectData(object):

    supported_build_vers = ['hg19', 'hg18', 'hg38']

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None):

        """
        Base class for the project instantiation classes that will run the necessary functions. Mostly a container.
        """

        self.input_dir = input_dir
        self.output_csv_path = output_dir
        self.vcf_files = self.find_vcfs()
        self.design_file = design_file
        self.annovar = annovar_path
        self.project_data = project_data
        self.buildver = build_ver
        self.mapping = self.get_mapping() #[(k, v['vcf']) for k, v in self.get_mapping().items()]

    def get_mapping(self):

        if self.design_file:
            design_df = pandas.read_csv(self.design_file)
            sample_to_vcf_map = self.check_file_existance_return_mapping(design_df)

        else:
            sample_to_vcf_map = {os.path.splitext(os.path.basename(vcf))[0]: vcf for vcf in self.vcf_files}

        return sample_to_vcf_map
        # return {os.path.join(self.input_dir, vcf): self.output_csv_path + os.path.splitext(os.path.basename(vcf))[0] +
        #        '_annotated' for vcf in sample_to_vcf_map.values()}

    def find_vcfs(self):
        vcf_files = os.listdir(self.input_dir)
        return [vcf for vcf in vcf_files if vcf.endswith('vcf')]

    def check_file_existance_return_mapping(self, design_df):

        design_file_mapping = design_df.set_index('Sample_Names').T.to_dict()
        for sample in design_file_mapping.keys():
            if not os.path.exists(os.path.join(self.input_dir, sample)):
                raise NameError('Could not find directory named %s as provided in design file' % sample)

            vcf_files = [i for i in os.listdir(os.path.join(self.input_dir, sample)) if i.endswith('.vcf')]
            logging.info('Found %i unique vcf files for sample %s' % (len(set(vcf_files)), sample))

            design_file_mapping[sample]['vcf_csv'] = [(vcf_file, os.path.splitext(os.path.basename(vcf_file))[0] +
                                                      '_annotated') for vcf_file in vcf_files]  # God bless listcomps

        return design_file_mapping

    def check_ver(self, build_ver):

        if not build_ver:
            self.buildver = 'hg19'  # Default genome build version

        if build_ver not in self.supported_build_vers:
            raise ValueError('Build version must not recognized. Supported builds are'
                             ' %s, %s, %s' % (self.supported_build_vers[0],
                                              self.supported_build_vers[1],
                                              self.supported_build_vers[2]))
        else:
            self.buildver = build_ver

        return self.buildver


class AnnotationProject(ProjectData):

    def __init__(self,
                 input_dir,
                 output_dir,
                 annovar_path,
                 project_data,
                 design_file=None,
                 build_ver=None):

        """ Class that implements the API and the major annotation/saving methods  """

        super(AnnotationProject, self).__init__(input_dir,
                                                output_dir,
                                                annovar_path,
                                                project_data,
                                                design_file=design_file,
                                                build_ver=build_ver)

        from annovar import AnnovarWrapper
        self.annovar_wrapper = AnnovarWrapper(self.input_dir,
                                              self.output_csv_path,
                                              self.annovar,
                                              self.project_data,
                                              design_file=self.design_file,
                                              build_ver=self.buildver)

        from parsers import VariantParsing
        self.annotator_wrapper = VariantParsing(self.input_dir,
                                                self.output_csv_path,
                                                self.annovar,
                                                self.project_data,
                                                design_file=self.design_file,
                                                build_ver=self.buildver)

    def download_dbs(self, all_dbs=True, dbs=None):
        """ Wrapper around Annovar database downloading function """
        self.annovar_wrapper.download_dbs(all_dbs=all_dbs, dbs=dbs)

    def run_annovar(self):
        """ Wrapper around multiprocess Annovar annotation  """
        self.annovar_wrapper.run_annovar()

    def annotate_and_save(self):
        """ Wrapper around annotation runner  """
        self.annotator_wrapper.annotate_and_save(buffer=False)

    def parallel_annotation_and_saving(self, n_processes=4, verbose=1):
        """ Wrapper around parallel annotation multiprocess runner  """
        self.annotator_wrapper.parallel_annotation(n_processes=n_processes, verbose=verbose)
