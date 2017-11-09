import sys
import os
import shlex
import subprocess
import logging
import pandas
from VAPr.vcf_mappings_maker import SingleVcfFileMappingMaker

__author__ = 'Adam Mark<a1mark@ucsd.edu>'

logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class MergeVcfs:
    def __init__(self, input_dir, output_dir, analysis_name, design_file, vcf_file_extension):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.vcf_name = analysis_name
        self.output_vcf_path = os.path.join(self.output_dir, self.vcf_name + ".vcf")
        self.vcf_file_extension = vcf_file_extension
        self.raw_vcf_path_list = self.get_vcf_file_paths_list(input_dir, design_file)

    def get_vcf_file_paths_list(self, input_dir, design_file):
        if design_file is not None:
            design_df = pandas.read_csv(design_file)
            vcf_file_paths_list = design_df[0].tolist()
        else:
            vcf_file_paths_list = self._get_vcf_file_paths_list_in_directory(input_dir)

        return vcf_file_paths_list

    def _get_vcf_file_paths_list_in_directory(self, base_dir):
        vcf_file_paths_list = []
        walker = os.walk(base_dir)
        for folder, _, files in walker:
            for curr_file in files:
                if curr_file.endswith(self.vcf_file_extension):
                    full_path_single_file = os.path.join(os.path.abspath(folder), curr_file)
                    vcf_file_paths_list.append(full_path_single_file)

        return vcf_file_paths_list

    def merge_vcfs(self):
        """Merge vcf files into single multisample vcf, bgzip and index merged vcf file."""
        try:
            os.mkdir(self.output_dir)
        except OSError:
            logging.info('Output directory %s for analysis already exists; using existing directory' %
                         self.output_dir)

        if len(self.raw_vcf_path_list) > 1:
            bgzipped_vcf_path_list = set([self.bgzip_index_vcf(vcf) for vcf in self.raw_vcf_path_list])
            single_vcf_path = self.execute_merge(bgzipped_vcf_path_list, self.output_vcf_path)
        else:
            single_vcf_path = self.raw_vcf_path_list[0]

        return SingleVcfFileMappingMaker(single_input_file_path=single_vcf_path,
                                          input_dir=self.input_dir,
                                          out_dir=self.output_dir,
                                          sample_id='infer',
                                          sample_id_type='files',
                                          extra_data=None).vcf_mapping_dict

    def execute_merge(self, vcf_list, out_vcf):
        merge_cmd_string = self._build_merge_vcf_command_str(vcf_list)
        merge_args = shlex.split(merge_cmd_string)
        with open(out_vcf, 'w') as outfile:
            p=subprocess.Popen(merge_args, stdout=outfile, stderr=subprocess.PIPE)
            p.communicate()
        return out_vcf

    def bgzip_index_vcf(self, vcf_path):
        """bgzip and index each vcf so it can be verged with bcftools."""

        if vcf_path.endswith(".vcf.gz"):
            return vcf_path
        else:
            vcf_gz = vcf_path + ".gz"
            bgzip_cmd_string = self._build_bgzip_vcf_command_str(vcf_path)
            bgzip_args = shlex.split(bgzip_cmd_string)
            with open(vcf_gz, "w") as outfile:
                p=subprocess.Popen(bgzip_args, stdout=outfile, stderr=subprocess.PIPE)
                p.communicate()

            index_cmd_string = self._build_index_vcf_command_str(vcf_gz)
            index_args = shlex.split(index_cmd_string)
            subprocess.call(index_args)
        return vcf_gz

    def _build_merge_vcf_command_str(self, raw_vcf_path_list):
        """Generate command string to merge vcf files into single multisample vcf."""

        command = " ".join([
            'bcftools merge', " ".join(raw_vcf_path_list)
            ])
        return command

    def _build_bgzip_vcf_command_str(self, vcf_path):
        """Generate command string to bgzip vcf file."""

        command = " ".join([
            'bgzip -c', vcf_path
            ])
        return command


    def _build_index_vcf_command_str(self, bgzipped_vcf):
        """Generate command string to index vcf file."""

        command = " ".join([
            'tabix -p vcf', bgzipped_vcf
            ])
        return command
