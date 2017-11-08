import sys
import os
import shlex
import subprocess
import logging
from VAPr.vcf_mappings_maker import VcfMappingsMaker, SingleVcfFileMappingMaker

__author__ = 'Adam Mark<a1mark@ucsd.edu>'

logger = logging.getLogger()
logger.setLevel(logging.INFO)
try:
    logger.handlers[0].stream = sys.stdout
except:
    pass


class MergeVcfs:
    def __init__(self, input_dir, output_dir, list_of_vcf_mapping_dicts, analysis_name):
        self.input_dir = input_dir
        self.output_dir = output_dir
        self.vcf_name = analysis_name
        # User must provide analysis_name without .vcf extension: if one single vcf file
        # If multiple vcf files, analysis name will be name of merged vcf file
        self.output_vcf_path = os.path.join(self.output_dir, self.vcf_name + '.vcf')
        self.list_of_vcf_mapping_dicts = list_of_vcf_mapping_dicts
        self.raw_vcf_path_list =  [vcf['raw_vcf_file_full_path'] for vcf in list_of_vcf_mapping_dicts]

    def merge_vcfs(self):
        """Merge vcf files into single multisample vcf, bgzip and index merged vcf file."""

        try:
            os.mkdir(self.output_dir)
        except OSError:
            logging.info('Output directory %s for analysis already exists; using existing directory' %
                         self.output_dir)

        if len(self.raw_vcf_path_list) > 1:
            bgzipped_vcf_path_list = set([self.bgzip_index_vcf(vcf) for vcf in self.raw_vcf_path_list])
            self.execute_merge(bgzipped_vcf_path_list, self.output_vcf_path)
            return [SingleVcfFileMappingMaker(single_input_file_path=self.output_vcf_path + ".gz",
                                         input_dir=self.input_dir,
                                         out_dir=self.output_dir,
                                         sample_id='infer',
                                         sample_id_type='files',
                                         extra_data=None).vcf_mapping_dict]
        else:
            return self.list_of_vcf_mapping_dicts


    def execute_merge(self, vcf_list, out_vcf):
        merge_cmd_string = self._build_merge_vcf_command_str(vcf_list)
        merge_args = shlex.split(merge_cmd_string)
        with open(out_vcf, 'w') as outfile:
            p=subprocess.Popen(merge_args, stdout=outfile, stderr=subprocess.PIPE)
            p.communicate()
        bgzip_index_vcf(out_vcf)

    def bgzip_index_vcf(self, vcf_path):
        """bgzip and index each vcf so it can be verged with bcftools."""
        if os.path.isfile(vcf_path + ".gz"):
            return vcf_path + ".gz"
        else:
            bgzip_cmd_string = self._build_bgzip_vcf_command_str(vcf_path)
            bgzip_args = shlex.split(bgzip_cmd_string)
            with open(vcf_path + ".gz", "w") as outfile:
                p=subprocess.Popen(bgzip_args, stdout=outfile, stderr=subprocess.PIPE)
                p.communicate()

            index_cmd_string = self._build_index_vcf_command_str(vcf_path + '.gz')
            index_args = shlex.split(index_cmd_string)
            subprocess.call(index_args)

            return vcf_path + ".gz"

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
