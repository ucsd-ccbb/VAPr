import os
import shlex
import subprocess
import logging
import pandas
from vcf import Reader
from shutil import copyfile

__author__ = 'Adam Mark<a1mark@ucsd.edu>'


def merge_vcfs(input_dir, output_dir, design_file, project_name, vcf_file_extension):
    """Merge vcf files into single multisample vcf, bgzip and index merged vcf file."""

    output_vcf_path = os.path.join(output_dir, project_name + ".vcf")
    raw_vcf_path_list = _get_vcf_file_paths_list(input_dir, design_file, vcf_file_extension)

    try:
        os.mkdir(output_dir)
    except OSError:
        logging.info('Output directory %s for analysis already exists; using existing directory' %
                     output_dir)

    if len(raw_vcf_path_list) > 1:
        bgzipped_vcf_path_list = set([_bgzip_index_vcf(vcf) for vcf in raw_vcf_path_list])
        single_vcf_path = _execute_merge(bgzipped_vcf_path_list, output_vcf_path)
    else:
        single_vcf_path = os.path.join(output_dir, project_name + vcf_file_extension)
        copyfile(raw_vcf_path_list[0], single_vcf_path)

    reader = Reader(open(single_vcf_path, 'r'))
    annovar_output_basename = os.path.splitext(os.path.basename(single_vcf_path))[0] + '_annotated'

    return single_vcf_path, annovar_output_basename, reader.samples


def _get_vcf_file_paths_list(input_dir, design_file, vcf_file_extension):
    if design_file is not None:
        design_df = pandas.read_csv(design_file)
        vcf_file_paths_list = design_df['Sample_Names'].tolist()
    else:
        vcf_file_paths_list = _get_vcf_file_paths_list_in_directory(input_dir, vcf_file_extension)

    return vcf_file_paths_list


def _get_vcf_file_paths_list_in_directory(base_dir, vcf_file_extension):
    vcf_file_paths_list = []
    walker = os.walk(base_dir)
    for folder, _, files in walker:
        for curr_file in files:
            if curr_file.endswith(vcf_file_extension):
                full_path_single_file = os.path.join(os.path.abspath(folder), curr_file)
                vcf_file_paths_list.append(full_path_single_file)

    return vcf_file_paths_list


def _execute_merge(vcf_list, out_vcf):
    merge_cmd_string = _build_merge_vcf_command_str(vcf_list)
    merge_args = shlex.split(merge_cmd_string)
    with open(out_vcf, 'w') as outfile:
        p = subprocess.Popen(merge_args, stdout=outfile, stderr=subprocess.PIPE)
        p.communicate()
    return out_vcf


def _bgzip_index_vcf(vcf_path):
    """bgzip and index each vcf so it can be verged with bcftools."""

    if vcf_path.endswith(".vcf.gz"):
        return vcf_path
    else:
        vcf_gz = vcf_path + ".gz"
        bgzip_cmd_string = _build_bgzip_vcf_command_str(vcf_path)
        bgzip_args = shlex.split(bgzip_cmd_string)
        with open(vcf_gz, "w") as outfile:
            p = subprocess.Popen(bgzip_args, stdout=outfile, stderr=subprocess.PIPE)
            p.communicate()

        index_cmd_string = _build_index_vcf_command_str(vcf_gz)
        index_args = shlex.split(index_cmd_string)
        subprocess.call(index_args)
    return vcf_gz


def _build_merge_vcf_command_str(raw_vcf_path_list):
    """Generate command string to merge vcf files into single multisample vcf."""

    command = " ".join([
        'bcftools merge', " ".join(raw_vcf_path_list)
        ])
    return command


def _build_bgzip_vcf_command_str(vcf_path):
    """Generate command string to bgzip vcf file."""

    command = " ".join([
        'bgzip -c', vcf_path
        ])
    return command


def _build_index_vcf_command_str(bgzipped_vcf):
    """Generate command string to index vcf file."""

    command = " ".join([
        'tabix -p vcf', bgzipped_vcf
        ])
    return command
