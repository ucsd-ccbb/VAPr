import os
import shlex
import subprocess
import shutil

__author__ = 'Adam Mark<a1mark@ucsd.edu>'

VCF_EXTENSION = ".vcf"
BGZIP_EXTENSION = ".gz"
BGZIPPED_VCF_EXTENSION = VCF_EXTENSION + BGZIP_EXTENSION

def bgzip_and_index_vcf(vcf_path):
    """bgzip and index each vcf so it can be merged with bcftools."""

    if vcf_path.endswith(BGZIPPED_VCF_EXTENSION):
        # TODO: someday: check that the input is *really* bgzipped, rather than just gzipped
        # TODO: someday: check to for an index file for this compressed file, rather than assuming is there
        # something like if not os.path.isfile(vcf_path + ".gz.tbi")): raise ValueError
        bgzipped_vcf_path = vcf_path
    else:
        bgzipped_vcf_path = vcf_path + BGZIP_EXTENSION
        bgzip_cmd_string = _build_bgzip_vcf_command_str(vcf_path)
        bgzip_args = shlex.split(bgzip_cmd_string)
        with open(bgzipped_vcf_path, "w") as outfile:
            p = subprocess.Popen(bgzip_args, stdout=outfile, stderr=subprocess.PIPE)
            p.communicate()

        index_cmd_string = _build_index_vcf_command_str(bgzipped_vcf_path)
        index_args = shlex.split(index_cmd_string)
        subprocess.call(index_args)

    return bgzipped_vcf_path


# TODO: someday: refactor since raw_vcf_path_list and vcfs_gzipped are really mutually exclusive
def merge_vcfs(input_dir, output_dir, project_name, raw_vcf_path_list=None, vcfs_gzipped=False):
    """Merge vcf files into single multisample vcf, bgzip and index merged vcf file."""

    if raw_vcf_path_list is None:
        vcf_file_extension = BGZIPPED_VCF_EXTENSION if vcfs_gzipped else VCF_EXTENSION
        raw_vcf_path_list = _get_vcf_file_paths_list_in_directory(input_dir, vcf_file_extension)
        if len(raw_vcf_path_list) == 0:
            raise ValueError("No VCFs found with extension '{0}'.".format(vcf_file_extension))
    elif len(raw_vcf_path_list) == 0:
            raise ValueError("Input list of VCF files is empty.")

    if len(raw_vcf_path_list) > 1:
        bgzipped_vcf_path_list = set([bgzip_and_index_vcf(vcf_fp) for vcf_fp in raw_vcf_path_list])
        single_vcf_path = os.path.join(output_dir, project_name + VCF_EXTENSION)
        _merge_bgzipped_indexed_vcfs(bgzipped_vcf_path_list, single_vcf_path)
    else:
        single_vcf_path = os.path.join(output_dir, project_name + VCF_EXTENSION)
        try:
            shutil.copyfile(raw_vcf_path_list[0], single_vcf_path)
        except shutil.SameFileError:
            # I ran into a case where there was a single input file, AND it was already named after the project, AND
            # the input and output dirs were the same so it was already where it needed to be.  In this case, an
            # error is thrown because you can't copy a file to itself, but that's cool, so just ignore it.
            pass

    return single_vcf_path


def _get_vcf_file_paths_list_in_directory(base_dir, vcf_file_extension):
    vcf_file_paths_list = []
    walker = os.walk(base_dir)
    for folder, _, files in walker:
        for curr_file in files:
            if curr_file.endswith(vcf_file_extension):
                full_path_single_file = os.path.join(os.path.abspath(folder), curr_file)
                vcf_file_paths_list.append(full_path_single_file)

    return sorted(vcf_file_paths_list)


def _merge_bgzipped_indexed_vcfs(bgzipped_vcf_path_list, output_vcf_fp):
    merge_cmd_string = _build_merge_vcf_command_str(bgzipped_vcf_path_list)
    merge_args = shlex.split(merge_cmd_string)
    with open(output_vcf_fp, 'w') as outfile:
        p = subprocess.Popen(merge_args, stdout=outfile, stderr=subprocess.PIPE)
        p.communicate()


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
