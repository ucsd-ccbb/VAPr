# standard libraries
import os
import unittest

# project-specific libraries
import VAPr.annovar_running as ns_test


class TestAnnovarWrapper(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls._base_dir = os.getcwd()
        cls._annovar_install_path = os.path.join(cls._base_dir, 'test_files/annovar_dir')
        cls._genome_build_version = "hg19"

    # region _get_annovar_dbs_to_use_for_build_version tests
    def test__get_annovar_dbs_to_use_for_build_version_hg19(self):
        dbs_ordered_dict = ns_test.AnnovarWrapper._get_annovar_dbs_to_use_for_build_version("hg19")
        self.assertDictEqual(ns_test.AnnovarWrapper.hg_19_databases, dbs_ordered_dict)

    def test__get_annovar_dbs_to_use_for_build_version_hg38(self):
        dbs_ordered_dict = ns_test.AnnovarWrapper._get_annovar_dbs_to_use_for_build_version("hg38")
        self.assertDictEqual(ns_test.AnnovarWrapper.hg_19_databases, dbs_ordered_dict)

    def test__get_annovar_dbs_to_use_for_build_version_error(self):
        with self.assertRaises(ValueError):
            ns_test.AnnovarWrapper._get_annovar_dbs_to_use_for_build_version("hg08")

    # endregion

    def test__build_table_annovar_command_str(self):
        expected_output = ('perl {0}/table_annovar.pl /my/vcf/dir/X45.raw.22.vcf {0}/humandb/ --buildver hg19 '
                           '-out /my/output/dir -remove -protocol knownGene,1000g2015aug_all -operation g,f '
                           '-nastring . -otherinfo -vcfinput'.format(self._annovar_install_path))
        wrapper = ns_test.AnnovarWrapper(self._annovar_install_path, self._genome_build_version)
        real_output = wrapper._build_table_annovar_command_str("/my/vcf/dir/X45.raw.22.vcf", "/my/output/dir")
        self.maxDiff = None
        self.assertEqual(expected_output, real_output)

    def test__build_annovar_database_download_command_str(self):
        expected_output = ["perl {0}/annotate_variation.pl -build hg19 -downdb -webfrom annovar knownGene "
                           "{0}/humandb/".format(self._annovar_install_path)]
        wrapper = ns_test.AnnovarWrapper(self._annovar_install_path, self._genome_build_version, ['knownGene'])
        real_output = wrapper._build_annovar_database_download_command_str()
        self.assertListEqual(expected_output, real_output)
