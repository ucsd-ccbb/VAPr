import unittest
import io

import VAPr.annovar_output_parsing as ns_test


# Most test info taken from test_files/test_out_csv_path/real_files/X7/X7.raw.11_annotated.hg19_multianno.txt'
class TestAnnovarTxtParser(unittest.TestCase):
    def test__normalize_header(self):
        input_raw_headers = ["Chr", "Start", "End", "Ref", "Alt", "Func.knownGene", "Gene.knownGene",
                             "GeneDetail.knownGene", "ExonicFunc.knownGene", "AAChange.knownGene", "tfbsConsSites",
                             "cytoBand", "targetScanS", "genomicSuperDups", "esp6500siv2_all", "1000g2015aug_all",
                             "PopFreqMax", "1000G_ALL", "1000G_AFR", "1000G_AMR", "1000G_EAS", "1000G_EUR", "1000G_SAS",
                             "ExAC_ALL", "ExAC_AFR", "ExAC_AMR", "ExAC_EAS", "ExAC_FIN", "ExAC_NFE", "ExAC_OTH",
                             "ExAC_SAS", "ESP6500siv2_ALL", "ESP6500siv2_AA", "ESP6500siv2_EA", "CG46", "CLINSIG",
                             "CLNDBN", "CLNACC", "CLNDSDB", "CLNDSDBID", "cosmic70", "nci60", "Otherinfo"]
        expected_output = ["chr", "start", "end", "ref", "alt", "func_knowngene", "gene_knowngene",
                           "genedetail_knowngene", "exonicfunc_knowngene", "aachange_knowngene", "tfbsconssites",
                           "cytoband", "targetscans", "genomicsuperdups", "esp6500siv2_all", "1000g2015aug_all",
                           "popfreqmax", "1000g_all", "1000g_afr", "1000g_amr", "1000g_eas", "1000g_eur", "1000g_sas",
                           "exac_all", "exac_afr", "exac_amr", "exac_eas", "exac_fin", "exac_nfe", "exac_oth",
                           "exac_sas", "esp6500siv2_all", "esp6500siv2_aa", "esp6500siv2_ea", "cg46", "clinsig",
                           "clndbn", "clnacc", "clndsdb", "clndsdbid", "cosmic70", "nci60", "otherinfo"]
        real_output = ns_test.AnnovarTxtParser._normalize_header(input_raw_headers)
        self.assertListEqual(expected_output, real_output)

    # region _rewrite_value_if_special_header tests
    def test__rewrite_value_if_special_header_unchanged(self):
        input_val = expected_output = "moon"
        real_output = ns_test.AnnovarTxtParser._rewrite_value_if_special_header("blue", input_val)
        self.assertEqual(expected_output, real_output)

    def test__rewrite_value_if_special_header_chrom_unchanged(self):
        input_val = expected_output = "chr1"
        real_output = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(ns_test.AnnovarTxtParser.CHR_HEADER,
                                                                                input_val)
        self.assertEqual(expected_output, real_output)

    def test__rewrite_value_if_special_header_chrom_cleaned(self):
        real_output = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.CHR_HEADER, ns_test.AnnovarTxtParser.RAW_CHR_MT_VAL)
        self.assertEqual(ns_test.AnnovarTxtParser.STANDARDIZED_CHR_MT_VAL, real_output)

    def test__rewrite_value_if_special_header_to_float(self):
        input_val = "12"
        expected_output = 12.0
        real_output1 = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.THOU_G_2015_ALL_HEADER, input_val)
        self.assertEqual(expected_output, real_output1)

        real_output2 = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.ESP6500_ALL_HEADER, input_val)
        self.assertEqual(expected_output, real_output2)

        real_output3 = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.NCI60_HEADER, input_val)
        self.assertEqual(expected_output, real_output3)

    def test__rewrite_value_if_special_header_to_int(self):
        input_val = "12"
        expected_output = 12
        real_output1 = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.START_HEADER, input_val)
        self.assertEqual(expected_output, real_output1)

        real_output2 = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.END_HEADER, input_val)
        self.assertEqual(expected_output, real_output2)

    def test__rewrite_value_if_special_header_to_dict_w_score(self):
        input_val = "Score=890;Name=V$FOXJ2_01"
        expected_output = {"Score": 890.0,
                           "Name": "V$FOXJ2_01"}
        real_output1 = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.GENOMIC_SUPERDUPS_HEADER, input_val)
        self.assertDictEqual(expected_output, real_output1)

        real_output2 = ns_test.AnnovarTxtParser._rewrite_value_if_special_header(
            ns_test.AnnovarTxtParser.TFBS_CONS_SITES_HEADER, input_val)
        self.assertDictEqual(expected_output, real_output2)

    # endregion

    def test__parse_to_dict_with_score_key(self):
        input_str = "Score=890;Name=V$FOXJ2_01"  # taken from a genomicsuperdups value for a real variant
        expected_output = {"Score": 890.0,
                           "Name": "V$FOXJ2_01"}
        real_output = ns_test.AnnovarTxtParser._parse_to_dict_with_score_key(input_str)
        self.assertDictEqual(expected_output, real_output)

    def test_read_chunk_of_annotations_to_dicts_list(self):
        # NB: first 4 variants are ignored because test starts on chunk index = 1 (i.e., the second chunk) and the
        # chunk size is 4.  Likewise, since the test only processes a single 4-variant chunk, the last 2 variants
        # are ignored.  This is as expected and demonstrates that the chunking function is working correctly.
        annovar_txt_str = """Chr	Start	End	Ref	Alt	Func.knownGene	Gene.knownGene	GeneDetail.knownGene	ExonicFunc.knownGene	AAChange.knownGene	tfbsConsSites	cytoBand	targetScanS	genomicSuperDups	esp6500siv2_all	1000g2015aug_all	PopFreqMax	1000G_ALL	1000G_AFR	1000G_AMR	1000G_EAS	1000G_EUR	1000G_SAS	ExAC_ALL	ExAC_AFR	ExAC_AMR	ExAC_EAS	ExAC_FIN	ExAC_NFE	ExAC_OTH	ExAC_SAS	ESP6500siv2_ALL	ESP6500siv2_AA	ESP6500siv2_EA	CG46	CLINSIG	CLNDBN	CLNACC	CLNDSDB	CLNDSDBID	cosmic70	nci60	Otherinfo
chrM	516	517	CA	-	upstream	DQ582201,JB137816	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	43.70	2	chrM	515	.	GCA	G	43.70	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=21.85;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:80,6,0	./.:0,0
chrM	1890	1890	G	A	ncRNA_exonic	TVAS5	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	56.74	2	chrM	1890	.	G	A	56.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=43.64;MQ0=0;QD=28.37;SOR=2.303	GT:AD:DP:GQ:PL	1/1:0,2:2:6:84,6,0	./.:0,0
chrM	6262	6262	G	A	ncRNA_exonic	BC018860	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	185.90	5	chrM	6262	.	G	A	185.90	.	AC=2;AF=1.00;AN=2;DP=5;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=29.97;SOR=1.981	GT:AD:DP:GQ:PL	1/1:0,5:5:15:214,15,0	./.:0,0
chrM	8698	8698	G	A	ncRNA_exonic	OK/SW-cl.16	.	.	.	Score=908;Name=V$FOXD3_01	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	107.28	3	chrM	8698	.	G	A	107.28	.	AC=2;AF=1.00;AN=2;DP=4;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=49.57;MQ0=0;QD=26.82;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,3:3:9:135,9,0	./.:0,0
chrM	146	146	T	C	upstream;downstream	JB137816;DQ582201	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	62.74	2	chrM	146	rs370482130	T	C	62.74	.	AC=2;AF=1.00;AN=2;DB;DP=3;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=20.91;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:90,6,0	./.:0,0
chrM	150	150	T	C	upstream;downstream	JB137816;DQ582201	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	62.74	2	chrM	150	.	T	C	62.74	.	AC=2;AF=1.00;AN=2;DP=3;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=20.91;SOR=0.693	GT:AD:DP:GQ:PL	./.:0,0	1/1:0,2:2:6:90,6,0
chr1	195	195	C	T	upstream;downstream	JB137816;DQ582201	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	61.74	2	chrM	195	.	C	T	61.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=30.87;SOR=0.693	GT:AD:DP:GQ:PL	1/1:0,2:2:6:89,6,0	1/1:0,3:3:9:135,9,0
chrM	410	410	A	T	upstream	DQ582201,JB137816	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	62.74	2	chrM	410	.	A	T	62.74	.	AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=31.37;SOR=0.693	GT:AD:DP:GQ:PL	./.:0,0	./.:0,0
chr2	626	626	G	A	ncRNA_exonic	BC018860	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	185.90	5	chrM	6262	.	G	A	185.90	.	AC=2;AF=1.00;AN=2;DP=5;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=60.00;MQ0=0;QD=29.97;SOR=1.981	GT:AD:DP:GQ:PL	1/1:0,5:5:15:214,15,0	./.:0,0
chr4	698	698	G	A	ncRNA_exonic	OK/SW-cl.16	.	.	.	Score=908;Name=V$FOXD3_01	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	.	1	107.28	3	chrM	8698	.	G	A	107.28	.	AC=2;AF=1.00;AN=2;DP=4;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=49.57;MQ0=0;QD=26.82;SOR=1.179	GT:AD:DP:GQ:PL	1/1:0,3:3:9:135,9,0	./.:0,0
"""

        expected_hgvs_list = ["chrMT:g.146T>C", "chrMT:g.150T>C", "chr1:g.195C>T", "chrMT:g.410A>T"]
        # The first variant exists only in sample 1, the second only in sample 2, the third in both, and the fourth
        # in neither.
        expected_dicts_list = [
            {'chr': 'chrMT', 'start': 146, 'end': 146, 'ref': 'T', 'alt': 'C', 'func_knowngene': 'upstream;downstream',
             'gene_knowngene': 'JB137816;DQ582201', 'hgvs_id': 'chrMT:g.146T>C', 'samples': [
                {'sample_id': 'test_sample1', 'genotype': '1/1', 'genotype_subclass_by_class': {'homozygous': 'alt'},
                 'filter_passing_reads_count': 2, 'genotype_likelihoods': [90.0, 6.0, 0.0], 'AD': [0, 2]}]},
            {'chr': 'chrMT', 'start': 150, 'end': 150, 'ref': 'T', 'alt': 'C', 'func_knowngene': 'upstream;downstream',
             'gene_knowngene': 'JB137816;DQ582201', 'hgvs_id': 'chrMT:g.150T>C', 'samples': [
                {'sample_id': 'test_sample2', 'genotype': '1/1', 'genotype_subclass_by_class': {'homozygous': 'alt'},
                 'filter_passing_reads_count': 2, 'genotype_likelihoods': [90.0, 6.0, 0.0], 'AD': [0, 2]}]},
            {'chr': 'chr1', 'start': 195, 'end': 195, 'ref': 'C', 'alt': 'T', 'func_knowngene': 'upstream;downstream',
             'gene_knowngene': 'JB137816;DQ582201', 'hgvs_id': 'chr1:g.195C>T', 'samples': [
                {'sample_id': 'test_sample1', 'genotype': '1/1', 'genotype_subclass_by_class': {'homozygous': 'alt'},
                 'filter_passing_reads_count': 2, 'genotype_likelihoods': [89.0, 6.0, 0.0], 'AD': [0, 2]},
                {'sample_id': 'test_sample2', 'genotype': '1/1', 'genotype_subclass_by_class': {'homozygous': 'alt'},
                 'filter_passing_reads_count': 3, 'genotype_likelihoods': [135.0, 9.0, 0.0], 'AD': [0, 3]}]},
            {'chr': 'chrMT', 'start': 410, 'end': 410, 'ref': 'A', 'alt': 'T', 'func_knowngene': 'upstream',
             'gene_knowngene': 'DQ582201,JB137816', 'hgvs_id': 'chrMT:g.410A>T', 'samples': []}]

        input_txt_stream = io.StringIO(annovar_txt_str)
        real_hgvs_list, real_dict_list = ns_test.AnnovarTxtParser.read_chunk_of_annotations_to_dicts_list(
            input_txt_stream, ['test_sample1', 'test_sample2'], 1, 4)
        self.assertListEqual(expected_hgvs_list, real_hgvs_list)
        self.assertListEqual(expected_dicts_list, real_dict_list)

    def test__parse_single_variant_record(self):
        input_headers_list = ['chr', 'start', 'end', 'ref', 'alt', 'func_knowngene', 'gene_knowngene',
                              'genedetail_knowngene', 'exonicfunc_knowngene', 'aachange_knowngene', 'tfbsconssites',
                              'cytoband', 'targetscans', 'genomicsuperdups', 'esp6500siv2_all', '1000g2015aug_all',
                              'popfreqmax', '1000g_all', '1000g_afr', '1000g_amr', '1000g_eas', '1000g_eur',
                              '1000g_sas', 'exac_all', 'exac_afr', 'exac_amr', 'exac_eas', 'exac_fin', 'exac_nfe',
                              'exac_oth', 'exac_sas', 'esp6500siv2_all', 'esp6500siv2_aa', 'esp6500siv2_ea', 'cg46',
                              'clinsig', 'clndbn', 'clnacc', 'clndsdb', 'clndsdbid', 'cosmic70', 'nci60', 'otherinfo']
        input_fields_list = ['chrM', '1890', '1890', 'G', 'A', 'ncRNA_exonic', 'TVAS5', '.', '.', '.', '.', '.', '.',
                             '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.',
                             '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '.', '1', '56.74', '2', 'chrM', '1890',
                             '.', 'G', 'A', '56.74', '.',
                             'AC=2;AF=1.00;AN=2;DP=2;FS=0.000;MLEAC=2;MLEAF=1.00;MQ=43.64;MQ0=0;QD=28.37;SOR=2.303',
                             'GT:AD:DP:GQ:PL', '1/1:0,2:2:6:84,6,0', './.:0,0', '0/1:0,40:40:99:1494,119,0']
        input_sample_names_list = ["test_sample1", "test_sample2", "test_sample3"]
        expected_hgvs_id = "chrMT:g.1890G>A"
        expected_dict = {'chr': 'chrMT',
                         'start': 1890,
                         'end': 1890,
                         'ref': 'G',
                         'alt': 'A',
                         'func_knowngene': 'ncRNA_exonic',
                         'gene_knowngene': 'TVAS5',
                         'hgvs_id': 'chrMT:g.1890G>A',
                         'samples': [{'sample_id': 'test_sample1',
                                      'genotype': '1/1',
                                      'genotype_subclass_by_class': {'homozygous': 'alt'},
                                      'filter_passing_reads_count': 2,
                                      'genotype_likelihoods': [84.0, 6.0, 0.0],
                                      'AD': [0, 2]
                                      },
                                     {'sample_id': 'test_sample3',
                                      'genotype': '0/1',
                                      'genotype_subclass_by_class': {'heterozygous': 'reference'},
                                      'filter_passing_reads_count': 40,
                                      'genotype_likelihoods': [1494.0, 119.0, 0.0],
                                      'AD': [0, 40]
                                      }]
                         }

        real_hgvs_id, real_dict = ns_test.AnnovarTxtParser._parse_single_variant_record(
            input_headers_list, input_fields_list, input_sample_names_list)
        self.assertEqual(expected_hgvs_id, real_hgvs_id)
        self.assertDictEqual(expected_dict, real_dict)


class TestAnnovarAnnotatedVariant(unittest.TestCase):
    # region _list_has_valid_content tests
    def test__list_has_valid_content_true(self):
        self.assertTrue(ns_test.AnnovarAnnotatedVariant._list_has_valid_content(["a", 1, None]))

    def test__list_has_valid_content_false_all_none(self):
        self.assertFalse(ns_test.AnnovarAnnotatedVariant._list_has_valid_content([None, None]))

    def test__list_has_valid_content_false_empty(self):
        self.assertFalse(ns_test.AnnovarAnnotatedVariant._list_has_valid_content([]))

    # endregion

    def test_make_per_variant_annotation_dict(self):
        input_annovar_fields = {'chr': 'chrMT',
                                'start': 1890,
                                'end': 1890,
                                'ref': 'G',
                                'alt': 'A',
                                'func_knowngene': 'ncRNA_exonic',
                                'gene_knowngene': 'TVAS5'}
        input_hgvs_id = "chrMT:g.1890G>A"
        input_format_string = "GT:AD:DP:GQ:PL"
        input_fields_strings_by_sample_name = {"test_sample1": "1/1:0,2:2:6:84,6,0",
                                               "test_sample2": "./.:0,0",
                                               "test_sample3": "0/1:0,40:40:99:1494,119,0"}

        expected_output = {'chr': 'chrMT',
                           'start': 1890,
                           'end': 1890,
                           'ref': 'G',
                           'alt': 'A',
                           'func_knowngene': 'ncRNA_exonic',
                           'gene_knowngene': 'TVAS5',
                           'hgvs_id': 'chrMT:g.1890G>A',
                           'samples': [{'sample_id': 'test_sample1',
                                        'genotype': '1/1',
                                        'genotype_subclass_by_class': {'homozygous': 'alt'},
                                        'filter_passing_reads_count': 2,
                                        'genotype_likelihoods': [84.0, 6.0, 0.0],
                                        'AD': [0, 2]
                                        },
                                       {'sample_id': 'test_sample3',
                                        'genotype': '0/1',
                                        'genotype_subclass_by_class': {'heterozygous': 'reference'},
                                        'filter_passing_reads_count': 40,
                                        'genotype_likelihoods': [1494.0, 119.0, 0.0],
                                        'AD': [0, 40]
                                        }]
                           }

        real_output = ns_test.AnnovarAnnotatedVariant.make_per_variant_annotation_dict(
            input_annovar_fields, input_hgvs_id, input_format_string, input_fields_strings_by_sample_name)
        self.assertDictEqual(expected_output, real_output)

    # region _make_per_sample_annotation_dict tests
    def test__make_per_sample_annotation_dict_none(self):
        """Ensure that if input genotype fields string doesn't contain any real info, None is returned."""

        real_output = ns_test.AnnovarAnnotatedVariant._make_per_sample_annotation_dict("test_sample1", "GT:AD:DP:GQ:PL",
                                                                                       "./.:.:.:.:.")
        self.assertIsNone(real_output)

    def test__make_per_sample_annotation_dict_maximum_content(self):
        """Ensure that if info exists in the format string for all relevant fields, it is put into the dictionary."""

        expected_output = {ns_test.AnnovarAnnotatedVariant.SAMPLE_ID_KEY: "test_sample1",
                           ns_test.AnnovarAnnotatedVariant.GENOTYPE_KEY: "1/1",
                           ns_test.AnnovarAnnotatedVariant.GENOTYPE_SUBCLASS_BY_CLASS_KEY: {'homozygous': 'alt'},
                           ns_test.AnnovarAnnotatedVariant.FILTER_PASSING_READS_COUNT_KEY: 2,
                           ns_test.AnnovarAnnotatedVariant.GENOTYPE_LIKELIHOODS_KEY: [89.0, 6.0, 0.0],
                           ns_test.AnnovarAnnotatedVariant.ALLELE_DEPTH_KEY: [0, 2]}
        real_output = ns_test.AnnovarAnnotatedVariant._make_per_sample_annotation_dict("test_sample1", "GT:AD:DP:GQ:PL",
                                                                                       "1/1:0,2:2:6:89,6,0")
        self.assertDictEqual(expected_output, real_output)

    def test__make_per_sample_annotation_dict_minimum_content(self):
        """Ensure that if info isn't in the format string for relevant fields, a minimal dictionary is made."""

        expected_output = {ns_test.AnnovarAnnotatedVariant.SAMPLE_ID_KEY: "test_sample1",
                           ns_test.AnnovarAnnotatedVariant.GENOTYPE_KEY: None}
        real_output = ns_test.AnnovarAnnotatedVariant._make_per_sample_annotation_dict("test_sample1", "GQ", "6")
        self.assertDictEqual(expected_output, real_output)

        # endregion

# class TestCytoBand(unittest.TestCase):
#     def test_fill(self):
#         cytoBand = ['1p36.33', '16p11.1', '16q11.2', '16q21', 'Xp22.32', 'Xp22.2', 'Xp22.11',
#                     'Xq12', 'Xq13.1']
#         cytoBand_dict = [
#                                 {'Band': 'p', 'Name': '1p36.33', 'Chromosome': '1', 'Region': 36, 'Sub_Band': 33},
#                                 {'Band': 'p', 'Name': '16p11.1', 'Chromosome': '16', 'Region': 11, 'Sub_Band': 1},
#                                 {'Band': 'q', 'Name': '16q11.2', 'Chromosome': '16', 'Region': 11, 'Sub_Band': 2},
#                                 {'Band': 'q', 'Name': '16q21', 'Chromosome': '16', 'Region': 21},
#                                 {'Band': 'p', 'Name': 'Xp22.32', 'Chromosome': 'X', 'Region': 22, 'Sub_Band': 32},
#                                 {'Band': 'p', 'Name': 'Xp22.2', 'Chromosome': 'X', 'Region': 22, 'Sub_Band': 2},
#                                 {'Band': 'p', 'Name': 'Xp22.11', 'Chromosome': 'X', 'Region': 22, 'Sub_Band': 11},
#                                 {'Band': 'q', 'Name': 'Xq12', 'Chromosome': 'X', 'Region': 12},
#                                 {'Band': 'q', 'Name': 'Xq13.1', 'Chromosome': 'X', 'Region': 13, 'Sub_Band': 1}
#                              ]
#
#         cytoband_list = []
#         for i in cytoBand:
#             cyto = ns_test.CytoBand(i)
#             cytoband_list.append(cyto.processed)
#
#         for idx,_ in enumerate(cytoBand_dict):
#             print(cytoband_list[idx])
#             self.assertEqual(cytoband_list[idx].keys(), cytoBand_dict[idx].keys())
