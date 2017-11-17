# TODO: someday: I think these should probably be refactored to be in an external file; yaml, maybe?

SAMPLE_ID_SELECTOR = 'samples.sample_id'


def get_sample_id_filter(sample_name):
    return {SAMPLE_ID_SELECTOR: sample_name}


def get_any_of_sample_ids_filter(sample_names_list):
    return {SAMPLE_ID_SELECTOR: {'$in': sample_names_list}}


def make_de_novo_variants_filter(proband, ancestor1, ancestor2):
    """
    Function for de novo variant analysis. Can be performed on multisample files or or on data coming
    from a collection of files. In the former case, every sample contains the same variants, although they have
    differences in their allele frequency and read values. A de novo variant is defined as a variant that
    occurs only in the specified sample (sample1) and not on the other two (sample2, sample3). Occurrence is
    defined as having allele frequencies greater than [0, 0] ([REF, ALT]).
    """

    return {
            "$and":
                    [
                        get_sample_id_filter(proband),
                        {
                            "$and":
                                [
                                    {SAMPLE_ID_SELECTOR: {"$ne": ancestor1}},
                                    {SAMPLE_ID_SELECTOR: {"$ne": ancestor2}}
                                ]
                        }
                    ]
            }


def make_deleterious_compound_heterozygous_variants_filter(sample_ids_list=None):
    and_list = [
                    {"genotype_subclass_by_class.heterozygous": "compound"},
                    {"cadd.phred": {"$gte": 10}}
               ]

    result = _append_sample_id_constraint_if_needed(and_list, sample_ids_list)
    return result


def make_known_disease_variants_filter(sample_ids_list=None):
    """ Function for retrieving known disease variants by presence in Clinvar and Cosmic."""

    result = {
            "$or":
                    [
                        {
                            "$and":
                                [
                                    {"clinvar.rcv.accession": {"$exists": True}},
                                    {"clinvar.rcv.clinical_significance": {"$nin": ["Benign", "Likely benign"]}}
                                ]
                        },
                        {"cosmic.cosmic_id": {"$exists": True}}
                    ]
        }

    if sample_ids_list is not None:
        result = _append_sample_id_constraint_if_needed([result], sample_ids_list)

    return result


def make_rare_deleterious_variants_filter(sample_ids_list=None):
    """ Function for retrieving rare, deleterious variants """

    and_list = [
                    {
                        "$or":
                            [
                                {"cadd.esp.af": {"$lt": 0.051}},
                                {"cadd.esp.af": {"$exists": False}}
                            ]
                    },
                    {
                        "$or":
                            [
                                {"func_knowngene": "exonic"},
                                {"func_knowngene": "splicing"}
                            ]
                    },
                    {"cadd.phred": {"$gte": 10}},
                    {"exonicfunc_knowngene": {"$ne": "synonymous SNV"}},
                    {"1000g2015aug_all": {"$lt": 0.051}}
                ]

    result = _append_sample_id_constraint_if_needed(and_list, sample_ids_list)
    return result


def _append_sample_id_constraint_if_needed(and_list, sample_ids_list):
    if sample_ids_list is not None:
        and_list.append(get_any_of_sample_ids_filter(sample_ids_list))
    return {"$and": and_list}
