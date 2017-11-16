# TODO: someday: I think these should probably be refactored to be in an external file; yaml, maybe?

SAMPLE_ID_SELECTOR = 'samples.sample_id'


def get_sample_id_filter(sample_name):
    return {SAMPLE_ID_SELECTOR: sample_name}


def get_any_of_sample_ids_filter(sample_names_list):
    return {SAMPLE_ID_SELECTOR: {'$in': sample_names_list}}


def make_rare_deleterious_variants_filter(sample_ids_list):
    """ Function for retrieving rare, deleterious variants """

    return {
            "$and":
                [
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
                    {"exonicfunc_knowngene": {"$ne": "synonymous SNV"}},
                    {"1000g2015aug_all": {"$lt": 0.051}},
                    get_any_of_sample_ids_filter(sample_ids_list)
                ]
        }


# TODO: the function this came from took sample ids, but didn't use them; should it?
def make_known_disease_variants_filter(sample_ids_list):
    """ Function for retrieving known disease variants by presence in Clinvar and Cosmic."""

    return {
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


def make_deleterious_compound_heterozygote_variants_filter(sample_ids_list):
    return {
            "$and":
                [
                    {"genotype_subclass_by_class.heterozygous": "compound"},
                    {"cadd.phred": {"$gte": 10}},
                    get_any_of_sample_ids_filter(sample_ids_list)
                ]
        }


# TODO: can we come up with some more informative names for these samples, in case users don't read docstring?
def make_de_novo_variants_filter(sample1, sample2, sample3):
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
                        get_sample_id_filter(sample1),
                        {
                            "$and":
                                [
                                    {SAMPLE_ID_SELECTOR: {"$ne": sample2}},
                                    {SAMPLE_ID_SELECTOR: {"$ne": sample3}}
                                ]
                        }
                    ]
            }


# def get_rare_deleterious_variants(collection, sample_names=None):
#     """ Function for retrieving rare, deleterious variants """
#
#     sample_ids_list = _construct_sample_ids_list(collection, sample_names)
#     filtered = collection.find(
#         {
#             "$and":
#                 [
#                     {
#                         "$or":
#                             [
#                                 {"cadd.esp.af": {"$lt": 0.051}},
#                                 {"cadd.esp.af": {"$exists": False}}
#                             ]
#                     },
#                     {
#                         "$or":
#                             [
#                                 {"func_knowngene": "exonic"},
#                                 {"func_knowngene": "splicing"}
#                             ]
#                     },
#                     {"exonicfunc_knowngene": {"$ne": "synonymous SNV"}},
#                     {"1000g2015aug_all": {"$lt": 0.051}},
#                     {_SAMPLE_ID_SELECTOR: {"$in": sample_ids_list}}
#
#                 ]
#         }
#     )
#
#     filtered = list(filtered)
#     logging.info('Variants found that match rarity criteria: {}'.format(len(filtered)))
#     return filtered
#
#
# def get_known_disease_variants(collection, samples=None):
#     """ Function for retrieving known disease variants by presence in Clinvar and Cosmic."""
#
#     if not samples:
#         samples = collection.distinct(_SAMPLE_ID_SELECTOR)
#     if not isinstance(samples, list):
#         samples = [samples]
#
#     filtered = collection.find(
#         {
#             "$or" :
#                     [
#                         {
#                             "$and":
#                                 [
#                                     {"clinvar.rcv.accession": {"$exists": True}},
#                                     {"clinvar.rcv.clinical_significance": {"$nin": ["Benign", "Likely benign"]}}
#                                 ]
#                         },
#                         {"cosmic.cosmic_id": {"$exists": True}}
#                     ]
#         }
#     )
#
#     filtered = list(filtered)
#     logging.info ('Variants found that match rarity criteria: {}'.format(len(filtered)))
#     return filtered
#
#
# def get_deleterious_compound_heterozygote_variants(collection, sample_names=None):
#     """ Function for retrieving deleterious compound heterozygote variants  """
#
#     sample_ids_list = _construct_sample_ids_list(collection, sample_names)
#     filtered = collection.find(
#         {
#             "$and":
#                 [
#                     {"genotype_subclass_by_class.heterozygous": "compound"},
#                     {"cadd.phred": {"$gte": 10}},
#                     {'samples.sample_id': {'$in': sample_ids_list}}
#                 ]
#         }
#     )
#
#     filtered = list(filtered)
#     logging.info('Variants found that match rarity criteria: {}'.format(len(filtered)))
#     return filtered
#
#
# def get_de_novo_variants(collection, sample1, sample2, sample3):
#     """
#     Function for de novo variant analysis. Can be performed on multisample files or or on data coming
#     from a collection of files. In the former case, every sample contains the same variants, although they have
#     differences in their allele frequency and read values. A de novo variant is defined as a variant that
#     occurs only in the specified sample (sample1) and not on the other two (sample2, sample3). Occurrence is
#     defined as having allele frequencies greater than [0, 0] ([REF, ALT]).
#     """
#
#     de_novo = collection.find(
#                     {
#                         "$and":
#                                 [
#                                     {_SAMPLE_ID_SELECTOR: sample1},
#                                     {
#                                         "$or":
#                                             [
#                                                 {_SAMPLE_ID_SELECTOR: {"$ne": sample2}},
#                                                 {_SAMPLE_ID_SELECTOR: {"$ne": sample3}}
#                                             ]
#                                     }
#                                 ]
#                     })
#     de_novo = list(de_novo)
#     logging.info('Variants found that match de novo criteria: {}'.format(len(de_novo)))
#     return list(de_novo)

