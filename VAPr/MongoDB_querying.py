from pymongo import MongoClient


class Filters(object):

    """
    Class for variant filtering. DB and collection names can be passed directly to the class and will be used in each
    method.

    Description:

        Rare variants filters are aimed at finding variants that occur with low frequency in the general population.
        Further, some of the filters implemented also target variants whose gene function is known. Higher selectivity
        filters such as the rare_high_impact_variant filter also screens variant data function using CADD Phred scores.

    """

    def __init__(self, db_name, collection_name):

        self.collection_name = collection_name
        self.db_name = db_name

    def rare_cancer_variant(self):
        """ Function for retrieving rare cancer variants """

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        filtered = collection.find(
            {
                "$and":
                    [
                        {
                            "$or":
                                [
                                    {"esp6500siv2_all": {"$lt": 0.05}},
                                    {"esp6500siv2_all": {"$exists": False}}
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
                        {"filter_passing_reads_count": {"$gte": 10}},
                        {"cosmic70": {"$exists": True}},
                        {"1000g2015aug_all": {"$lt": 0.1}}
                    ]
            }
        )

        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def rare_disease_variant(self):
        """ Function for retrieving rare disease variants. Includes filter on Clinvar data  """

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)

        filtered = collection.find(
            {
                "$and":
                    [
                        {
                            "$or":
                                [
                                   {"esp6500siv2_all": {"$lt": 0.05}},
                                   {"esp6500siv2_all": {"$exists": False}}
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
                        {"filter_passing_reads_count": {"$gte": 10}},
                        {
                            "$or":
                                [
                                    {"cosmic70": {"$exists": True}},
                                    {"clinvar": {"$exists": True}}
                                ]
                        },
                        {"1000g2015aug_all": {"$lt": 0.1}}
                    ]
            }
        )

        filtered = list(filtered)
        print ('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def rare_high_impact_variants(self):
        """ Function for retrieving rare disease variants. Includes filter on Clinvar as well CADD Phred data  """

        client = MongoClient()

        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)

        filtered = collection.find(
            {
                "$and":
                    [
                        {
                            "$or":
                                [
                                    {"esp6500siv2_all": {"$lt": 0.05}},
                                    {"esp6500siv2_all": {"$exists": False}}
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
                        {"filter_passing_reads_count": {"$gte": 10}},
                        {
                            "$or":
                                [
                                    {"cosmic70": {"$exists": True}},
                                    {"clinvar": {"$exists": True}}
                                ]
                        },
                        {"1000g2015aug_all": {"$lt": 0.1}},
                        {"cadd.phred": {"$gte": 15}}
                    ]
            }
        )

        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def deleterious_compound_heterozygote(self):
        """ Function for retrieving deleterious compound heterozygote variants  """

        client = MongoClient()

        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)

        filtered = collection.find(
            {
                "$and":
                    [
                        {"genotype_subclass_by_class.heterozygous": "compound"},
                        {"cadd.phred": {"$gte": 10}}
                    ]
            }
        )

        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def get_de_novo(self, sample1, sample2, sample3, multisample=True):
        """
        Function for de novo variant analysis. Can be performed on multisample files or or on data coming
        from a collection of files. In the former case, every sample contains the same variants, although they have
        differences in their allele frequency and read values. A de novo variant is defined as a variant that
        occurs only in the specified sample (sample1) and not on the other two (sample2, sample3). Occurrence is
        defined as having allele frequencies greater than [0, 0] ([REF, ALT]).

        For the case in which the variants come from multiple files, de novo variants are the ones that occur only
        once on the sample specified (usually found in one file), and that do not occur on the other samples (found
        in the other files).
        """

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        if multisample:

            filtered_hgvs = collection.distinct("hgvs_id",
                                                {
                                                    "$or": [
                                                            {"$and": [
                                                                {"alleles": [0, 0]},
                                                                {"sample_id": sample2}
                                                            ]
                                                            },
                                                            {"$and": [
                                                                {"alleles": [0, 0]},
                                                                {"sample_id": sample3}
                                                            ]
                                                            },
                                                        ]
                                                }
                                                )

            de_novo = collection.find(
                {
                    "$and":
                        [
                            {"hgvs_id": {"$in": filtered_hgvs}},
                            {"alleles": {"$ne": [0, 0]}},
                            {"sample_id": sample1}
                        ]

                }
            )
        else:
            filtered_hgvs = collection.aggregate(
                [
                    {"$group": {
                        "_id": "$hgvs_id",
                        "samples":
                            {"$addToSet": {
                                "$cond": [  # works with a if ... then ... else syntax

                                    {  # if these are satisfied,
                                        "$or": [
                                            {"$eq": ["$sample_id", sample1]},
                                            {"$eq": ["$sample_id", sample2]},
                                            {"$eq": ["$sample_id", sample3]}
                                        ]
                                    },
                                    # then add sample_id to the set of variants
                                    "$sample_id",
                                    # else, add zero
                                    0
                                ]

                            }
                            }
                    }
                    },
                    {"$redact": {
                        "$cond": {
                            "if": {
                                "$and": [
                                    {"$in": [sample1, "$samples"]},
                                    {"$eq": [{"$size": "$samples"}, 1]}  # Occurs only once, hence a unique HGVS_ID
                                ]
                            },
                            "then": "$$KEEP",
                            "else": "$$PRUNE"
                        }
                    }
                    }
                ]
            )
            as_hgvs_list = [i['_id'] for i in filtered_hgvs]
            de_novo = collection.find({'hgvs_id': {'$in': as_hgvs_list}})

        print('Variants found that match de novo criteria: {}'.format(len(de_novo)))
        return list(de_novo)
