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

    def variants_from_sample(self, sample_name):

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        filtered = collection.find({'samples.sample_id': sample_name})

        return list(filtered)

    def variants_from_samples(self, sample_list):

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        filtered = collection.find({'samples.sample_id': {'$in': sample_list}})

        return list(filtered)

    def rare_deleterious_variants(self, samples=None):
        """ Function for retrieving rare, deleterious variants """

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        if not samples:
            samples = collection.distinct('samples.sample_id')
        if not isinstance(samples, list):
            samples = [samples]

        filtered = collection.find(
            {
                "$and":
                    [
                        {
                            "$or":
                                [
                                    {"cadd.esp.af": {"$lt": 0.05}},
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
                        {"1000g2015aug_all": {"$lt": 0.05}},
                        {'samples.sample_id': {'$in': samples}}

                    ]
            }
        )

        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def known_disease_variants(self, samples=None):
        """ Function for retrieving known disease variants by presence in Clinvar and Cosmic."""

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        if not samples:
            samples = collection.distinct('samples.sample_id')
        if not isinstance(samples, list):
            samples = [samples]

        filtered = collection.find(
            {
                "$and":
                        [
                            {"clinvar.rcv.accession": {"$exists": True}},
                            {"clinvar.rcv.accession": {"$nin": ["Benign", "Likely benign"]}},
                            {"cosmic.cosmic_id": {"$exists": True}}
                        ]
            }
        )

        filtered = list(filtered)
        print ('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def deleterious_compound_heterozygote_variants(self, samples=None):
        """ Function for retrieving deleterious compound heterozygote variants  """

        client = MongoClient()

        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        if not samples:
            samples = collection.distinct('samples.sample_id')
        if not isinstance(samples, list):
            samples = [samples]

        filtered = collection.find(
            {
                "$and":
                    [
                        {"genotype_subclass_by_class.heterozygous": "compound"},
                        {"cadd.phred": {"$gte": 10}},
                        {'samples.sample_id': {'$in': samples}}
                    ]
            }
        )

        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def de_novo_variants(self, sample1, sample2, sample3):
        """
        Function for de novo variant analysis. Can be performed on multisample files or or on data coming
        from a collection of files. In the former case, every sample contains the same variants, although they have
        differences in their allele frequency and read values. A de novo variant is defined as a variant that
        occurs only in the specified sample (sample1) and not on the other two (sample2, sample3). Occurrence is
        defined as having allele frequencies greater than [0, 0] ([REF, ALT]).
        """

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        de_novo = collection.find(
                        {
                            "$and":
                                    [
                                        {"samples.sample_id": sample1},
                                        {"samples.sample_id" : {"$nin": [sample2, sample3]}}
                                    ]})
        de_novo = list(de_novo)
        print('Variants found that match de novo criteria: {}'.format(len(de_novo)))
        return list(de_novo)
