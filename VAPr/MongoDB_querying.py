from pymongo import MongoClient


class Filters(object):

    """
    Class for variant filtering. DB and collection names can be passed directly to the class and will be used in each
    method.
    """

    def __init__(self, db_name, collection_name):

        self.collection_name = collection_name
        self.db_name = db_name

    def rare_cancer_variant(self):
        """
        Function for retrieving rare cancer variants.
        :return:
        """

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)
        #
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
        """Rare high impact (CADD) variant filter"""

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

    def get_de_novo(self, de_novo_sample):

        client = MongoClient()
        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)

        filtered = collection.find(
            {
                "$and":
                    [
                        {"alleles": [0, 0]},
                        {"sample_id": de_novo_sample}
                    ]
            }
        )

        de_novo_vars = []

        for var in filtered:
            de_novo_counter = 0
            others = list(collection.find({'hgvs_id': var['hgvs_id']}))
            for other in others:
                if (other['sample_id'] != de_novo_sample) and (other['alleles'] != [0, 0]):
                    de_novo_counter += 1
                if de_novo_counter == 2:
                    de_novo_vars.append(var)
        print('Variants found that match criteria: {}'.format(len(de_novo_vars)))

        return de_novo_vars

