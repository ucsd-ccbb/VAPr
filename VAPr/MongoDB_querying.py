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

        filtered = collection.find({"$and": [
                                     {"$or": [{"ThousandGenomeAll": {"$lt": 0.05}}, {"ThousandGenomeAll": {"$exists": False}}]},
                                     {"$or": [{"ESP6500siv2_all": { "$lt": 0.05}}, {"ESP6500siv2_all": {"$exists": False}}]},
                                     {"$or": [{"Func_knownGene": "exonic"}, {"Func_knownGene": "splicing"}]},
                                     {"ExonicFunc_knownGene": {"$ne": "synonymous SNV"}},
                                     {"Genotype_Call.DP": {"$gte": 10}},
                                     {"cosmic70": {"$exists": True}}
        ]})

        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def rare_disease_variant(self):
        client = MongoClient()

        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)

        filtered = collection.find({"$and": [
                                   {"$or": [{"ThousandGenomeAll": {"$lt": 0.05}}, {"ThousandGenomeAll": {"$exists": False}}]},
                                   {"$or": [{"ESP6500siv2_all": {"$lt": 0.05}}, {"ESP6500siv2_all": {"$exists": False}}]},
                                   {"$or": [{"Func_knownGene": "exonic"}, {"Func_knownGene": "splicing"}]},
                                   {"ExonicFunc_knownGene": {"$ne": "synonymous SNV"}},
                                   {"Genotype_Call.DP": {"$gte": 10}},
                                   {"cosmic70": {"$exists": True}},
                                   {"clinvar": {"$exists": True}} ##This is the change

       ]})

        filtered = list(filtered)
        print ('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered

    def rare_high_impact_variants(self):
        """Rare high impact (CADD) variant filter"""

        client = MongoClient()

        db = getattr(client, self.db_name)
        collection = getattr(db, self.collection_name)

        filtered = collection.find({"$and": [
                                   {"$or": [{"ThousandGenomeAll": {"$lt": 0.05}}, {"ThousandGenomeAll": {"$exists": False}}]},
                                   {"$or": [{"ESP6500siv2_all": {"$lt": 0.05}}, {"ESP6500siv2_all": {"$exists": False}}]},
                                   {"$or": [{"Func_knownGene": "exonic"}, {"Func_knownGene": "splicing"}]},
                                   {"ExonicFunc_knownGene": {"$ne": "synonymous SNV"}},
                                   {"Genotype_Call.DP": {"$gte": 10}},
                                   {"cadd.phred": {"$gte": 15}}, ##This is the change CADD Phred score >= 15
                                   {"cosmic70": {"$exists": True}} ,
                                   {"clinvar": {"$exists": True}}  ##This is the change

         ]})


        filtered = list(filtered)
        print('Variants found that match rarity criteria: {}'.format(len(filtered)))
        return filtered



