import logging
import os
import shlex
import subprocess
import time

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
                            {"clinvar.rcv.clinical_significance": {"$nin": ["Benign", "Likely benign"]}},
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
                                        {
                                            "$or":
                                                [
                                                    {"samples.sample_id" : {"$ne": sample2}},
                                                    {"samples.sample_id": {"$ne": sample3}}
                                                ]
                                        }
                                    ]
                        })
        de_novo = list(de_novo)
        print('Variants found that match de novo criteria: {}'.format(len(de_novo)))
        return list(de_novo)


def store_annotations_to_db(annotation_dicts_list, db_name, collection_name, client=None, mongod_cmd=None):
    if client is None:
        client = MongoClient(maxPoolSize=None, waitQueueTimeoutMS=200)

    db = getattr(client, db_name)
    collection = getattr(db, collection_name)

    logging.info('Parsing Buffer...')
    if len(annotation_dicts_list) == 0:
        logging.info('Empty list of documents trying to be parsed, skipping and continuing operation...')
        return None
    else:
        try:
            collection.insert_many(annotation_dicts_list, ordered=False)
        except Exception as error:
            if "Connection refused" in str(error):
                if mongod_cmd:
                    logging.info('MongoDB Server seems to be off. Attempting restart...')
                    args = shlex.split(mongod_cmd)
                    subprocess.Popen(args, stdout=subprocess.PIPE)
                else:
                    logging.error('Error connecting to mongodb.' + str(error))

                logging.info('Retrying to parse...')
                time.sleep(4)

                store_annotations_to_db(annotation_dicts_list, collection, client, mongod_cmd)  # Recurse!
            else:
                logging.error('Unrecoverable error: ' + str(error))

    try:
        client.close()
    except:
        pass  # if the client is already closed, just move along


def generate_output_files_by_sample(db_name, collection_name, output_dir):
    # TODO: finish refactor of this method
    raise NotImplementedError("Not refactored yet")

    client = MongoClient()
    db = getattr(client, db_name)
    collection = getattr(db, collection_name)

    # Get all distinct sample_ids
    samples = collection.distinct('sample_id')

    fwriter = Writer(db_name, collection_name)
    filt = Filters(db_name, collection_name)

    # Generate files for each sample
    for sample in samples:
        q = filt.variants_from_sample(sample)
        list_docs = list(q)
        out_path_by_sample = os.path.join(output_dir, 'csv_by_sample')
        if not os.path.exists(out_path_by_sample):
            os.makedirs(out_path_by_sample)
        fname = os.path.join(out_path_by_sample, 'annotated_csv_' + sample + '_all_vars.csv')
        fwriter.generate_annotated_csv(list_docs, fname)
