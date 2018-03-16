.. _downstream-analysis:

Downstream Analysis
===================

For notes on how to implement these features, refer to the :ref:`tutorial` and the :ref:`API Reference <api-reference>`


Filtering Variants
------------------
Four different pre-made filters that allow for the retrieval of specific variants have been implemented. These allow
the user to query in an easy and efficient manner variants of interest

.. _rare-del-variants:

1. Rare Deleterious Variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* criteria 1: 1000 Genomes (ALL) allele frequency (Annovar) < 0.05 or info not available
* criteria 2: ESP6500 allele frequency (MyVariant.info - CADD) < 0.05 or info not available
* criteria 3: cosmic70 (MyVariant.info) information is present
* criteria 4: Func_knownGene (Annovar) is exonic, splicing, or both
* criteria 5: ExonicFunc_knownGene (Annovar) is not "synonymous SNV"

.. _known-disease:

2. Known Disease Variants
~~~~~~~~~~~~~~~~~~~~~~~~~

* criteria: cosmic70 (MyVariant.info) information is present or ClinVar data is present and clinical significance is
not Benign or Likely Benign

.. _del-compound:

3. Deleterious Compound Heterozygous Variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* criteria 1: genotype_subclass_by_class (VAPr) is compound heterozygous
* criteria 2: CADD phred score (MyVariant.info - CADD) > 10

.. _de-novo:

4. De novo Variants
~~~~~~~~~~~~~~~~~~~

* criteria 1: Variant present in proband
* criteria 2: Variant not present in either ancestor-1 or ancestor-2

.. _custom-filter:

Create your own filter
~~~~~~~~~~~~~~~~~~~~~~

As long as you have a MongoDB instance running and an annotation job ran successfully, filtering can be performed
through `pymongo <https://api.mongodb.com/python/current/>`_ as shown by the code below.
Running the query will return a :code:`cursor` object, which can be iterated upon.

If instead a list is intended to be created from it, simply add: filter2 = list(filter2)

.. WARNING:: If the number of variants in the database is large and the filtering is not set up correctly,
returning a list will be probably crash your computer since lists are kept in memory. Iterating over the cursor
object perform `lazy evaluations` (i.e., one item is returned at a time instead of in bulk) which are much more memory
efficient.

Further, if you'd like to customize your filters, a good idea would be to look at the available fields to be filtered.
Looking at the myvariant.info documentation, you can see what are all the fields available and can be used for filtering.


.. code-block:: python

    from pymongo import MongoClient

    client = MongoClient()
    db = getattr(client, mongodb_name)
    collection = getattr(db, mongo_collection_name)

    filtered = collection.find({"$and": [
                                       {"$or": [{"func_knowngene": "exonic"},
                                                {"func_knowngene": "splicing"}]},
                                       {"cosmic70": {"$exists": True}},
                                       {"1000g2015aug_all": {"$lt": 0.05}}
                             ]})

    # filtered = list(filtered) Uncomment this if you'd like to return them as a list
    for var in filtered:
        print(var)

Output Files
------------
Although iterating over variants can be interesting for cursory analyses, we provide functionality to retrieve as well
csv files for downstream analysis. A few options are available:

Unfiltered Variants
~~~~~~~~~~~~~~~~~~~

:code:`write_unfiltered_annotated_csv(out_file_path)`

* All variants will be written to a CSV file.


Filtered Variants
~~~~~~~~~~~~~~~~~

:code:`write_filtered_annotated_csv(variant_list, out_file_path)`

* A list of filtered variants will be written to a CSV file.


Unfiltered Variants VCF
~~~~~~~~~~~~~~~~~~~~~~~

:code:`write_unfiltered_annotated_vcf(vcf_out_path)`

* All variants will be written to a VCF file.


Write Options #4: Filtered Variants VCF

:code:`write_filtered_annotated_csv(variant_list, vcf_out_path)`

* A List of filtered variants will be written to a VCF file.
