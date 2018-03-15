.. -*- mode: rst -*-

Introduction
============

This package is aimed at providing a way of retrieving variant information using ANNOVAR and myvariant.info.
In particular, it is suited for bioinformaticians interested in aggregating variant information into a single
NoSQL database (MongoDB solely at the moment).


Installation
------------

VAPr
~~~~

VAPr is compatible with Python 2.7 or later, but it is preferred to use Python 3.5 or later to take full advantage of all functionality.
The simplest way to install VAPr is from PyPI_ with pip_, Python's preferred package installer.

.. code-block:: bash

    $ pip install VAPr

.. NOTE:: Jupyter, Pandas, and other ancillary libraries are not installed with VAPr and must be installed separately. These can be conveniently install using `Anaconda <https://conda.io/docs/user-guide/install/download.html>`_:


.. code-block:: bash

    $ conda install python=3 pandas mongodb pymongo jupyter notebook


MongoDB
~~~~~~~

Further, MongoDB needs to installed and running. See the `official instructions <https://docs.mongodb.com/manual/installation>`_
for details on how to install it on your platform. Alternatively, using Docker (possibly the easiest way):


.. code-block:: bash

    $ docker run --name vapr -it -p 27017:27017 mongo

This may take a few seconds since the mongodb image needs to be downloaded from a repository. The process
should produce logs that tells you that mongodb is running inside the container. Feel free to leave this
terminal window open and proceed with the rest of the installation.

.. _PyPI: https://pypi.python.org/pypi/yellowbrick
.. _pip: https://docs.python.org/3/installing/


Tabix, ANNOVAR, BCFtools
~~~~~~~~~~~~~~~~~~~~~~~~

**BCFtools** will be used for VCF file merging between samples. To download and install:

.. code-block:: bash
    $ wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2
    $ tar -vxjf bcftools-1.6.tar.bz2
    $ cd bcftools-1.6 && make && make install
    $ export PATH=/where/to/install/bin:$PATH


**Tabix** and **bgzip** binaries are available through the HTSlib project:

.. code-block:: bash
    $ wget https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2
    $ tar -vxjf htslib-1.6.tar.bz2
    $ cd htslib-1.6 && make && make install
    $ export PATH=/where/to/install/bin:$PATH


Refer `here <https://github.com/samtools/htslib/blob/develop/INSTALL>`_ for installation debugging.


**ANNOVAR**

(It is possible to proceed without installing ANNOVAR. Variants will only be annotated with MyVariant.info. In that case,
users can skip the next steps and go straight to the section Known Variant Annotation and Storage)

Users who wish to annotate novel variants will also need to have a local installation of the popular command-line
software ANNOVAR(1), which VAPr wraps with a Python interface. If you use ANNOVAR's functionality through VAPr, please
remember to cite the ANNOVAR publication (see #1 in Citations)!

The base ANNOVAR program must be installed by each user individually, since its license agreement does not permit
redistribution. Please visit the ANNOVAR download form here, ensure that you meet the requirements for a free license,
and fill out the required form. You will then receive an email providing a link to the latest ANNOVAR release file.
Download this file (which will usually have a name like annovar.latest.tar.gz) and place it in the location on your
machine in which you would like the ANNOVAR program and its data to be installed--the entire disk size of the databases
will be around 25 GB, so make sure you have such space available!


Annotation Quickstart
---------------------
An annotation project can be started by providing the API with a small set of information and then running the core
methods provided to spawn annotation jobs. This is done in the following manner:


.. code-block:: python

    # Import core module
    from VAPr import vapr_core
    import os

    # Start by specifying the project information
    IN_PATH = "/path/to/vcf"
    OUT_PATH = "/path/to/out"
    ANNOVAR_PATH = "/path/to/annovar"
    MONGODB = 'VariantDatabase'
    COLLECTION = 'Cancer'

    annotator = vapr_core.VaprAnnotator(input_dir=IN_PATH,
                                       output_dir=OUT_PATH,
                                       mongo_db_name=MONGODB,
                                       mongo_collection_name=COLLECTION,
                                       build_ver='hg19',
                                       vcfs_gzipped=False,
                                       annovar_install_path=ANNOVAR_PATH)

    annotator.download_databases()
    annotator.annotatee()

This will download the required databases from ANNOVAR for annotation and will kickstart the annotation
process, storing the variants in MongoDB.

Filtering Variants
------------------
Four different pre-made filters that allow for the retrieval of specific variants have been implemented. These allow
the user to query in an easy and efficient manner variants of interest

1. Rare Deleterious Variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* criteria 1: 1000 Genomes (ALL) allele frequency (Annovar) < 0.05 or info not available
* criteria 2: ESP6500 allele frequency (MyVariant.info - CADD) < 0.05 or info not available
* criteria 3: cosmic70 (MyVariant.info) information is present
* criteria 4: Func_knownGene (Annovar) is exonic, splicing, or both
* criteria 5: ExonicFunc_knownGene (Annovar) is not "synonymous SNV"


2. Known Disease Variants
~~~~~~~~~~~~~~~~~~~~~~~~~

* criteria: cosmic70 (MyVariant.info) information is present or ClinVar data is present and clinical significance is
not Benign or Likely Benign

3. Deleterious Compound Heterozygous Variants
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

* criteria 1: genotype_subclass_by_class (VAPr) is compound heterozygous
* criteria 2: CADD phred score (MyVariant.info - CADD) > 10

4. De novo Variants
~~~~~~~~~~~~~~~~~~~

* criteria 1: Variant present in proband
* criteria 2: Variant not present in either ancestor-1 or ancestor-2


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


Tutorials
---------

Jupyter Notebooks
