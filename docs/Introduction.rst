.. -*- mode: rst -*-

Introduction
============

This package is aimed at providing a way of retrieving variant information using ANNOVAR and myvariant.info.
In particular, it is suited for bioinformaticians interested in aggregating variant information into a single
NoSQL database (MongoDB solely at the moment).

.. _install:

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
    dataset = annotator.annotate(num_processes=8)

This will download the required databases from ANNOVAR for annotation and will kickstart the annotation
process, storing the variants in MongoDB.
