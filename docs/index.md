# Welcome to VAPr's Documentation


VAPr (Variant Annotation and PRioritization package) package is aimed at providing a way of retrieving variant information 
 using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) and [myvariant.info](http://myvariant.info/). 
 In particular, it is suited for bioinformaticians interested in aggregating variant information into a single NoSQL 
 database (MongoDB at the moment).  
 
The package provides a high-level python API to perform batch annotation jobs [efficiently](link to parallelism), alongside 
 a powerful set of [querying and filtering protocols](FILTERS). 
 
This documentation is aimed at exposing the main methods available, and how they work from a user perspective. To better understand
what happens under the hood, we advise the user reading the source code and the docstrings provided. There 
is quite a bit of information there.

## Authors

* **Carlo Mazzaferro** (cmazzafe@ucsd.edu)
* **Adam Mark, M.S.** (a1mark@ucsd.edu)
* **Amanda Birmingham, B.A.** (abirmingham@ucsd.edu) 
* **Guorong Xu, Ph.D** 
* **Kathleen Fisch, Ph.D** (kfisch@ucsd.edu)

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

<a id='toc'></a>
## Table of contents
1. [Getting started](#getting_started)  
  1.1. [Pre-requisites and setup](#setup)
2. [Quick-start](#examples)
3. [Supplemental Information](#SI)  
  3.1. [Workflow Overview](#workflow)    
  3.2. [VaprAnnotator - Tips on usage](#usage_tips)    
    &nbsp;&nbsp;&nbsp;&nbsp;3.2.1 [Arguments](#required)  
  3.3. [Core Methods](#core)  
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.1 [Annovar Databases](#anno)   
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.2 [Annotation](#annotation)  
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.3 [Filtering](#filt)  
  3.3. [Filtering](#filt)   
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.1 [Available Filters](#filts)   
  3.4. [File Creation](#write)   


[Table of contents](#toc)
<a id='getting_started'></a>
## Getting started

These instructions will get you a copy of the package up and running on your local machine, and will enable you to run annotation jobs 
 on any number of vcf files while storing the data in MongoDB. See the [workflow](#workflow) 

<a id='setup'></a>
### Prerequisites

- MongoDB Community Edition. [Instsallation instructions](https://goo.gl/TpBkcb)
- Python (2.7 and 3.5 currently supported and tested)
- [BCFtools](http://www.htslib.org/download/)
- [Tabix](http://www.htslib.org/download/)
- [Annovar scripts](http://annovar.openbioinformatics.org/en/latest/user-guide/download/) (optional)

#### Python 3 and MongoDB

VAPr is written in Python and stores variant annotations in NoSQL database, using a locally-installed instance of MongoDB.
  If you wish to run the Jupyter notebooks provided with VAPr locally, it is also necessary to have a local Jupyter 
  notebook server installed. To set up these requirements in your machine's default environment using `conda`, run the 
  following command:

    conda install python=3 pandas mongodb pymongo jupyter notebook

MongoDB also needs a location to store its data, so create a directory for this in the location of your choice, e.g.:

    mkdir -p /Temp/MongoDbData
    
(or, on Windows systems, `md /Temp/MongoDbData`). 

Then start the MongoDb with the command

    mongod --dbpath /Temp/MongoDbData

    
To check if mongodb is currently running, run:

    service mongod status
    
which should return    
    
    mongod (pid 591) is running...
    
519 is the process number so it may differ every time it is run. In case it is not running, the command will return:

    mongod dead but subsys locked

### BCFtools

BCFtools will be used for VCF file merging between samples. To download and install:
    
    wget https://github.com/samtools/bcftools/releases/download/1.6/bcftools-1.6.tar.bz2
    tar -vxjf bcftools-1.6.tar.bz2
    cd bcftools-1.6
    make
    make install
    export PATH=/where/to/install/bin:$PATH

Refer [here](https://github.com/samtools/bcftools/blob/develop/INSTALL) for installation debugging.

### Tabix

Tabix and bgzip binaries are available through the HTSlib project: 

    wget https://github.com/samtools/htslib/releases/download/1.6/htslib-1.6.tar.bz2
    tar -vxjf htslib-1.6.tar.bz2
    cd bcftools-1.6
    make
    make install
    export PATH=/where/to/install/bin:$PATH

Refer [here](https://github.com/samtools/htslib/blob/develop/INSTALL) for installation debugging.

### ANNOVAR
(It is possible to proceed without installing ANNOVAR. In that case, the variants that will be annotated and sent to 
Mongo are the ones found in MyVariant.info. In that case, users can skip the next steps and go straight to the section 
**Known Variant Annotation and Storage**)

Users who wish to annotate novel variants will also need to have a local installation of the popular command-line software 
ANNOVAR([1](#Citations)), which VAPr wraps with a Python interface.  If you use ANNOVAR's functionality through VAPr, 
please remember to cite the ANNOVAR publication (see #1 in [Citations](#Citations))!

The base ANNOVAR program must be installed by each user individually, since its license agreement does not permit 
redistribution.  Please visit the ANNOVAR download form [here](http://www.openbioinformatics.org/annovar/annovar_download_form.php), 
ensure that you meet the requirements for a free license, and fill out the required form. You will then receive an email 
providing a link to the latest ANNOVAR release file. Download this file (which will usually have a name like 
`annovar.latest.tar.gz`) and place it in the location on your machine in which you would like the ANNOVAR program and 
its data to be installed--the entire disk size of the databases will be around 25 GB, so make sure you have such space 
available!  

#### VAPr

VAPr is available from PyPi.  Once the above requirements have been installed, VAPr itself can be installed by just running:

    pip install VAPr

## Quick-Start
See this [jupyter notebook](https://github.com/ucsd-ccbb/VAPr/blob/master/VAPr%20Quick-Start%20Guide.ipynb) to create your first annotation job

[Table of contents](#toc)
<a id='SI'></a>
## Supplemental Information


<a id='workflow'></a>
### Workflow Overview
![Workflow](https://raw.githubusercontent.com/ucsd-ccbb/VAPr/master/pipeline.png)

<a id='usage_tips'></a>
### Tips on Usage
An annotation project can be started by providing the API with a small set of information and then running the 
 [core methods](#core) provided to spawn annotation jobs. This is done in the following manner:
 
 ```python
 
# Import core module
from VAPr import vapr_core

# Start by specifying the project information
input_directory = "/path/to/my/vcf"
output_director = "/path/tp/desired/output"
annovar_directory = '/path/to/annovar'
mongodb = 'VariantDatabase'
collection = 'Cancer'

annotator = vapr_core.VaprAnnotator(input_directory,
                                    output_director,
                                    mongodb,
                                    collection,
                                    annovar_install_path=annovar_directory,
                                    build_ver='hg19')
 ```
 
This will allow you to use any of [core methods](#core) in the package. 

`dataset = annotator.annotate(num_processes=8)`


<a id='required'></a>
### Required Arguments
The first four arguments are required to run the full annotation pipeline. These are:

- `input_dir`: the path to the **directory** where all vcf files live. They may also be inside subdirectories, and VAPr will
 find them.
- `output_dir`: the path to the **directory** where the annotated csv files will be written to. It will be used in two 
 different instances: writing the file outputs from Annovar, and writing the file outputs from VAPr, in case these are 
 needed

- `build_ver `: Human genome build. VAPr currently supports the two human genome builds, `hg19`, `hg38`.

- `mongo_db_name`: Database used for variant storage as well as the output of annovar i.e. 'VariantDatabase'.

- `mongo_collection_name`: Collection name for this analysis i.e. 'cancer_analysis_01012018'.

<a id='optional'></a>
### Optional Arguments

- `annovar_install_path`: Path location of annovar installtion. **NOTE**: we can't provide the ANNOVAR package as it requires registration to be downloaded. It is, however, freely available [here](http://annovar.openbioinformatics.org/en/latest/user-guide/download/). Download, unzip and place it in whatever directory you'd like. Make sure you have enough space on disk (~15 GB for the datasets used for annotation). It is required that you specify the location to which you downloaded annovar. The folder where annovar lives looks like this:

    ... /annovar/
             annotate_variation.pl
             coding_change.pl
             convert2annovar.pl
             example/    
             humandb/
             retrieve_seq_from_fasta.pl
             table_annovar.pl
             variants_reduction.pl 

- `vcfs_gzipped`: Boolean. Only files with one vcf extension will be processed. If you are only analyzing one vcf, the file will not be bgzipped. If you are providing a directory or design file with multiple vcf files, they will be bgzipped and merged. If they are already bgzipped, please specify `vcfs_gzipped=True` and the bgzip step will be skipped.

- `design_file`: Path to design file. The purpose of an optional design file is to accommodate VCF files scattered throughout a file system. The design file must be set up as a CSV file with the first field name as "Sample_Names", where the column should be populated with full file paths to each VCF you wish to include in the analysis. We anticipate in the future to be able to accommodate meta-data as successive columns which would be included as sample information in each variant document. A sample design file:

            Sample_Names
            /path/to/file1.vcf
            /path/to/file2.vcf
            /path/to/file3.vcf

 
[Table of contents](#toc)
<a id='core'></a>
## Core Methods
Once the `VaprAnnotator` instance has been initiated, you can call a variety of methods from it. In particular, 
the API lets you call any core method for the annotation part. These include:
 - Annovar methods: `download_annovar_databases` and `_run_annovar_annotation` and `update_dbs`
 - Annotation and storing methods: `annotate` and `annotate_lite`

The differences and nuances of each will be discussed next.

<a id='anno'></a>
### Annovar

#### `download_annovar_databases`
`download_annovar_databases()`: this function downloads the databases required to run Annovar to the `.../annovar/humandb/` directory. 
It will download the databases according to the genome version specified. If your databases are out-of-date, re-running
this command will download the latest version of them. If you currently have the required databases, you may get an error.

**Args**: 

_Required_: 
  - None

_Optional_: 
  - None

<a id='annotation'></a>
### Annotation

#### `annotate`

`annotate()`: this requires running Annovar beforehand, and will kick-start the main functionality
of this package. Namely, it will collect all the variant data from Annovar annotations, combine it with data coming
from MyVariant.info, and parse it to MongoDB, in the database and collection specified in `project_data`.

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `n_processes`: An integer value that specifies the number of parallel processing jobs to be used. Larger vcf files 
  will always benefit from a larger number of parallel processes, but that may not be the case for smaller vcf files. 
  As a rule of thumb, use at most `number of CPU cores - 1`, and for smaller vcf files (less than 50 thousand variants)
  4-5 cores at most. Default: 4.
  - `verbose`: An integer value from 0 to 3 that specifies the verbosity level. Default: 0.

#### `annotate_lite` (not recommended)

`annotate_lite()`: Execution will skip annotating with Annovar. It will grab the HGVS ids from the
vcf files and query the variant data from MyVariant.info. It is subject to the issue of potentially having completely
empty data for some of the variants, and inability to run native VAPr queries on the data. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `n_processes`: An integer value that specifies the number of parallel processing jobs to be used. Use here 
  `number of CPU cores - 1`. 
  

<a id='filt'></a>
## Filtering  Methods

Four different pre-made filters that allow for the retrieval of specific variants have been implemented. Refer to the
[README.md](link to readme/filters) file for more more information about the fields and thresholds used.

### Usage
 In order to use the filters, proceed as follows:
 
 ```python
rare_deleterious_variants = dataset.get_rare_deleterious_variants()
```

This will return a list of dictionaries, where each dictionary is contains variant containing annotations. 

<a id='filts'></a>
### Available filters

#### Filter #1: Rare Deleterious Variants

`get_rare_deleterious_variants()`: this will retrieve all the variants in your collection matching the thresholds specified in the
README.md file. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `sample_names_list`: A list of strings specifying the sample names from which you'd like to extract your variants. If this is
  not used, all variants are queried (that is, variants from all sample sin your collection). Default: `None`
  
 
#### Filter #2: Known Disease Variants

`get_known_disease_variants()`: this will retrieve all the variants in your collection matching the thresholds specified in the
README.md file. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `sample_names_list`: A list of strings specifying the sample names from which you'd like to extract your variants. If this is
  not used, all variants are queried (that is, variants from all sample sin your collection). Default: `None`

#### Filter #3: Deleterious Compound Heterozygous Variants

`get_deleterious_compound_heterozygous_variants()`: this will retrieve all the variants in your collection matching the thresholds specified in the
README.md file. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `sample_names_list`: A list of strings specifying the sample names from which you'd like to extract your variants. If this is
  not used, all variants are queried (that is, variants from all sample sin your collection). Default: `None`

#### Filter #4: De novo Variants

`get_de_novo_variants()`: this will retrieve all the variants in your collection matching the thresholds 
specified in the README.md file. 

**Args**: 

_Required_:
  - `proband`: first sample name as string. 
  - `ancestor1`: second sample name as string. 
  - `ancestor2`: third sample name as string. De novo variants will be looked for in the proband sample, i.e. variants that
  occur in the first sample but not in either ancestor1 or ancestor2.
  

#### Create your own filter

As long as you have a MongoDB instance running, filtering can be performed through pymongo as shown by the code below. 
If a list is intended to be created from it, simply add: `filter2 = list(filter2)`

If you'd like to customize your filters, a good idea would be to look at the available fields to be filtered. Looking at 
the myvariant.info [documentation](http://docs.myvariant.info/en/latest/doc/data.html), you can see what are all the 
fields available and can be used for filtering.

```python
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
filtered = list(filtered)
```

<a id='write'></a>
## Output Files

#### Write Options #1: Unfiltered Variants CSV

`write_unfiltered_annotated_csv()`: All variants will be written to a CSV file.

**Args**: 

_Required_:
  - `output_fp`: Name of output file path

#### Write Options #2: Filtered Variants CSV

`write_filtered_annotated_csv()`: List of filtered variants will be written to a CSV file.

**Args**: 

_Required_:
  - `filtered_variants`: List of filtered variants retrieved by VAPr filters or custom filters.

  - `output_fp`: Name of output file path.

#### Write Options #3: Unfiltered Variants VCF

`write_unfiltered_annotated_vcf()`: All variants will be written to a VCF file.

**Args**: 

_Required_:
  - `vcf_out_path`: Name of output file path

_Optional_:
  - `info_out`: if set to true (Default), will write all annotation data to INFO column, else, it won't.

#### Write Options #4: Filtered Variants VCF

`write_filtered_annotated_csv()`: List of filtered variants will be written to a VCF file.

**Args**: 

_Required_:
  - `filtered_variants`: List of filtered variants retrieved by VAPr filters or custom filters.
  
  - `vcf_out_path`: Name of output file path.

```python

# List of rare deleterious variants
filtered_variants = dataset.get_rare_deleterious_variants()
# Write variants to csv file
dataset.write_filtered_annotated_csv(filtered_variants, output_dir + “/myfile.csv”)

```


**Citations**

* Wang K, Li M, Hakonarson H. [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research](http://nar.oxfordjournals.org/content/38/16/e164) , 38:e164, 2010
* Xin J, Mark A, Afrasiabi C, Tsueng G, Juchler M, Gopal N, Stupp GS, Putman TE, Ainscough BJ, Griffith OL, Torkamani A, Whetzel PL, Mungall CJ, Mooney SD, Su AI, Wu C (2016) [High-performance web services for querying gene and variant annotation.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9) Genome Biology 17(1):1-7

