# Welcome to VAPr's Documentation


VAPr (Variant Annotation and PRioritization package) package is aimed at providing a way of retrieving variant information 
 using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) and [myvariant.info](http://myvariant.info/). 
 In particular, it is suited for bioinformaticians interested in aggregating variant information into a single NoSQL 
 database (MongoDB solely at the moment).  
 
The package provides a high-level python API to perform batch annotation jobs [efficiently](link to parallelism), alongside 
 with a powerful set of [querying and filtering protocols](FILTERS). 
 
This documentation is aimed at exposing the main methods available, and how the work from an user perspective. To understand
more deeply what happens under the hood, we advise the user reading the source code and the docstrings provided. There 
is quite a bit of information there.

## Authors

* **Carlo Mazzaferro** (cmazzafe@ucsd.edu)
* **Kathleen Fisch, Ph.D** (kfisch@ucsd.edu)
* **Amanda Birmingham, Ph.D** 
* **Guorong Xu, Ph.D** 

## License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details

<a id='toc'></a>
## Table of contents
1. [Getting started](#getting_started)  
  1.1. [Pre-requisites and setup](#setup)
2. [Quick-start](#examples)
3. [Supplemental Information](#SI)  
  3.1. [Workflow Overview](#workflow)    
  3.2. [Tips on usage](#usage_tips)    
    &nbsp;&nbsp;&nbsp;&nbsp;3.2.1 [Required Arguments](#required)  
    &nbsp;&nbsp;&nbsp;&nbsp;3.2.2 [Optional Arguments](#optional)  
  3.3. [Core Methods](#core)  
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.1 [Annovar Methods](#anno)   
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.2 [Parallel Annotation](#parallel)  
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.3 [Filtering](#filt)  
  3.3. [Filtering](#filt)   
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.1 [Available Filters](#filts)   
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.2 [Parallel Annotation](#parallel)  
    &nbsp;&nbsp;&nbsp;&nbsp;3.3.3 [Filtering](#filt)  
  3.4. [File Creation](#fwrite)   
    &nbsp;&nbsp;&nbsp;&nbsp;3.4.1 [Available Writers](#required)   
    &nbsp;&nbsp;&nbsp;&nbsp;3.4.2 [Extra](#node_specific)   


[Table of contents](#toc)
<a id='getting_started'></a>
## Getting started

These instructions will get you a copy of the package up and running on your local machine, and will enable you to run annotation jobs 
 on any number of vcf files while storing the data in MongoDB. See the [workflow](#workflow) 

<a id='setup'></a>
### Prerequisites

- MongoDB Community Edition. [Instsallation instructions](https://goo.gl/TpBkcb)
- Python (2.7 and 3.5 currently supported and tested)
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
    
    
#### VAPr

VAPr is available from PyPi.  Once the above requirements have been installed, VAPr itself can be installed by just running:

    pip install VAPr


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

[Table of contents](#toc)
<a id='examples'></a>
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
from VAPr.base import AnnotationProject

# Start by specifying the project information
input_dir = '/path/to/vcf/dir'
output_dir = '/path/to/desired/output'
annovar_path = '/path/to/annovar'  
project_data = {'db_name': 'VariantDatabase',
                'collection_name': 'VarCollection'}  # Database and Collection names (optional)


Project = AnnotationProject(input_dir,
                            output_dir,
                            annovar_path,
                            project_data)
 ```
 
This will allow you to use any of [core methods](#core) in the package. They are accessible as methods using dot notation.
That is, if you'd like to use the method `run_annovar()`, you'll do so by simply running:

`Project.run_annovar()`


<a id='required'></a>
### Required Arguments
The first four arguments are required to run the full annotation pipeline. These are:

- `input_dir`: the path to the **directory** where all vcf files live. They may also be inside subdirectories, and VAPr will
 find them.
- `output_dir`: the path to the **directory** where the annotated csv files will be written to. It will be used in two 
 different instances: writing the file outputs from Annovar, and writing the file outputs from VAPr, in case these are 
 needed
- `annovar_path`: path to the Annovar download folder. Once the directory is downloaded from the Annovar website, unzip and place 
 it somewhere that has a good amount of free disk memory. The Annovar databases can take up to ~30GB. 
- `project_data`: a dictionary with keys `db_name` and `collection_name` and values strings of a name of your choosing.

<a id='optional'></a>
### Optional Arguments

- `design_file`: if you are running a job on multiple samples, and you'd like to include extra data about those samples
 in the annotation job, you can add those to the design file. The format for a typical design file can be found 
 [here](design_file). In particular, each line in the design file must refer to either a single vcf file (most common), 
 or a directory containing vcf files from that sample (less common: this would happen if, for instance, your vcf files
 are split by chromosome and are contained in a directory). The design file's requirement is simply having a column
 named 'Sample_Names', and all the other columns are optional. 

- `build_ver`: one of `hg19`, `hg18` or `hg38`. It specifies the genomic build version to be used for annotation
 
[Table of contents](#toc)
<a id='core'></a>
## Core Methods
Once the `AnnotationProject` instance has been initiated, you can call a variety of methods from it. In particular, 
the API lets you call any core method for the annotation part. These include:
 - Annovar methods: `download_dbs` and `run_annovar` and `update_dbs`
 - Annotation and storing methods: `parallel_annotation_and_saving` and `quick_annotate_and_save`

The differences and nuances of each will be discussed next.

<a id='anno'></a>
### Annovar Methods

#### `download_dbs`
`download_dbs()`: this function downloads the databases required to run Annovar to the `.../annovar/humandb/` directory. 
It will download the databases according to the genome version specified. If your databases are out-of-date, re-running
this command will download the latest version of them.

**Args**: 

_Required_: 
  - None

_Optional_: 
  - None

#### `run_annovar`

`run_annovar()`: will spawn Annovar jobs for the files found in the specified directory. This will run Annovar's command
line script on every file and generate output files in the specified `output_dir`. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `batch_jobs`: An integer value that specifies the number of subprocesses to be spawned. When you are working on many,
  many vcf files, running Annovar processes in parallel may improve compute time significantly. Each annotation job 
  occurs as a separate subprocess, so spawning multiple subprocesses makes the annotation parallel in a truly GIL-free 
  fashion. Use always a reasonable number (`number of CPU cores - 1` is usually a good number). Default: 10


<a id='parallel'></a>
### Parallel Annotation

#### `parallel_annotation_and_saving`

`parallel_annotation_and_saving()`: this requires running Annovar beforehand, and will kick-start the main functionality 
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
  - `csv_only`: Determines if only csv data from Annovar will be parsed to MongoDB. This will skip annotation with
  MyVariant.info and will result in much faster annotation. It takes a Boolean True of False. Default: False


#### `quick_annotate_and_save` (not recommended)

`quick_annotate_and_save()`: this can run without having ran Annovar beforehand. It will grab the variant names from the
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

Five different pre-made filters that allow for the retrieval of specific variants have been implemented. Refer to the
[README.md](link to readme/filters) file for more more information about the fields and thresholds used.

### Usage
 In order to use the filters, proceed as follows:
 
 ```python
from VAPr.queries import Filters

db_name = 'VariantDatabase'
collection_name = 'CollectionName'
filter_collection = Filters(db_name, collection_name)
```

This will allow you to use any of the pre-made filters in the package. They are accessible as methods using dot notation.
That is, if you'd like to use the `rare_cancer_variant` filter, you simply need to call:

`rare_vars = filter_collcetion.rare_cancer_variant()`

This will return a list of dictionaries, where each dictionary is contains data about a variant. 

<a id='filts'></a>
### Available filters

#### Filter #1: `rare_cancer_variant`

`rare_cancer_variant()`: this will retrieve all the variants in your collection matching the thresholds specified in the
README.md file. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `samples`: A list of strings specifying the sample names from which you'd like to extract your variants. If this is
  not used, all variants are queried (that is, variants from all sample sin your collection). Default: `None`
  
 
#### Filter #2: `rare_disease_variants`

`rare_disease_variants()`: this will retrieve all the variants in your collection matching the thresholds specified in the
README.md file. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `samples`: A list of strings specifying the sample names from which you'd like to extract your variants. If this is
  not used, all variants are queried (that is, variants from all sample sin your collection). Default: `None`

#### Filter #3: `rare_high_impact_variants`

`rare_high_impact_variants()`: this will retrieve all the variants in your collection matching the thresholds specified in the
README.md file. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `samples`: A list of strings specifying the sample names from which you'd like to extract your variants. If this is
  not used, all variants are queried (that is, variants from all sample sin your collection). Default: `None`

#### Filter #4: `deleterious_compound_heterozygote_variants`

`deleterious_compound_heterozygote_variants()`: this will retrieve all the variants in your collection matching the thresholds 
specified in the README.md file. 

**Args**: 

_Required_: 
  - None

_Optional_: 
  - `samples`: A list of strings specifying the sample names from which you'd like to extract your variants. If this is
  not used, all variants are queried (that is, variants from all sample sin your collection). Default: `None`


#### Filter #5: `de_novo_variants`

`de_novo_variants()`: this will retrieve all the variants in your collection matching the thresholds 
specified in the README.md file. 

**Args**: 

_Required_:
  - `sample1`: first sample name as string. 
  - `sample2`: second sample name as string. 
  - `sample3`: third sample name as string. De novo variants will be looked for in the third sample, i.e. variants that
  occur in the third sample but not in the other two.
  
_Optional_: 
  - `multisample`. Takes a Boolean `True` or `False` value. Set to true if the query is to be performed on
   samples coming from multi-sample vcf files. Default: `False`


#### Create your own filter

As long as you have a MongoDB instance running, filtering can be performed through pymongo as shown by the code below. 
If a list is intended to be created from it, simply add: `filter2 = list(filter2)`

If you'd like to customize your filters, a good idea would be to look at the available fields to be filtered. Looking at 
the myvariant.info [documentation](http://docs.myvariant.info/en/latest/doc/data.html), you can see what are all the 
fields available and can be used for filtering.

```python
from pymongo import MongoClient

client = MongoClient()
db = client.My_Variant_Database
collection = db.ANNOVAR_MyVariant_chunks

filter2 = collection.find({ "$and": [
                                 {"$or": [{"ThousandGenomeAll": {"$lt": 0.05}}, {"ThousandGenomeAll": {"$exists": False}}]},
                                 {"$or": [{"ESP6500siv2_all": { "$lt": 0.05}}, {"ESP6500siv2_all": { "$exists": False}}]},
                                 {"$or": [{"Func_knownGene": "exonic"}, {"Func_knownGene": "splicing"}]},
                                 {"ExonicFunc_knownGene": {"$ne": "synonymous SNV"}},
                                 {"Genotype_Call.DP": {"$gte": 10}},
                                 {"cosmic70": { "$exists": True}}
                         ]})
```

<a id='fwrite'></a>
### File Creation

```python

from VAPr.queries import Filters
from VAPr.writes import Writer

db_name = 'VariantDatabase'
collection_name = 'VariantCollection'
filepath = '.../data/files'

# Create output files (if needed): specify name of files and path 
rare_cancer_variants_csv = filepath + "/tumor_rna_rare_cancer_vars_csv.csv"

# Apply filter
filter_collection = Filters(db_name, collection_name)
rare_cancer_variants = filter_collection.rare_cancer_variants()

# Create writer object for filtered lists
my_writer = Writer(db_name, collection_name)

#cancer variants filtered files
my_writer.generate_annotated_csv(rare_cancer_variants, rare_cancer_variants_csv)
```


## Output Files

### Create unfiltered annotated vcf and csv files 
This output file will contains all annotation data. This may be useful for researchers interested in obtaining a full description of their files.

```python
# Create output files (if needed): specify name of files and path 
out_unfiltered_vcf_file = filepath + "/out_unfiltered_rnaseq_vcf.vcf"
out_unfiltered_csv_file = filepath + "/out_unfiltered_rnaseq_csv.csv"
input_vcf_compressed = filepath + '/test_vcf/Tumor_RNAseq_variants.vcf.gz'

# Create writer object
# db and collection name must be specified since no list is given. The entire collection will be queried.
my_writer = create_output_files.FileWriter(db_name, collection_name)

# Write collection to csv and vcf
# The in_vcf_file must be the .vcf.gz file and it needs to have an associated .tbi file.

my_writer.generate_unfiltered_annotated_csv(out_unfiltered_csv_file)
my_writer.generate_unfiltered_annotated_vcf(input_vcf_compressed, out_unfiltered_vcf_file)
```

### Create filtered annotated vcf and csv files

Further, the package allows the user to parse these variants into an annotated csv or vcf file. 
If needed, annotated, unfiltered vcf and csv files can also be created. They will have the same length 
(number of variants) as the original files, but will contain much more complete annotation data coming from MyVariant.info
 and ANNOVAR databases. 

To create a csv file, just the filtered output is needed. To create an annotated vcf file, a tab indexed file (.tbi) 
file is needed (see comments in  section Create unfiltered annotated vcf and csv files at the end of this page). This 
can be created using tabix.  

First, the file needs to be compressed:

From the command line, running:

`bgzip -c input_file.vcf > input_file.vcf.gz`

returns `input_vcf_file.vcf.gz`

and running 

`tabix input_vcf_file.vcf.gz`

will return: `input_vcf_file.vcf.gz.tbi`

**Citations**

* Wang K, Li M, Hakonarson H. [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research](http://nar.oxfordjournals.org/content/38/16/e164) , 38:e164, 2010
* Xin J, Mark A, Afrasiabi C, Tsueng G, Juchler M, Gopal N, Stupp GS, Putman TE, Ainscough BJ, Griffith OL, Torkamani A, Whetzel PL, Mungall CJ, Mooney SD, Su AI, Wu C (2016) [High-performance web services for querying gene and variant annotation.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9) Genome Biology 17(1):1-7

