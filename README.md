## VAPr 
### Variant Annotation and Prioritization package

This package is aimed at providing a way of retrieving variant information using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) and [myvariant.info](http://myvariant.info/). In particular, it is suited for bioinformaticians interested in aggregating variant information into a single NoSQL database (MongoDB solely at the moment). 


### Requirements

- [MongoDB installation](https://docs.mongodb.com/getting-started/python/)
- Python 3 and pip3 for installation 
- [Annovar](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)

### Background

VAPr was developed to simplify the steps required to get mutation data from a VCF file to a downstream analysis process.
 A query system was implemented allowing users to quickly slice the genomic variant (GV) data and select variants 
 according to their characteristics, allowing researchers to focus their analysis only on the subset of data that contains meaningful 
 information. Further, this query system allows the user to select the format in which the data can be retrieved. 
 Most notably, CSV or VCF files can be retrieved from the database, allowing any researcher to quickly filter variants 
 and retrieve them in commonly used formats. 
The package can also be installed and used without having to download ANNOVAR. In that case, variant data can be 
retrieved solely by MyVariant.info and rapidly parsed to the MongoDB instance. 

#### Notes


ANNOVAR, alongside with some of their data sets, should to be installed. We can not provide it with this package since
 downloading ANNOVAR requires the user to register on [their website](http://www.openbioinformatics.org/annovar/annovar_download_form.php).
 The data sets required to benefit from the full functionalities of the package are the following:

- knownGene
- tfbsConsSites
- cytoBand
- genomicSuperDups
- esp6500siv2_all
- 1000g2015aug_all
- popfreq_all
- clinvar_20140929
- cosmic70
- nci60

A wrapper is provided that will download automatically these databases. See the complete guide for how to proceed.

More on how to download and troubleshoot: 
ANNOVAR standard script (this is run automatically, as long as the required databases, as well as ANNOVAR, are installed.

After these steps have been taken, you can proceed and install VAPr as follows:

```pip install VAPr```

Alternatively, clone repository on your desktop. Unzip archive and run, on your terminal:


```
cd VAPr
python setup.py build
python setup.py install
```

Sample code for all the methods supplied by the package is provided. 

### Quick-Startup Guide
See this [iPython notebook](https://github.com/ucsd-ccbb/VAPr/blob/master/Quick%20Start-Up%20Guide.ipynb) for a complete tour of the functionalities and installation procedures. 


### Workflow Overview
![Workflow](https://github.com/ucsd-ccbb/VAPr/blob/master/pipeline.png)

### Available Methods
The package offers 2 different ways to obtain variant data. One requires annovar, while the other is based solely on the 
use of MyVariant.info service. The latter can be used without having to worry about downloading and installing annovar 
databases, but it tends to return partial or no information for some variants. 
Calling the methods is extremely straightforward since the syntax is the same. The main difference is in the arguments 
passed to the class VariantParsing.

### Data Models
A sample entry in the Mongo Database will look like [this](https://github.com/ucsd-ccbb/VAPr/blob/master/sample_variant_document). 
The variaty of data that can be retrieved from the sources results from the richness of databases that can be accessed 
through MyVariant.info. However, not every variant will have such data readily available. In some cases, the data will 
be restricted to what can be inferred from the vcf file and the annotation carried out with ANNOVAR. In that case, the 
entries that will be found in the document will be the following:

    {
      "_id": ObjectId("5887c6d3bc644e51c028971a"),
      "chr": "chr19",
      "cytoband": {
        "Region": "13",
        "Sub_Band": "41",
        "Chromosome": 19,
        "Band": "q",
        "Name": "19q13.41"
            },
            
      "genotype": {
        "alleles": [
            0,
            2
            ],
        "genotype_lieklihoods": [
              89.0,
              6.0,
              0.0
                ],
        "filter_passing_reads_count": 2,
        "genotype": "1/1"
          },
          
      "end": 51447065,
      "vcf": {
        "alt": "G",
        "position": "25194768",
        "ref": "A"
            }
        
      "hgvs_id": "chr20:g.25194768A>G",
      "esp6500siv2_all": "0.71"
    }


### Filtering 

Three different pre-made filters that allow for the retrieval of specific variants have been implemented. The filters 
are in the form of MongoDB queries, and are designed to provide the user with a set of relevant variants. In case the 
user would like to define its own querying, a template is provided. 
The output of the queries is a list of dictionaries (JSON documents), where each dictionary contains data relative to 
one variant. 
The data can be retrieved by querying any field that is contained in MyVariant.info or ANNOVAR annotation databases. 
Since there are more than 700 possible keys that variant can have data associated with, it is out of the scope of this
document to list them all. They are however fully available online. To see a list of the full-fledged possible query 
fields can be found at the bottom of [this page](http://docs.myvariant.info/en/latest/doc/data.html). If the queries are
done instead on the ANNOVAR data fields, them the possible keys to be queried upon are the following:

- `chr`   
- `cytoband.Region` 
- `cytoband.Sub_Band`
- `cytoband.Chromosome`
- `cytoband.Band`
- `cytoband.Name` 
- `genotype.alleles`
- `genotype.genotype_lieklihoods`
- `genotype.filter_passing_reads_count`
- `genotype.genotype`
- `start`
- `end`
- `vcf.alt`
- `vcf.position`
- `vcf.ref`
- `hgvs_id`
- `esp6500siv2_all`

These are entries that are guaranteed to be present for every variant. Their typical values can be retrieved from the 
dictionary displayed above, which can serve as a guide to build effective queries that target the desired values. 


### Filter #1: specifying cancer-specifc rare variants

- filter 1: ThousandGenomeAll < 0.05 or info not available
- filter 2: ESP6500siv2_all < 0.05 or info not available
- filter 3: cosmic70 information is present
- filter 4: Func_knownGene is exonic, splicing, or both
- filter 5: ExonicFunc_knownGene is not "synonymous SNV"
- filter 6: Read Depth (DP) > 10

```python
filepath = '.../data/files'
# Create output files (if needed): specify name of files and path 
rare_cancer_variants_csv = filepath + "/tumor_rna_rare_cancer_vars_csv.csv"
rare_cancer_variants_vcf = filepath + "/tumor_rna_rare_cancer_vars_vcf.vcf"
input_vcf_compressed = filepath + '/test_vcf/Tumor_RNAseq_variants.vcf.gz'

# Apply filter
filter_collection = MongoDB_querying.Filters(db_name, collection_name)
rare_cancer_variants = filter_collection.rare_cancer_variant()

# Crete writer object for filtered lists
my_writer = create_output_files.FileWriter(db_name, collection_name)

#cancer variants filtered files
my_writer.generate_annotated_csv(rare_cancer_variants, rare_cancer_variants_csv)
my_writer.generate_annotated_vcf(rare_cancer_variants,input_vcf_compressed, rare_cancer_variants_vcf)
```

### Filter #2: specifying rare disease-specifc (rare) variants

- filter 1: ThousandGenomeAll < 0.05 or info not available
- filter 2: ESP6500siv2_all < 0.05 or info not available
- filter 3: cosmic70 information is present
- filter 4: Func_knownGene is exonic, splicing, or both
- filter 5: ExonicFunc_knownGene is not "synonymous SNV"
- filter 6: Read Depth (DP) > 10
- filter 7: Clinvar data is present 

```python
filter_collection = MongoDB_querying.Filters(db_name, collection_name)
rare_disease_variants = filter_collection.rare_disease_variant()
```


### Filter #3: specifying rare disease-specifc (rare) variants with high impact

- filter 1: ThousandGenomeAll < 0.05 or info not available
- filter 2: ESP6500siv2_all < 0.05 or info not available
- filter 3: cosmic70 information is present
- filter 4: Func_knownGene is exonic, splicing, or both
- filter 5: ExonicFunc_knownGene is not "synonymous SNV"
- filter 6: Read Depth (DP) > 10
- filter 7: Clinvar data is present 
- filter 8: cadd.phred > 10


    filter_collection = MongoDB_querying.Filters(db_name, collection_name)
    rare_high_impact_variants = filter_collection.rare_high_impact_variants()


### Create your own filter

As long as you have a MongoDB instance running, filtering can be perfomed trough pymongo as shown by the code below. If a list is intended to be created from it, simply add: `filter2 = list(filter2)`

If you'd like to customize your filters, a good idea would be to look at the available fields to be filtered. Looking at the myvariant.info [documentation](http://docs.myvariant.info/en/latest/doc/data.html), you can see what are all the fields avaialble and can be used for filtering. 

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



## Output Files

### Create unfiltered annotated vcf and csv files 
This output file will contains all annotation data. This may be useful for researchers interested in obtaining a full description of their files.

```python
# Create output files (if needed): specify name of files and path 
out_unfiltered_vcf_file = filepath + "/out_unfiltered_rnaseq_vcf.vcf"
out_unfiltered_csv_file = filepath + "/out_unfiltered_rnaseq_csv.csv"
input_vcf_compressed = filepath + '/test_vcf/Tumor_RNAseq_variants.vcf.gz'

#Create writer object
#db and collection name must be specified since no list is given. The entire collection will be queried. 
my_writer = create_output_files.FileWriter(db_name, collection_name)

#Write collection to csv and vcf
#The in_vcf_file must be the .vcf.gz file and it needs to have an associated .tbi file. 

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

Wang K, Li M, Hakonarson H. [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research](http://nar.oxfordjournals.org/content/38/16/e164) , 38:e164, 2010
Xin J, Mark A, Afrasiabi C, Tsueng G, Juchler M, Gopal N, Stupp GS, Putman TE, Ainscough BJ, Griffith OL, Torkamani A, Whetzel PL, Mungall CJ, Mooney SD, Su AI, Wu C (2016) [High-performance web services for querying gene and variant annotation.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9) Genome Biology 17(1):1-7