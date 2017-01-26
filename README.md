##VAPr 
###Variant Annotation and Prioritization package

This package is aimed at providing a way of retrieving variant information using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) and [myvariant.info](http://myvariant.info/). In particular, it is suited for bioinformaticians interested in aggregating variant information into a single NoSQL database (MongoDB solely at the moment). 

###Requirements

- [MongoDB installation](https://docs.mongodb.com/getting-started/python/)
- Python 3 and pip3 for installation 
- [Annovar](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)

###Background

VAPr was developed to simplify the steps required to get mutation data from a VCF file to a downstream analysis process. A query system was implemented allowing users to quickly slice the GV data and select variants according to their characteristics, allowing researchers to focus their analysis only on the subset of data that contains meaningful information. Further, this query system allows the user to select the format in which the data can be retrieved. Most notably, CSV or VCF files can be retrieved from the database, allowing any researcher to quickly filter variants and retrieve them in commonly used formats. 

####Notes

ANNOVAR, alongside with some of their data sets, needs to be installed. The required data sets are the following:

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

They can be downloaded after ANNOVAR has been installed. Head to the directory where ANNOVAR has been installed and run these commands:



```
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb
...
```

**NOTE**: automated database download script is on its way. 

More on how to download and troubleshoot: 
ANNOVAR standard script (this is run automatically, as long as the required databases, as well as ANNOVAR, are installed.

After these steps have been taken, you can proceed and install VAPr as follows:

```pip install VAPr```

Alternatively, clone repository on your desktop. Unzip archive and run, on your terminal:


```
cd variantannotation
python setup.py build
python setup.py install
```

For method 1 and 2, a processed csv file from annovar is required, and it will provide the user with integrated data from annovar and myvariant.info. For methods 3 and 4, the VCF file is enough, and the functions used will create a list of dictionaries with information retrieved from myvariant.info query service.

Here is some sample code for all the methods supplied by the package. It is possible to retrieve the HGVS ID from all the variants contained in a VCF file, thanks to a functionality offered by myvariant.info. It is possible then to integrate the data supplied by myvariant.info databases with ANNOVAR's data.

###Quick-Startup Guide
See this [iPython notebook](https://github.com/ucsd-ccbb/VAPr/blob/master/VAPr_sample_usage.ipynb) for sample usage. 

###Workflow Overview
![Workflow](https://github.com/ucsd-ccbb/VAPr/blob/master/simpler.jpg)

###Available Methods
The package offers 4 different methods to obtain variant data. Two of them require annovar, while the other two are based solely on the use of myvariant.info service. The latter can be used without having to worry about downloading and installing annovar databases, but it tends to return partial or no information for some variants. 

The different methods also enable the user to decide how the data will be parsed to MongoDB. 1 and 3 parse the data by chunks: the user specifies a number of variants (usually 1000), and the data from the vcf and csv files are parsed as soon as those 1000 variants are processed and integrated. This enables huge files to be processed without having to hold them in memory and potentially cause a Memory Overflow error. 

Methods 2 and 4, on the other hand, process the files on their entirety and send them to MongoDB at once.

####Method 1 & 3

####Method 2 & 4


###Data Models
A sample entry in the Mongo Database will look like [this](https://github.com/ucsd-ccbb/VAPr/blob/master/sample_variant_document). The variaty of data that can be retrieved from the sources results from the richness of databases that can be accessed through myvariant.info. However, not every variant will have such data readily available. In some cases, the data will be restricted to what can be inferred from the vcf file and the annotation carried out with annovar. In that case, the entries that will be found in the document will be the following:

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


###Filtering 

Three different pre-made filters that allow for the retrieval of specific variants have been implemented. The filters are in the form of MongoDB queries, and are designed to provide the user with a set of relevant variants. In case the user would like to define its own querying, a template is provided. 
The output of the queries is a list of dictionaries (JSON documents), where each dictionary contains data relative to one variant. 

Further, the package allows the user to parse these variants into an annotated csv or vcf file. 
If needed, annotated, unfiltered vcf and csv files can also be created. They will have the same length (number of variants) as the original files, but will contain much more complete annotation data coming from myvariant.info and ANNOVAR databases. 

To create a csv file, just the filtered output is needed. To create an annotated vcf file, a tab indexed file (.tbi) file is needed (see comments in  section Create unfiltered annotated vcf and csv files at the end of this page). This can be created using tabix.  

First, the file needs to be compressed:

From the command line, running:

`bgzip -c input_file.vcf > input_file.vcf.gz`

returns `input_vcf_file.vcf.gz`

and running 

`tabix input_vcf_file.vcf.gz`

will return: `input_vcf_file.vcf.gz.tbi`

#### Filter #1: specifying cancer-specifc rare variants

 - filter 1: ThousandGenomeAll < 0.05 or info not available
 - filter 2: ESP6500siv2_all < 0.05 or info not available
 - filter 3: cosmic70 information is present
 - filter 4: Func_knownGene is exonic, splicing, or both
 - filter 5: ExonicFunc_knownGene is not "synonymous SNV"
 - filter 6: Read Depth (DP) > 10


        filepath = '.../data/files'
        #Create output files (if needed): specify name of files and path 
        rare_cancer_variants_csv = filepath + "/tumor_rna_rare_cancer_vars_csv.csv"
        rare_cancer_variants_vcf = filepath + "/tumor_rna_rare_cancer_vars_vcf.vcf"
        input_vcf_compressed = filepath + '/test_vcf/Tumor_RNAseq_variants.vcf.gz'
        
        #Apply filter.
        filter_collection = MongoDB_querying.Filters(db_name, collection_name)
        rare_cancer_variants = filter_collection.rare_cancer_variant()
        
        #Crete writer object for filtered lists:
        my_writer = create_output_files.FileWriter(db_name, collection_name)
        
        #cancer variants filtered files
        my_writer.generate_annotated_csv(rare_cancer_variants, rare_cancer_variants_csv)
        my_writer.generate_annotated_vcf(rare_cancer_variants,input_vcf_compressed, rare_cancer_variants_vcf)


### Filter #2: specifying rare disease-specifc (rare) variants

- filter 1: ThousandGenomeAll < 0.05 or info not available
- filter 2: ESP6500siv2_all < 0.05 or info not available
- filter 3: cosmic70 information is present
- filter 4: Func_knownGene is exonic, splicing, or both
- filter 5: ExonicFunc_knownGene is not "synonymous SNV"
- filter 6: Read Depth (DP) > 10
- filter 7: Clinvar data is present 

        filter_collection = MongoDB_querying.Filters(db_name, collection_name)
        rare_disease_variants = filter_collection.rare_disease_variant()

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


### Create unfiltered annotated vcf and csv files 
This output file will contains all annotation data. This may be useful for researchers interested in obtaining a full description of their files.

        #Create output files (if needed): specify name of files and path 
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


**Citations**

Wang K, Li M, Hakonarson H. [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research](http://nar.oxfordjournals.org/content/38/16/e164) , 38:e164, 2010
Xin J, Mark A, Afrasiabi C, Tsueng G, Juchler M, Gopal N, Stupp GS, Putman TE, Ainscough BJ, Griffith OL, Torkamani A, Whetzel PL, Mungall CJ, Mooney SD, Su AI, Wu C (2016) [High-performance web services for querying gene and variant annotation.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9) Genome Biology 17(1):1-7