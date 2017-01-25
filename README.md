##VAPr 
####Variant Annotation and Prioritization package

This package is aimed at providing a way of retrieving variant information using [ANNOVAR](http://annovar.openbioinformatics.org/en/latest/) and [myvariant.info](http://myvariant.info/). In particular, it is suited for bioinformaticians interested in aggregating variant information into a single NoSQL database (MongoDB solely at the moment). 

###Requirements

- [MongoDB installation](https://docs.mongodb.com/getting-started/python/)
- Python 3 and pip3 for installation 
- [Annovar](http://annovar.openbioinformatics.org/en/latest/user-guide/download/)

###Background

VAPr was developed to simplify the steps required to get mutation data from a VCF file to a downstream analysis process. A query system was implemented allowing users to quickly slice the GV data and select variants according to their characteristics, allowing researchers to focus their analysis only on the subset of data that contains meaningful information. Further, this query system allows the user to select the format in which the data can be retrieved. Most notably, CSV or VCF files can be retrieved from the database, allowing any researcher to quickly filter variants and retrieve them in commonly used formats. 

####Notes

ANNOVAR, with some of their datasets needs to be installed. The required datasets are the following:

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
See this [iPython notebook](https://github.com/ucsd-ccbb/VAPr/blob/master/VAPr_sample_usage.ipynb) for sample usag


###Workflow Overview
![Workflow](https://github.com/ucsd-ccbb/VAPr/blob/master/simpler.jpg)


**Citations**

Wang K, Li M, Hakonarson H. [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research](http://nar.oxfordjournals.org/content/38/16/e164) , 38:e164, 2010
Xin J, Mark A, Afrasiabi C, Tsueng G, Juchler M, Gopal N, Stupp GS, Putman TE, Ainscough BJ, Griffith OL, Torkamani A, Whetzel PL, Mungall CJ, Mooney SD, Su AI, Wu C (2016) [High-performance web services for querying gene and variant annotation.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9) Genome Biology 17(1):1-7