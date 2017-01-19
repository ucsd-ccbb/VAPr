**variantannotation**

This package is aimed at providing a way of retrieving variant information using ANNOVAR http://annovar.openbioinformatics.org/en/latest/ and [myvariant.info](http://myvariant.info/). In particular, it is suited for bioinformaticians interested in aggregating variant information into a single database. The data structure used is a python dictionaty, so the data can be easily parsed to a MongoDBinstance.

Required library: myvariant, pymongo. The other ones should be installed natively. Run:

```pip install myvariant```

To get the package. For more info, https://github.com/SuLab/myvariant.py

Basic steps to set up pymongo: install mongodb for your operating system, create a directory /data/db wherever you'd like to store your data.
Install pymongo:

```pip install pymongo```

Then run, on your terminal:

```mongod --dbpath /data/db```

And as soon as you run the scripts from variantannotaiton the data will automatically be stored to it. Database and collection name should be specified (refer to sample code below).
For pymongo, and more information on how to set up a Mongo Database: https://docs.mongodb.com/getting-started/python/

If you don't have pymongo, you can still get the results in a list, where each item in the list will be a dictionary (JSON formatted) containing the information for the variants coming from both ANNOVAR and myvariant.info.

Other required software tools: ANNOVAR, with some of their datasets installed. The required datasets are the following:

- knownGene
- tfbsConsSites
- cytoBand
- genomicSuperDups
- gwasCatalog
- esp6500siv2_all
- 1000g2015aug_all
- snp138
- ljb26_all
- cg46
- cg69
- popfreq_all
- clinvar_20140929
- cosmic70
- nci60

They can be downloaded after ANNOVAR has been installed. Head to the directory where ANNOVAR has been installed and run these commands


```
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar knownGene humandb
perl annotate_variation.pl -buildver hg19 -downdb -webfrom annovar esp6500siv2_all humandb
...
```

More on how to download and troubleshoot: http://annovar.openbioinformatics.org/en/latest/user-guide/download/

ANNOVAR standard script (this is run automatically, as long as the required databases, as well as ANNOVAR, are installed.

```sudo perl /database/annovar/table_annovar.pl (/filepath) /database/annovar/humandb/ -buildver hg19 -out (/output filepath and name) -remove -protocol knownGene,tfbsConsSites,cytoBand,targetScanS,genomicSuperDups,gwasCatalog,esp6500siv2_all,1000g2015aug_all,snp138,ljb26_all,cg46,cg69,popfreq_all,clinvar_20140929,cosmic70,nci60 -operation g,r,r,r,r,r,f,f,f,f,f,f,f,f,f,f -nastring . -vcfinput -csvout```

After these steps have been taken, you can proceed and install variantannotation as follows:

```pip install variantannotation```

Alternatively, clone repository on your desktop. Unzip archive and run, on your terminal:


```
cd variantannotation
python setup.py build
python setup.py install
```
For method 1 and 2, a processed csv file from annovar is required, and it will provide the user with integrated data from annovar and myvariant.info. For methods 3 and 4, the VCF file is enough, and the functions used will create a list of dictionaries with information retrieved from myvariant.info query service.

Here is some sample code for all the methods supplied by the package. It is possible to retrieve the HGVS ID from all the variants contained in a VCF file, thanks to a functionality offered by myvariant.info. It is possible then to integrate the data supplied by myvariant.info databases with ANNOVAR's data.

See this [iPython notebook](https://github.com/Mazzafish/VAPr/blob/master/variantannotation_sample_usage.ipynb) for sample usage.


**Citations**

Wang K, Li M, Hakonarson H. [ANNOVAR: Functional annotation of genetic variants from next-generation sequencing data Nucleic Acids Research](http://nar.oxfordjournals.org/content/38/16/e164) , 38:e164, 2010
Xin J, Mark A, Afrasiabi C, Tsueng G, Juchler M, Gopal N, Stupp GS, Putman TE, Ainscough BJ, Griffith OL, Torkamani A, Whetzel PL, Mungall CJ, Mooney SD, Su AI, Wu C (2016) [High-performance web services for querying gene and variant annotation.](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0953-9) Genome Biology 17(1):1-7