import pandas as pd
""" 
formatting.py: formats the output from VAPr such that the output matches MAF format, allowing for downstream processing 
and analysis in Maftools 
"""

__author__ = 'John David Lin', 'Kriti Agrawal'
__date__ = 'Sept. 11, 2018'

### 1. EXTRACT SAMPLES


# Take each sample out from the sample column and append it in its individual row
def extract_samples(dataset_list_in):
    dataset_list_out = []  # output list

    # for every row
    for x in range(len(dataset_list_in)):
        temp_sample_list = dataset_list_in[x]['samples'].copy()  # temp copy of samples list

        # for every sample in samples list
        for y in range(len(temp_sample_list)):
            temp_row_dict = dataset_list_in[x].copy()  # copy over the entire original row
            temp_row_dict['samples'] = temp_sample_list[y].copy()  # over-write with an individual sample
            dataset_list_out.append(temp_row_dict)  # add the new dict to output list

    return dataset_list_out


### 2. UNNEST DICTIONARIES

# Call unnest_dict_core and re-assemble dictionaries
def unnest_dict(dataset_list_in):
    dataset_list_out = []
    for x in range(len(dataset_list_in)):
        temp_dict = {}
        temp_dict = unnest_dict_core(dataset_list_in[x], "", True, temp_dict)
        dataset_list_out.append(temp_dict)
    return dataset_list_out


# unnest dictionaries recursively
def unnest_dict_core(nested_dict_in, prev_key, first_level, unnested_dict_out):
    keys = list(nested_dict_in.keys())  # a list of keys

    for key in keys:

        # set key names
        if first_level:
            key_name = str(key)
        else:
            key_name = prev_key + "." + str(key)

        # check type
        if type(nested_dict_in[str(key)]) == dict:
            unnest_dict_core(nested_dict_in[str(key)], key_name, False, unnested_dict_out)
        else:
            unnested_dict_out[key_name] = nested_dict_in[str(key)]
    return unnested_dict_out


### 3. UNNEST LIST

# A wrapper that unnnest lists by calling find_list
def unnest_list(dataset_list_in):
    dataset_list_out = dataset_list_in.copy()  # make a copy for output

    # for every row
    for x in range(len(dataset_list_in)):
        # print(x)
        dataset_list_out[x] = find_list(dataset_list_in[x])  # find any key whose value is a list

    return dataset_list_out


# Finds field that are dictionaries. Find column in the table that contain a list whose keys that would be identical
# when unnested and calls rename list. Look for key whose value is a list, call rename_list_content.
def find_list(dict_in):
    dict_out = dict_in.copy()  # make a copy of the row
    dict_keys = list(dict_in.keys())  # a list of keys

    # for every key in that row
    for dict_key in dict_keys:

        # find fields that are lists
        if type(dict_in[dict_key]) == list:
            # print(dict_key + " is a list")

            has_dict = True  # contains dictionaries

            # check if the list contains dictionaries
            for key_values in dict_in[dict_key]:
                if type(key_values) != dict:
                    has_dict = False

            # proceed when the list contains dictionaries
            if has_dict:
                # print(dict_key + " contains dictionaries")
                dict_out[dict_key] = rename_list_content(dict_in[dict_key])  # rename list content

    return dict_out


# Since the keys can have the same names in the dictionary, this function renames the items by appending a number at
# the end of each item. Having keys with the same name in a dictionary throws an error.
# rename dictionaries in the list such that they can be unnested in the same level and and keys can be unique.
def rename_list_content(same_keys_list):
    new_keys_dict_out = {}  # output

    # for every dictionary in the list
    for x in range(len(same_keys_list)):
        old_keys = list(same_keys_list[x].keys())  # original keys
        new_keys = old_keys.copy()  # keys to be modified

        # for every key
        for y in range(len(new_keys)):
            new_keys[y] += str(x + 1)  # modify the key with index

            # add the modified key and its value to a new dict
            new_keys_dict_out[str(new_keys[y])] = same_keys_list[x][old_keys[y]]

    return new_keys_dict_out


### 4. CHANGE 'samples.AD' TO 't_ref_count' AND 't_alt_count'

# Call change_name to change 'samples.AD' to 't_ref_count' and 't_alt_count', because they are required fields in MAf.
def samples_AD(dataset_list_in):
    dataset_list_out = []  # output list

    # for every row
    for x in range(len(dataset_list_in)):
        dataset_list_out.append(change_name(dataset_list_in[x]))
    return dataset_list_out


# Change 'samples.AD' to 't_ref_count' and 't_alt_count'
def change_name(sample_dict_in):
    sample_dict_out = {}  # output dict
    keys = list(sample_dict_in.keys())
    for key in keys:

        # change 'samples.AD'
        if key == 'samples.AD':
            sample_dict_out['t_ref_count'] = sample_dict_in['samples.AD'][0]
            sample_dict_out['t_alt_count'] = sample_dict_in['samples.AD'][1]
        else:
            sample_dict_out[key] = sample_dict_in[key]  # copy over other fields
    return sample_dict_out


### 5. RECORD Variant_Type

# Record variant type according to ref and alt alleles. Determines the variant type based on the
# number of ref and alt alleles
def varType(dataset_list_in):
    for x in range(len(dataset_list_in)):
        ref = len(dataset_list_in[x]['ref'])
        alt = len(dataset_list_in[x]['alt'])
        if ref == alt == 1:
            dataset_list_in[x]['Variant_Type'] = 'SNP'
        if ref == alt == 2:
            dataset_list_in[x]['Variant_Type'] = 'DNP'
        if ref == alt == 3:
            dataset_list_in[x]['Variant_Type'] = 'TNP'
        if ref == alt and ref > 3:
            dataset_list_in[x]['Variant_Type'] = 'ONP'
        if ref < alt:
            dataset_list_in[x]['Variant_Type'] = 'INS'
        if ref > alt:
            dataset_list_in[x]['Variant_Type'] = 'DEL'

    return dataset_list_in


### 6. Unnest multiple values separated by ";" in a single field

# Fields with multiple values to separate: gene_knowngene, genedetail_knowngene, func_knowngene
# Takes in the whole table. Find fields with multiple values and attemp to separates those values and put each in a
# new row.
def unnest_semicolon_values(data_list_in):
    debug_genes = False
    debug_genes_and_func = False
    debug_genes_and_genedetails = False
    debug_genes_and_func_and_genedetails = False
    debug_lengths_check = False
    multiple_values_count = 0  # testing use

    data_list_out = []  # the output table to be returned

    # for every row in the whole table
    for row in data_list_in:
        if row['gene_knowngene'].find(';') is not -1:

            # testing use
            if debug_genes:
                print(row['gene_knowngene'])

            temp_rows_list = []  # holds row copies temporarily
            gene_knowngene_values = row['gene_knowngene'].split(';')

            # for every gene after splitting
            for gene in gene_knowngene_values:
                temp_row_dict = row.copy()  # copy the whole row
                temp_row_dict['gene_knowngene'] = gene  # replace with individual gene
                temp_rows_list.append(temp_row_dict)  # put rows in a list temporarily

            func_knowngene_values = row['func_knowngene'].split(';')

            # Match genes with corresponding funcs only when there are the same numbers of them
            # In the case of multiple genes with only 1 func, that 1 func is copied over to a new row
            #    with every gene, so that case is always taken care of.
            if len(func_knowngene_values) == len(gene_knowngene_values):
                # for every new row created after splitting genes
                for index in range(len(temp_rows_list)):
                    temp_rows_list[index]['func_knowngene'] = func_knowngene_values[
                        index]  # replace with individual func

                # testing use
                if debug_genes_and_func:
                    print(row['gene_knowngene'])
                    print(row['func_knowngene'] + "\n")
                    for x in temp_rows_list:
                        print(x['gene_knowngene'] + "  " + x['func_knowngene'])
                    print('\n')

            if 'genedetail_knowngene' in row.keys():
                genedetail_knowngene_values = row['genedetail_knowngene'].split(';')

                # Match genes with corresponding genedetails only when there are the same numbers of them
                # In the case of multiple genes with only 1 genedetail, that 1 genedetail is copied over to a new row
                #    with every gene, so that case is always taken care of.
                if len(genedetail_knowngene_values) == len(gene_knowngene_values):

                    # for every new row created after splitting genes
                    for index in range(len(temp_rows_list)):
                        temp_rows_list[index]['genedetail_knowngene'] = genedetail_knowngene_values[index]

                    # testing use
                    if debug_genes_and_genedetails:
                        print(row['gene_knowngene'])
                        print(row['func_knowngene'])
                        print(row['genedetail_knowngene'] + "\n")
                        for x in temp_rows_list:
                            print(x['gene_knowngene'] + "  " + x['func_knowngene'] + "  " + x['genedetail_knowngene'])
                        print('\n')

                    # testing use
                    if debug_genes_and_func_and_genedetails:
                        if len(genedetail_knowngene_values) == len(gene_knowngene_values) == len(func_knowngene_values):
                            print(row['gene_knowngene'])
                            print(row['func_knowngene'])
                            print(row['genedetail_knowngene'] + "\n")
                            for x in temp_rows_list:
                                print(
                                x['gene_knowngene'] + "  " + x['func_knowngene'] + "  " + x['genedetail_knowngene'])
                            print('\n')

            # testing use
            if debug_genes:
                for x in temp_rows_list:
                    print(x['gene_knowngene'] + "  " + x['func_knowngene'])

            # write the newly created rows to the output table
            for temp_row in temp_rows_list:
                data_list_out.append(temp_row)

                # testing use
                if debug_lengths_check:
                    multiple_values_count += 1

            # testing use
            if debug_lengths_check:
                multiple_values_count -= 1

        # write the row to the output table as it is
        else:
            data_list_out.append(row.copy())

    # testing use
    if debug_lengths_check:
        print("The lenght of data input: " + str(len(data_list_in)))
        print("The lenght of data output: " + str(len(data_list_out)))
        print("The difference in lengths: " + str(multiple_values_count))
        print("New rows created: " + str(len(data_list_in) - len(data_list_out)))

    return data_list_out


### 7. RE-ARRANGE AND RENAME COLUMNS

# Re-arrange and rename columns to match the MAF format
def change_cols(df):
    required_cols = ['gene_knowngene', 'chr', 'start', 'end', 'ref', 'alt', 'Variant_Type', 'func_knowngene',
                     'samples.sample_id', 'dbsnp.rsid', 't_ref_count', 't_alt_count', 'aachange_knowngene']
    for col in required_cols:
        if col not in df.columns.values:
            df[col] = ""

    df_req = df[required_cols]
    df_extra = df.drop(required_cols, axis=1)
    df = pd.concat([df_req, df_extra], axis=1)

    df = df.rename(columns={'gene_knowngene': 'Hugo_Symbol',
                            'chr': 'Chromosome',
                            'start': 'Start_Position',
                            'end': 'End_Position',
                            'ref': 'Reference_Allele',
                            'alt': 'Tumor_Seq_Allele2',
                            'func_knowngene': 'Variant_Classification',
                            'samples.sample_id': 'Tumor_Sample_Barcode',
                            'cadd.1000g.afr': 'AFR_MAF',
                            'cadd.1000g.amr': 'AMR_MAF',
                            'cadd.1000g.asn': 'ASN_MAF',
                            'cadd.1000g.eur': 'EUR_MAF',
                            'dbsnp.rsid': 'dbSNP_RS',
                            'aachange_knowngene': 'Protein_Change'
                            })

    df['Variant_Classification'] = df['Variant_Classification'].replace({
        'intronic': 'Intron',
        'intergenic': 'IGR',
        'UTR3': "3'UTR",
        "UTR5": "5'UTR",
        'downstream': "3'Flank",
        'upstream': "5'Flank",
        'splicing': 'Splice_Site', 'ncRNA_exonic': 'RNA',
        'ncRNA_intronic': 'RNA', 'ncRNA_UTR3': 'RNA',
        'ncRNA_UTR5': 'RNA', 'ncRNA': 'RNA', })

    mask = df.exonicfunc_knowngene == 'nonsynonymous SNV'
    df.loc[mask, 'Variant_Classification'] = "Missense_Mutation"
    mask = df.exonicfunc_knowngene == 'synonymous SNV'
    df.loc[mask, 'Variant_Classification'] = "Silent"
    mask = df.exonicfunc_knowngene == 'stopgain'
    df.loc[mask, 'Variant_Classification'] = "Nonsense_Mutation"
    mask = df.exonicfunc_knowngene == 'stoploss'
    df.loc[mask, 'Variant_Classification'] = "Nonstop_Mutation"
    mask = df.exonicfunc_knowngene == 'frameshift insertion'
    df.loc[mask, 'Variant_Classification'] = "Frame_Shift_Ins"
    mask = df.exonicfunc_knowngene == 'frameshift deletion'
    df.loc[mask, 'Variant_Classification'] = "Frame_Shift_Del"
    mask = df.exonicfunc_knowngene == 'nonframeshift insertion'
    df.loc[mask, 'Variant_Classification'] = "In_Frame_Ins"
    mask = df.exonicfunc_knowngene == 'nonframeshift deletion'
    df.loc[mask, 'Variant_Classification'] = "In_Frame_Del"

    return df


### 8. WRAPPER

# Formats the output from VAPr such that the output matches MAF format, allowing for downstream processing and
# analysis in Maftools A wrapper for the functions above
# Input: a VAPr output list
# Output: a formatted dataframe ready to be saved as a MAF file
def maf_formatter(dataset_list_in):
    import pandas as pd
    k0 = dataset_list_in.copy()
    k1 = extract_samples(k0)
    k2 = unnest_dict(k1)
    k3 = unnest_list(k2)
    k4 = unnest_dict(k3)
    k5 = samples_AD(k4)
    k6 = varType(k5)
    k7 = unnest_semicolon_values(k6)

    df = pd.DataFrame(data=k7)
    maf = change_cols(df)

    return maf

#########################################################################################
### CREATE A LIST OF THE WHOLE DATASET
def create_whole_dataset_list(MONGODB, COLLECTION):
     #access the mongodb database
    from pymongo import MongoClient
    c = MongoClient()
    c.test_database
    db = c[MONGODB]
    col = db[COLLECTION]
    
    #create a list of all the values in the database
    cursor = col.find({})
    wholeList = []
    for document in cursor:
        wholeList.append(document)

    return wholeList
