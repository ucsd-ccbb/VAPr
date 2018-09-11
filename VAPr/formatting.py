
""" 
formatting.py: formats the output from VAPr such that the output matches MAF format, allowing for downstream processing and analysis in Maftools 
"""

__author__ = 'John David Lin', 'Kriti Agrawal'
__date__ = 'Sept. 11, 2018'

### 1. EXTRACT SAMPLES

# Take each sample out from the sample column and append it in its individual row
def extract_samples(list_in):
    list_out = []  # output list
    
    # for every row 
    for x in range(len(list_in)):
        temp_sample_list = list_in[x]['samples'].copy() # temp copy of samples list
        
        # for every sample in samples list
        for y in range(len(temp_sample_list)):
            temp_row_dict = list_in[x].copy()  # copy over the entire original row
            temp_row_dict['samples'] = temp_sample_list[y].copy() # over-write with an individual sample
            list_out.append(temp_row_dict) # add the new dict to output list
    
    return list_out 


### 2. UNNEST DICTIONARIES

# unnest dictionaries recursively
def unnest_dict_core(dict_in, prev, first, dict_out):
    keys= list(dict_in.keys())   # a list of keys
    
    for key in keys:
        
        # set key names
        if first == True:
            key_name = str(key)
        else:
            key_name = prev + "." + str(key)
            
        # check type
        if type(dict_in[str(key)]) == dict:
            unnest_dict_core(dict_in[str(key)], key_name, False, dict_out)
        else:
            dict_out[key_name] = dict_in[str(key)]
    return dict_out

# Call unnest_dict_core and re-assemble dictionaries
def unnest_dict(list_in):
    list_out = []
    for x in range(len(list_in)):
        myDict = {}
        myDict = unnest_dict_core(list_in[x], "", True, myDict)
        list_out.append(myDict)
    return list_out


### 3. UNNEST LIST

# rename dictionaries in the list such that they can be unnested in the same level
def rename_list_content(sub_list_in):
    sub_list_in # make a reference
    out = {}    # output
            
    # for every dictionary in the list
    for x in range(len(sub_list_in)):
        old_keys = list(sub_list_in[x].keys()) # original keys
        new_keys = old_keys.copy() # keys to be modified

        # for every key
        for y in range(len(new_keys)):
            new_keys[y] = new_keys[y] + str(x + 1) # modify the key with index
            out[str(new_keys[y])] = sub_list_in[x][old_keys[y]] # add the modified key and its value to a new dict

    return out

# look for key whose value is a list, call rename_list_content
def find_list(dict_in):
    dict_out = dict_in.copy() # make a copy of the row
    top_keys = list(dict_in.keys())   # a list of keys   
    
    # for every key in that row
    for top_key in top_keys:
        
        # find fields that are lists
        if (type(dict_in[top_key]) == list):
            #print(top_key + " is a list")
            
            has_dict = True # contains dictionaries
            
            # check if the list contains dictionaries
            for item in dict_in[top_key]:
                if type(item) != dict:
                    has_dict = False
                    
            # proceed when the list contains dictionaries
            if has_dict == True:        
                #print(top_key + " contains dictionaries")
                dict_out[top_key] = rename_list_content(dict_in[top_key]) # rename list content
            
    return dict_out

# Unnest lists by calling find_list
def unnest_list(list_in):
    list_out = list_in.copy() # make a copy for output
    
    # for every row
    for x in range(len(list_in)):
        #print(x)
        list_out[x] = find_list(list_in[x]) # find any key whose value is a list

    return list_out


### 4. CHANGE 'samples.AD' TO 't_ref_count' AND 't_alt_count'

# chagne 'samples.AD' to 't_ref_count' and 't_alt_count'
def change_name(dict_in):
    dict_out = {}  # output dict
    keys = list(dict_in.keys())
    for key in keys:
        
        # change 'samples.AD'
        if key == 'samples.AD':
            dict_out['t_ref_count'] = dict_in['samples.AD'][0]
            dict_out['t_alt_count'] = dict_in['samples.AD'][1]            
        else:
            dict_out[key] = dict_in[key]  # copy over other fields
    return dict_out

# Call change_name to change 'samples.AD' to 't_ref_count' and 't_alt_count'
def samples_AD(list_in):
    list_out = []  # output list
    
    # for every row
    for x in range(len(list_in)):
        list_out.append(change_name(list_in[x]))
    return list_out


### 5. RECORD Variant_Type

# Record variant type according to ref and alt alleles
def varType(list_in):
    for x in range(len(list_in)):
        ref = len(list_in[x]['ref'])
        alt = len(list_in[x]['alt'])
        if(ref == alt == 1):
            list_in[x]['Variant_Type'] = 'SNP'
        if(ref == alt == 2):
            list_in[x]['Variant_Type'] = 'DNP'
        if(ref == alt == 3):
            list_in[x]['Variant_Type'] = 'TNP'
        if(ref == alt and ref > 3):
            list_in[x]['Variant_Type'] = 'ONP'
        if(ref < alt):
            list_in[x]['Variant_Type'] = 'INS'
        if(ref > alt):
            list_in[x]['Variant_Type'] = 'DEL'
            
    return list_in 


### 6. RE-ARRANGE AND RENAME COLUMNS

# Re-arrange and rename columns to match the MAF format
def change_cols(df):
    cols = list(df)
    cols.insert(0, cols.pop(cols.index('gene_knowngene')))
    cols.insert(1, cols.pop(cols.index('chr')))
    cols.insert(2, cols.pop(cols.index('start')))
    cols.insert(3, cols.pop(cols.index('end')))
    cols.insert(4, cols.pop(cols.index('ref')))
    cols.insert(5, cols.pop(cols.index('alt')))
    cols.insert(6, cols.pop(cols.index('Variant_Type')))
    cols.insert(7, cols.pop(cols.index('func_knowngene')))
    cols.insert(8, cols.pop(cols.index('samples.sample_id')))
    cols.insert(9, cols.pop(cols.index('dbsnp.rsid')))
    cols.insert(10, cols.pop(cols.index('t_ref_count')))
    cols.insert(11, cols.pop(cols.index('t_alt_count')))

    df=df.loc[:,cols]

    df = df.rename(columns={'gene_knowngene':'Hugo_Symbol', 
                            'chr':'Chromosome', 
                            'start':'Start_Position', 
                            'end':'End_Position',
                            'ref':'Reference_Allele', 
                            'alt':'Tumor_Seq_Allele2', 
                            'func_knowngene':'Variant_Classification', 
                            'samples.sample_id':'Tumor_Sample_Barcode',
                            'cadd.1000g.afr': 'AFR_MAF',
                            'cadd.1000g.amr': 'AMR_MAF',
                            'cadd.1000g.asn': 'ASN_MAF',
                            'cadd.1000g.eur': 'EUR_MAF',
                            'dbsnp.rsid': 'dbSNP_RS'
                           })
    
    df['Variant_Classification']=df['Variant_Classification'].replace({
        'intronic':'Intron', 
        'intergenic':'IGR',
        'UTR3':"3'UTR", 
        "UTR5":"5'UTR",
        'downstream':"3'Flank", 
        'upstream':"5'Flank",
        'splicing':'Splice_Site', 'ncRNA_exonic':'RNA',
        'ncRNA_intronic':'RNA', 'ncRNA_UTR3':'RNA',
        'ncRNA_UTR5':'RNA', 'ncRNA':'RNA',})
    
    mask = df.exonicfunc_knowngene == 'nonsynonymous SNV'
    df.loc[mask, 'Variant_Classification'] = "Missense_Mutation"
    mask = df.exonicfunc_knowngene == 'synonymous SNV'
    df.loc[mask, 'Variant_Classification'] = "Silent"

    return df


