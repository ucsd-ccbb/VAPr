
if __name__ == '__main__':


    mapping = [{'csv_file_basename': 'BC001.final_annotated',
                'csv_file_full_path': '/Volumes/Carlo_HD1/CCBB/VAPr_files/csv_benchmark/sample_BC001',
                'num_samples_in_csv': 1,
                'raw_vcf_file_full_path': '/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_benchmark/sample_BC001/BC001.final.vcf',
                'sample_names': ['BC001'],
                'vcf_file_basename': 'BC001.final.vcf',
                'vcf_sample_dir': '/Volumes/Carlo_HD1/CCBB/VAPr_files/vcf_benchmark/sample_BC001'}]

    def get_sample_csv_vcf_tuple(mapping):
        """ Locate files associated with a specific sample """

        list_tupls = []
        db = 'VariantDBBenchmarkWES'
        collection = 'collect'
        for _map in mapping:

            matching_csv = [i for i in os.listdir(_map['csv_file_full_path']) if i.startswith(_map['csv_file_basename'])
                            and i.endswith('txt')]

            matching_vcf = [i for i in os.listdir(_map['csv_file_full_path']) if i.startswith(_map['csv_file_basename'])
                            and i.endswith('vcf')]

            print(matching_vcf, matching_csv)

            if len(matching_csv) > 1 or len(matching_vcf) > 1:
                raise ValueError('Too many matching csvs')
            elif len(matching_csv) == 0 or len(matching_vcf) == 0:
                raise ValueError('Csv not found')
            else:
                csv_path = os.path.join(_map['csv_file_full_path'], matching_csv[0])
                vcf_path = os.path.join(_map['csv_file_full_path'], matching_vcf[0])
                list_tupls.append((_map['sample_names'],
                                   vcf_path,
                                   csv_path,
                                   db,
                                   collection))

        return list_tupls

    list_tupls = get_sample_csv_vcf_tuple(mapping)

    for tpl in list_tupls:
        hgvs = HgvsParser(tpl[1])
        csv_parsing = TxtParser(tpl[2], samples=hgvs.samples)
        num_lines = csv_parsing.num_lines
        print(num_lines, 5000)
        n_steps = int(num_lines / 5000) + 1
        print(n_steps)

    _tuple = list_tupls[0]
    new_tuple_list = []
    for i in range(n_steps):
        sample = _tuple[0]
        vcf_file = _tuple[1]
        csv_file = _tuple[2]
        db_name = _tuple[3]
        collection_name = _tuple[4]
        step = i
        new_tuple_list.append((sample,
                               vcf_file,
                               csv_file,
                               db_name,
                               collection_name,
                               step))

    # print(new_tuple_list[0])
    parse_by_step(new_tuple_list[0])