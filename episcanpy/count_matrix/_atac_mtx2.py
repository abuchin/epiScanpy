import pybedtools


# the codes are from EpiScanpy but adapt for these input files by using pybedtools

HUMAN = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', '2',
         '20', '21', '22', '3', '4', '5', '6', '7', '8', '9', 'X', 'Y']

MOUSE = ['1', '10', '11', '12', '13', '14', '15', '16', '17', '18', '19', 
        '2', '3', '4', '5', '6', '7', '8', '9','X', 'Y']


def load_features(file_features, chromosomes=HUMAN, path="", sort=False):
   
    features_chrom = {}
    for c in chromosomes:
        features_chrom[c] = []
    with open(path + file_features) as features:
        for line in features:
            ar = line.strip().split("\t")
            if ar[0][3:] in chromosomes:
                features_chrom[ar[0][3:]].append([int(ar[1]), int(ar[2])])
                         
    if sort == True:
        for c in chromosomes:
            sorted(features_chrom[c], key=lambda x: x[0])
            
    return(features_chrom)

def name_features(loaded_features):
    feat_names = []
    i = 0
    for c in loaded_features.keys():
        for name in loaded_features[c]:
            add_name = '_'.join([c, str(name[0]), str(name[1])])
            add_name ='chr' + add_name
            if add_name[-1] =='\n':
                add_name = add_name[:-1]
            feat_names.append(add_name)
            i += 1
    return(feat_names)

def bld_atac_mtx_pybedtools(list_bam_files, loaded_feat, output_file_name=None,
    path=None, writing_option='a', header=None, mode='rb',
    check_sq=True, chromosomes=HUMAN):
       
    if output_file_name==None:
        output_file_name='std_output_ct_mtx.txt'


    if path==None:
        path=''
    
    # open file to write
    output_file = open(path+output_file_name, writing_option)
    # write header if specified
    if header != None:
        output_file.write('sample_name\t')
        for feature in header:
            output_file.write(feature)
            output_file.write('\t')
        output_file.write('\n')
    # close file to write   
    output_file.close()

    # start going through the bam files
    for name_file in list_bam_files[0:]:

        ## important variables for output
        index_feat = {key: 0 for key in chromosomes}    
        val_feat = {key: [0 for x in range(len(loaded_feat[key]))] for key in chromosomes}
    
        ## PART 1 read the bam file
        keep_lines = []

        samfile = pybedtools.example_bedtool(name_file)

        co = 0
        for read in samfile:
            co = co + 1
            try:
                line = str(read).split('\t')
                #print(str(read))
                if line[2][3:] in chromosomes:
                    keep_lines.append(line[2:4])
            except:
                print("Error: from reading "+name_file+" at line "+str(co)+"\n")
                break
                
            ### print -- output
        print(name_file, len(keep_lines), 'mapped reads')
        #samfile.close()
        
        ## PART2 reads that fall into 
        for element in keep_lines:
            ## 2 things per line:
            chrom = element[0][3:]
            read_pos = int(element[1])
            max_value_index = len(loaded_feat[chrom])
            ## I want to check if the read map to a feature in the same chrom
            pointer_feat_pos = index_feat[chrom]
            for feat_pos in loaded_feat[chrom][pointer_feat_pos:]:
                pointer_feat_pos += 1
                # Go through all features that are smaller than the read position
                if read_pos > feat_pos[1]:
                    continue
                # if read_pos fall in a feature
                elif read_pos > feat_pos[0]:
                    # Update the pointer for the next read if the pointer isn't out of range
                    if pointer_feat_pos < max_value_index:
                        index_feat[chrom] = pointer_feat_pos
                        val_feat[chrom][pointer_feat_pos] += 1
                    else:
                        index_feat[chrom] = max_value_index
                    # Check the following features without updating the pointer. 
                    break
                else:
                    break
     
            for feat_pos in loaded_feat[chrom][pointer_feat_pos:]:
                # +1 if reads fall into more than one feature
                if feat_pos[0] < read_pos:
                    val_feat[chrom][pointer_feat_pos] += 1
                    pointer_feat_pos += 1
                # if read_pos > start position of the new feature break
                elif read_pos < feat_pos[0]:
                    break
                else:
                    print('error')
                    break
                    
        # PART 3
        # open
        output_file = open(path+output_file_name, 'a')
        # write down the result of the cell
        output_file.write(name_file)
        output_file.write('\t')
        for chrom in chromosomes:
            output_file.write('\t'.join([str(p) for p in val_feat[chrom]]))
            output_file.write('\t')
        output_file.write('\n')
        #close
        output_file.close()