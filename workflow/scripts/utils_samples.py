#########################################
# scripts/utils_samples.py
#########################################

import os
import glob
import sys
import pandas as pd


def load_metadata(file):
    df = pd.read_excel(file)
    mapping_metadata = dict(zip( df ["NGS_ID"], df["Sample_ID"]))
    #print(mapping_metadata)  # {'SRR6669849': 'B', 'SRR6669961': 'F'}
    return mapping_metadata


def make_samples_list(samples_file):
    read_1_suffixes = ["_R1_001", "_1", "_R1"]
    read_2_suffixes = ["_R2_002", "_2", "_R2"]
    #read_1_suffixes = ["_1", "_R1"]
    #read_2_suffixes = ["_2", "_R2"]
    extensions = [".fastq.gz", ".fq.gz"]
    
    samples_list_NGS = []

    list_files = os.listdir(samples_file)

    for i in range(0, len(list_files), 2):
        r1 = list_files[i]     
        r2 = list_files[i+1]   
        def detect_suffix_ext(filepath, read_suffixes): 
            base = os.path.basename(filepath).replace(" ", "")
            #print(base) # SRR6669849_1.fastq.gz
            for sfx in read_suffixes:
                for ext in extensions:
                    if base.endswith(sfx + ext):
                        return sfx, ext
            re_suf=", ".join(read_suffixes)
            re_ext=", ".join(extensions)
            raise ValueError(f"Filename does not follow the expected pattern: {filepath}. Expected suffixes: {re_suf} and extensions: {re_ext}. Or check if the raw data is completed.")

        r1_sfx, r1_ext = detect_suffix_ext(r1, read_1_suffixes)
        r2_sfx, r2_ext = detect_suffix_ext(r2, read_2_suffixes)

        sample_name = os.path.basename(r1)[: -len(r1_sfx + r1_ext)]

        samples_list_NGS.append([sample_name, r1_sfx, r2_sfx, r1_ext])
    
    return samples_list_NGS

def mapping_samples_list(samples_list_NGS, mapping_metadata):
    
    final_mapping = {}
  
    for sample in samples_list_NGS:
        ngs_id = sample[0]
        read_1 = sample[1]
        read_2 = sample[2]
        extensions = sample[3]
        suffixes = sample[1:]
        
        if ngs_id in mapping_metadata:
            sample_id = mapping_metadata[ngs_id]
            final_list = [sample_id] + suffixes
            final_mapping[ngs_id] = final_list
        else:
            continue
            #sys.exit(f"Error: NGS_ID {ngs_id} was not found in the metadata table.")

    #short_inverse = { v[0]: k for k, v in final_mapping.items() }   # {'B': 'SRR6669849', 'F': 'SRR6669961'}
    short_inverse = { v[0]: [k] + v[1:] for k, v in final_mapping.items() }
    
    return final_mapping, short_inverse  #{'SRR6669849': ['B', '_1', '_2', '.fastq.gz'], 'SRR6669961': ['F', '_1', '_2', '.fastq.gz']} --> final_mapping

