""" LIBRARIES """
import pandas as pd
import numpy as np
import re

""" FUNCTIONS """
# Precompile the regular expression
pattern = re.compile("ENSG...........")

def load_data(file_path):
    df = pd.read_csv(file_path, delimiter="\t", dtype= object, compression='gzip')
    return df

def process_data(data):
    pass

def summarize_data(summary_data):
    pass

def write_sumary_to_file(summary_data):
    pass 

# Eliminating repeats
def no_repeats(nparray):
    return [np.unique(arr) for arr in nparray]

# Mining gene ID
def extract(lists):
    return[pattern.findall(lst) for lst in lists]

# Looking for differences
def difference(array1,array2):
    return [np.setdiff1d(arr1, arr2, assume_unique=True) for arr1, arr2 in zip(array1, array2)]

# Intercepts
def inter(array1,array2):
    return [np.intersect1d(arr1, arr2, assume_unique=True) for arr1, arr2 in zip(array1, array2)]

# Counting
def counter(arrays):
    return sum(1 for arr in arrays if arr.size > 0)

def annotation_agreement_rate(uni_clean, AN, SN, VP):
    # Initialize a dictionary to store various counts.
    counters = {
        'all_agree': 0,
        'two_agree': 0,
        'one_agree': 0,
        'none_agree': 0,
        'AN_SN_agree_2': 0,
        'SN_VP_agree_2': 0,
        'AN_VP_agree_2': 0,
        'AN_agree_1': 0,
        'SN_agree_1': 0,
        'VP_agree_1': 0,
    }

    # Mapping for two-agree scenarios.
    two_agree_mapping = {
        'AS': 'AN_SN_agree_2',
        'SV': 'SN_VP_agree_2',
        'AV': 'AN_VP_agree_2'
    }
    
    # Mapping for one-agree scenarios.
    one_agree_mapping = {
        'A': 'AN_agree_1',
        'S': 'SN_agree_1',
        'V': 'VP_agree_1'
    }

    # Iterate through the lists
    for uni, an, sn, vp in zip(uni_clean, AN, SN, VP):
        temp = 0  # Temporary counter for the loop.
        des = []  # List to store which arrays agree with uni_clean.
        
        zize = len(uni)  # Length of the union set for this iteration.
        
        # Check if AN agrees with uni_clean and update counters.
        if zize == len(an):
            temp += 1
            des.append('A')
            
        # Check if SN agrees with uni_clean and update counters.
        if zize == len(sn):
            temp += 1
            des.append('S')
            
        # Check if VP agrees with uni_clean and update counters.
        if zize == len(vp):
            temp += 1
            des.append('V')
        
        # Convert list to a string for easier matching.
        des_str = ''.join(des)
        
        # Update the relevant counter based on how many arrays agreed with uni_clean.
        if temp == 3:  # All agree
            counters['all_agree'] += 1
        elif temp == 2:  # Exactly two agree
            counters['two_agree'] += 1
            counters[two_agree_mapping[des_str]] += 1
        elif temp == 1:  # Only one agrees
            counters['one_agree'] += 1
            counters[one_agree_mapping[des_str]] += 1
        elif temp == 0:  # None agree
            counters['none_agree'] += 1

    return counters  # Return the counters dictionary.




def main():

    """ DATA PROCESSING """
    
    """ File Input & Data Mining """
    # Reading file
    file_path = "https://github.com/quemeb/USC_research/raw/main/Huaiyu/AnnoQ/Test_data.txt.gz"
    cdata = load_data(file_path)
    
    # Selecting data
    AN_ID_genic = cdata["ANNOVAR_ensembl_Gene_ID"]  #Gene ID should be here
    AN_ID_intergenic = cdata["ANNOVAR_ensembl_Closest_gene(intergenic_only)"]    #alternative place for Gene ID
    AN_ID_tog = [AN_ID_intergenic[i] if gene_id == "." else gene_id for i, gene_id in enumerate(AN_ID_genic)]
    
    SN_ID_tog = cdata["SnpEff_ensembl_Gene_ID"]   #Gene ID for SnpEff
    SN_ID_intergenic = [SN_ID_tog[i] if gene_id == '.' else '.' for i, gene_id in enumerate(AN_ID_genic)]
    SN_ID_genic = [SN_ID_tog[i] if gene_id != "." else '.' for i, gene_id in enumerate(AN_ID_genic)]
    
    VP_ID_tog = cdata["VEP_ensembl_Gene_ID"]      #Gene ID for VEP
    VP_ID_intergenic = [VP_ID_tog[i] if gene_id == "." else "." for i, gene_id in enumerate(AN_ID_genic)]
    VP_ID_genic = [VP_ID_tog[i] if gene_id != "." else '.' for i, gene_id in enumerate(AN_ID_genic)]
    
    # Extracting info from file
    chrs = cdata["chr"] # chromosome number
    pos = cdata["pos"] # SNP position
    ref = cdata["ref"] # reference
    alt = cdata["alt"] # mutation
    rs = cdata["rs_dbSNP151"] # rs ID


    """ Extracting Gene IDs"""
    AN_ID = extract(AN_ID_tog)
    AN_ID_int = extract(AN_ID_intergenic)
    AN_ID_gen = extract(AN_ID_genic)
    
    # Extracting Gene ID in SnpEff
    SN_ID = extract(SN_ID_tog)
    SN_ID_int = extract(SN_ID_intergenic)
    SN_ID_gen = extract(SN_ID_genic)
    
    # Extracting Gene ID in VEP
    VP_ID = extract(VP_ID_tog)
    VP_ID_int = extract(VP_ID_intergenic)
    VP_ID_gen = extract(VP_ID_genic)

    
    # Converting lists to arrays
    AN_ID = np.array(AN_ID, dtype=object)
    SN_ID = np.array(SN_ID, dtype=object)
    VP_ID = np.array(VP_ID, dtype=object)
    
    # Getting rid of repeated Gene IDs
    AN_ID_clean = no_repeats(AN_ID)
    SN_ID_clean = no_repeats(SN_ID)
    VP_ID_clean = no_repeats(VP_ID)
    
    size = len(AN_ID)
    
    
    """ Comparing tools to 'master' annotiation  """
    
    # Merging all Gene IDs
    united = [np.concatenate((AN_ID[i], SN_ID[i], VP_ID[i]), axis=None) for i in range(size)]
    
    # Getting rid of repeats for each SNP
    unique_clean = no_repeats(united)
    
    AN_ID_check = []
    SN_ID_check = []
    VP_ID_check = []
    
    
    # Parsing through dataset
    for i in range(size):
        zize = len(unique_clean[i])  # number of Gene IDs
    
        # Helper function to check presence
        def check_presence(target, source):
            if zize == 1:
                return [zize] if target in source else []
            elif zize > 1:
                return [x + 1 for x in range(zize) if target[x] in source]
    
        AN_ID_check.append(check_presence(unique_clean[i], AN_ID_clean[i]))
        SN_ID_check.append(check_presence(unique_clean[i], SN_ID_clean[i]))
        VP_ID_check.append(check_presence(unique_clean[i], VP_ID_clean[i]))
    
    
    """ Giving table """
    # Creating results table
    d = {"Chr":chrs,
         "Position":pos,
         "Ref":ref,
         "Alt":alt,
         "rs ID":rs,
         "Gene_IDs":unique_clean,
         "ANNOVAR":AN_ID_check,
         "SnpEff":SN_ID_check,
         "VEP":VP_ID_check}
    
    results = pd.DataFrame(d)
    
    """ RESULTS FILE CREATION """
    # Naming output file assuming .gz compression
    # filename_1 = filename.split("/")
    # filename_1 = filename_1[-1]
    # filename_1 = filename_1.split(".")
    # result_filename = filename_1[0]+"_processed."+filename_1[1] #adding the file extension
    
    # Creating a .txt file with results inside
    #results.to_csv(path_or_buf=("/home1/queme/AnnoQ/processed_hrc_12_2019/"+result_filename), sep="\t", index=False)
    #/home1/queme/AnnoQ/processed_hrc_12_2019/  for personal_directory in my_HPC

    tool_agreement_tog = annotation_agreement_rate(unique_clean, AN_ID_check, SN_ID_check, VP_ID_check)
    
    """ SUMMARY FILE CREATION """
    # Summary results
    l_1 = "For Chromosome: " + chr[1] + "\n"
    l_8 = "There are: %i" % size + " SNPs\n\n"
    l_9 = "All tools agree with all GeneIDs: %i (%.2f%%) \n" % (all_agree, all_agree/size*100)
    l_10= """At least two tools agree with all GeneIDs: %i (%.2f%%)
          \t - Annovar: %i \n
          \t - SnpEff: %i  \n
          \t - VEP: %i \n""" % (two_agree, two_agree/size*100, AN_agree_2, SN_agree_2, VP_agree_2)
    l_11= """At least one tool agrees with all GeneIDs: %i (%.2f%%)
          \t - Annovar: %i \n
          \t - SnpEff: %i \n
          \t - VEP: %i \n""" % (one_agree, one_agree/size*100, AN_agree_1, SN_agree_1, VP_agree_1)

    
    # Creating a .txt file with the summary results
    summary_filename = filename_1[0]+"_summary."+filename_1[1]
    f = open("/home1/queme/AnnoQ/processed_hrc_12_2019/"+summary_filename,"w+")
    f.writelines([l_1, l_8, l_9, l_10, l_11, l_12, l_2, l_3, l_4, l_5, l_6, l_7])
    f.close()
