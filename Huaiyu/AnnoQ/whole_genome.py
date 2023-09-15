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

def save_processed_file(filename, results, output_directory="/home1/queme/AnnoQ/processed_hrc_12_2019/"):
    """
    Saves the processed results to a new file.

    Parameters:
    - filename: str
        The original filename.
    - results: pd.DataFrame
        The processed results to be saved.
    - output_directory: str, optional
        The directory where the processed file will be saved. Default is "/home1/queme/AnnoQ/processed_hrc_12_2019/"
    - Example usage
        save_processed_file("your/original/filename.txt", results_dataframe)
    """
    
    # Extracting the base filename and extension
    base_filename = filename.split("/")[-1]
    name_parts = base_filename.split(".")
    
    # Creating the processed filename
    result_filename = f"{name_parts[0]}_processed.{name_parts[1]}"
    
    # Saving the results to a .txt file
    output_path = f"{output_directory}{result_filename}"
    results.to_csv(path_or_buf=output_path, sep="\t", index=False)

def create_summary_file(chr, size, all_agree, two_agree, AN_agree_2, SN_agree_2, VP_agree_2, one_agree, AN_agree_1, SN_agree_1, VP_agree_1, filename_1):
    """Create and save a summary file."""
    
    summary_lines = [
        f"For Chromosome: {chr[1]}\n",
        f"There are: {size} SNPs\n\n",
        f"All tools agree with all GeneIDs: {all_agree} ({all_agree/size*100:.2f}%)\n",
        f"""At least two tools agree with all GeneIDs: {two_agree} ({two_agree/size*100:.2f}%)
          \t - Annovar: {AN_agree_2}\n
          \t - SnpEff: {SN_agree_2}\n
          \t - VEP: {VP_agree_2}\n""",
        f"""At least one tool agrees with all GeneIDs: {one_agree} ({one_agree/size*100:.2f}%)
          \t - Annovar: {AN_agree_1}\n
          \t - SnpEff: {SN_agree_1}\n
          \t - VEP: {VP_agree_1}\n"""
    ]
    
    summary_filename = f"{filename_1[0]}_summary.{filename_1[1]}"
    output_path = f"/home1/queme/AnnoQ/processed_hrc_12_2019/{summary_filename}"
    
    with open(output_path, "w+") as f:
        f.writelines(summary_lines)

# Example usage
#create_summary_file(chr, size, all_agree, two_agree, AN_agree_2, SN_agree_2, VP_agree_2, one_agree, AN_agree_1, SN_agree_1, VP_agree_1, filename_1)


def annotation_single_to_master(test, master):
    match_list = []
    complete_annotation_list = []
    partial_annotation_list = []
    no_annotation_list = []
    SNP_annotation_rate = 0
    
    # Loop through the entire source list
    for i in range(len(master)):
        zize = len(test[i])  # Number of Gene IDs in the current source element
        SNP_annotation_rate = SNP_annotation_rate + zize # total number of annotations 

        # Check for matches and append to temp accordingly
        if zize == 1:
            match_list.apprend([zize] if test[i] in master[i] else [])
            if all(x in test[i] for x in master[i]):
                complete_annotation_list.append(1)
                partial_annotation_list.append(0)
            elif 
                complete_annotation_list.append(1 if all(x in test[i] for x in master[i]) else 0)
                partial_annotation_list.append(1 if any(x in test[i] for x in master[i]) else 0)
            
        elif zize > 1:
            match_list.append([x+1 for x in range(len(test[i])) if master[i][x] in test[i]])
            complete_agree_list.append(1 if all(x in test[i] for x in master[i]) else 0)
            partial_agree_list.append(1 if any(x in test[i] for x in master[i]) else 0)

    return match_list, complete_annotation_list, partial_annotation_list, SNP_annotation_rate 

AN_ID_inter_check = annotation_single_to_master(AN_ID_inter, united_unique_inter)
AN_ID_genic_check = annotation_single_to_master(AN_ID_genic, united_unique_genic)



def main():

    # Reading file
    file_path = "https://github.com/quemeb/USC_research/raw/main/Huaiyu/AnnoQ/Test_data.txt.gz"
    cdata = load_data(file_path)
    
    # Selecting data
    AN_ID_genic = cdata["ANNOVAR_ensembl_Gene_ID"]  #Gene ID should be here
    AN_ID_intergenic = cdata["ANNOVAR_ensembl_Closest_gene(intergenic_only)"]    #alternative place for Gene ID
    AN_ID_tog = [AN_ID_intergenic[i] if gene_id == "." else gene_id for i, gene_id in enumerate(AN_ID_genic)]
    AN_ID_genic = [gene_id for gene_id in AN_ID_genic if gene_id != '.']
    AN_ID_intergenic = [gene_id for gene_id in AN_ID_intergenic if not gene_id.startswith('.')]


    SN_ID_tog = cdata["SnpEff_ensembl_Gene_ID"]   #Gene ID for SnpEff
    SN_ID_intergenic = [SN_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]
    SN_ID_genic = [SN_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]

    
    VP_ID_tog = cdata["VEP_ensembl_Gene_ID"]      #Gene ID for VEP
    VP_ID_intergenic = [VP_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]
    VP_ID_genic = [VP_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]
    
    # Extracting info from file
    chrs = cdata["chr"] # chromosome number
    pos = cdata["pos"] # SNP position
    ref = cdata["ref"] # reference
    alt = cdata["alt"] # mutation
    rs = cdata["rs_dbSNP151"] # rs ID


    """ Extracting Gene IDs"""
    AN_ID_inter = extract(AN_ID_intergenic)
    AN_ID_genic = extract(AN_ID_genic)
    
    # Extracting Gene ID in SnpEff
    SN_ID_inter = extract(SN_ID_intergenic)
    SN_ID_genic = extract(SN_ID_genic)
    
    # Extracting Gene ID in VEP
    VP_ID_inter = extract(VP_ID_intergenic)
    VP_ID_genic = extract(VP_ID_genic)
    
    size_inter = len(AN_ID_inter)
    size_genic = len(AN_ID_genic)
    
    
    """ ---------------  Comparing tools to 'master' annotiation  ------------------"""
    
    # Merging all Gene IDs
    united_inter = [np.concatenate((AN_ID_inter[i], SN_ID_inter[i], VP_ID_inter[i]), axis=None) for i in range(size_inter)]
    united_genic = [np.concatenate((AN_ID_genic[i], SN_ID_genic[i], VP_ID_genic[i]), axis=None) for i in range(size_genic)]
    
    # Converting lists to arrays
    AN_ID_inter = np.array(AN_ID_inter, dtype=object)
    SN_ID_inter = np.array(SN_ID_inter, dtype=object)
    VP_ID_inter = np.array(VP_ID_inter, dtype=object)
    
    AN_ID_genic = np.array(AN_ID_genic, dtype=object)
    SN_ID_genic = np.array(SN_ID_genic, dtype=object)
    VP_ID_genic = np.array(VP_ID_genic, dtype=object)
    
    united_inter = np.array(united_inter, dtype=object)
    united_genic = np.array(united_genic, dtype=object)
    
    # Getting rid of repeats for each SNP
    united_unique_inter = no_repeats(united_inter)
    united_unique_genic = no_repeats(united_genic)

    AN_ID_inter_check = annotation_single_to_master(AN_ID_inter, united_unique_inter)
    SN_ID_inter_check = annotation_single_to_master(SN_ID_inter, united_unique_inter)
    VP_ID_inter_check = annotation_single_to_master(VP_ID_inter, united_unique_inter)
    
    AN_ID_genic_check = annotation_single_to_master(AN_ID_genic, united_unique_genic)
    SN_ID_genic_check = annotation_single_to_master(SN_ID_genic, united_unique_genic)
    VP_ID_genic_check = annotation_single_to_master(VP_ID_genic, united_unique_genic)
    
    
    
    """ Giving table """
    # Creating results table
    column_names = ["Chr", "Position", "Ref", "Alt", "rs ID", "Gene_IDs", "ANNOVAR", "SnpEff", "VEP"]
    data_values = [chrs, pos, ref, alt, rs, unique_clean, AN_ID_check, SN_ID_check, VP_ID_check]
    
    results = pd.DataFrame({col: val for col, val in zip(column_names, data_values)})

    
    """ RESULTS FILE CREATION """

    save_processed_file(file_path, results)

    """  HYPOTHESIS TESTING  """

    tool_agreement_tog = annotation_agreement_rate(unique_clean, AN_ID_check, SN_ID_check, VP_ID_check)
    tool_agreement_intergenic = annotation_agreement_rate(unique_clean, AN, SN, VP)
    tool_agreement_genetic = annotation_agreement_rate(unique_clean, )
    
    """ SUMMARY FILE CREATION """
    # Summary results
    l_1 = "For Chromosome: " + chr[1] + "\n"
    l_2 = "There are: %i" % size + " SNPs\n\n"
    l_3 = f"The agreement between all tools is: {tool_agreement_tog}"


    
    # Creating a .txt file with the summary results
    summary_filename = filename_1[0]+"_summary."+filename_1[1]
    f = open("/home1/queme/AnnoQ/processed_hrc_12_2019/"+summary_filename,"w+")
    f.writelines([l_1, l_8, l_9, l_10, l_11, l_12, l_2, l_3, l_4, l_5, l_6, l_7])
    f.close()
