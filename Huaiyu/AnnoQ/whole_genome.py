""" LIBRARIES """
import pandas as pd
import numpy as np
import re

""" FUNCTIONS """
# Precompile the regular expression
pattern = re.compile("ENSG...........")

def load_data(file_path):
    columns_needed = ["chr", "pos", "ref", "alt", "rs_dbSNP151", "ANNOVAR_ensembl_Gene_ID","ANNOVAR_ensembl_Closest_gene(intergenic_only)", "SnpEff_ensembl_Gene_ID", "VEP_ensembl_Gene_ID"]
    df = pd.read_csv(file_path, usecols=columns_needed, delimiter="\t", dtype= object, compression='gzip')
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
    for i in range(len(uni_clean)):
        temp = 0  # Temporary counter for the loop.
        des = []  # List to store which arrays agree with uni_clean.
        
        zize = len(uni_clean[i])  # Length of the union set for this iteration.
        
        # Check if AN agrees with uni_clean and update counters.
        if zize == len(AN[i]):
            temp += 1
            des.append('A')
            
        # Check if SN agrees with uni_clean and update counters.
        if zize == len(SN[i]):
            temp += 1
            des.append('S')
            
        # Check if VP agrees with uni_clean and update counters.
        if zize == len(VP[i]):
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
    results.to_csv(path_or_buf=output_path, sep="/t", index=False)

def create_summary_file(chr, size, all_agree, two_agree, AN_agree_2, SN_agree_2, VP_agree_2, one_agree, AN_agree_1, SN_agree_1, VP_agree_1, filename_1):
    """Create and save a summary file."""
    
    summary_lines = [
        f"For Chromosome: {chr[1]}/n",
        f"There are: {size} SNPs/n/n",
        f"All tools agree with all GeneIDs: {all_agree} ({all_agree/size*100:.2f}%)/n",
        f"""At least two tools agree with all GeneIDs: {two_agree} ({two_agree/size*100:.2f}%)
          /t - Annovar: {AN_agree_2}/n
          /t - SnpEff: {SN_agree_2}/n
          /t - VEP: {VP_agree_2}/n""",
        f"""At least one tool agrees with all GeneIDs: {one_agree} ({one_agree/size*100:.2f}%)
          /t - Annovar: {AN_agree_1}/n
          /t - SnpEff: {SN_agree_1}/n
          /t - VEP: {VP_agree_1}/n"""
    ]
    
    summary_filename = f"{filename_1[0]}_summary.{filename_1[1]}"
    output_path = f"/home1/queme/AnnoQ/processed_hrc_12_2019/{summary_filename}"
    
    with open(output_path, "w+") as f:
        f.writelines(summary_lines)

# Example usage
#create_summary_file(chr, size, all_agree, two_agree, AN_agree_2, SN_agree_2, VP_agree_2, one_agree, AN_agree_1, SN_agree_1, VP_agree_1, filename_1)


def annotation_single_to_master(test, master):
    # Initialize lists to store matches and annotation rates.
    match_list = []
    complete_annotation_list = []
    partial_annotation_list = []
    no_annotation_list = []
    snp_annotation_rate = 0  # Initialize the rate counter
    
    # Iterate through both the test and master lists
    for i, (test_item, master_item) in enumerate(zip(test, master)):
        zize = len(test_item)  # Number of Gene IDs in the current test item
        snp_annotation_rate += zize  # Increment the annotation rate counter
        
        # Case 1: If there is only one Gene ID in the test item
        if zize == 1:
            # Check if the test item is in the master item
            is_match = test_item in master_item
            
            # Add to match_list based on the result
            match_list.append([zize] if is_match else [])
            
            # Add to annotation lists, casting boolean to integer (1 if True, 0 if False)
            complete_annotation_list.append(int(is_match))
            #complete_agreement += int(is_match)
            # There can't be partial agreement in a dichotomous outcome, this is just a placeholder
            # if we ever want to include something else and for list ordering purposes
            partial_annotation_list.append(0)
            no_annotation_list.append(int(not is_match))
            #no_agreement += int(not is_match)
        
        # Case 2: If there is more than one Gene ID in the test item
        elif zize > 1:
            # Get indices of matching IDs from the master item
            match_indices = [x + 1 for x, val in enumerate(master_item) if val in test_item]
            match_list.append(match_indices if match_indices else [])
            
            # Initialize both as False at the start of the loop
            complete_match = False
            partial_match = False
            
            # Check for complete match first
            if all(x in test_item for x in master_item):
                complete_match = True
            # Check for partial match only if complete_match is False
            elif any(x in test_item for x in master_item):
                partial_match = True
            
            # Append results to respective lists
            complete_annotation_list.append(int(complete_match))
            partial_annotation_list.append(int(partial_match))
            no_annotation_list.append(int(not partial_match and not complete_match))
            
    complete_agreement = sum(complete_annotation_list)
    partial_agreement = sum(partial_annotation_list)
    no_agreement = sum(no_annotation_list)
            
    # Return the results as a tuple
    return (match_list, 
    complete_annotation_list, 
    partial_annotation_list, 
    no_annotation_list, 
    snp_annotation_rate,
    complete_agreement,
    partial_agreement,
    no_agreement)


def run_and_store_results(func, args, output_names):
    """
    Runs a function with given arguments and stores its multiple outputs in a dictionary.
    
    Parameters:
    func (callable): The function to run.
    args (tuple): The arguments to pass to the function.
    base_name (str): The base name to prepend to each output variable.
    output_names (list of str): The names of the output variables.
    
    Returns:
    dict: A dictionary containing the named outputs.
    """
    # Initialize a dictionary to hold the results
    result_dict = {}
    
    # Function call
    outputs = func(*args)

    # Create dynamic dictionary keys and assign values
    for name, output in zip(output_names, outputs):
        key_name = f"{name}"
        result_dict[key_name] = output
    
    return result_dict




def main():

    # Reading file
    file_path = "https://github.com/quemeb/USC_research/raw/main/Huaiyu/AnnoQ/Test_data.txt.gz"
    cdata = load_data("C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Mi_lab\\AnnoQ\\AnnoQ_data\\21.annotated.snp.gz")
    
    # Selecting data
    AN_ID_genic = cdata["ANNOVAR_ensembl_Gene_ID"]  #Gene ID should be here
    AN_ID_intergenic = cdata["ANNOVAR_ensembl_Closest_gene(intergenic_only)"]    #alternative place for Gene ID
    AN_ID_tog = [AN_ID_intergenic[i] if gene_id.startswith('.') else gene_id for i, gene_id in enumerate(AN_ID_genic)]

    SN_ID_tog = cdata["SnpEff_ensembl_Gene_ID"]   #Gene ID for SnpEff
    SN_ID_intergenic = [SN_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if gene_id.startswith('.')]
    SN_ID_genic = [SN_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]

    
    VP_ID_tog = cdata["VEP_ensembl_Gene_ID"]      #Gene ID for VEP
    VP_ID_intergenic = [VP_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if gene_id.startswith('.')]
    VP_ID_genic = [VP_ID_tog[i] for i, gene_id in enumerate(AN_ID_genic) if not gene_id.startswith('.')]
    
    #Fixing lists after using them as reference
    AN_ID_genic = [gene_id for gene_id in AN_ID_genic if not gene_id.startswith('.')]
    AN_ID_intergenic = [gene_id for gene_id in AN_ID_intergenic if not gene_id.startswith('.')]
    
    # Extracting info from file
    chrs = cdata["chr"] # chromosome number
    pos = cdata["pos"] # SNP position
    ref = cdata["ref"] # reference
    alt = cdata["alt"] # mutation
    rs = cdata["rs_dbSNP151"] # rs ID


    """ Extracting Gene IDs"""
    AN_ID_inter = no_repeats(extract(AN_ID_intergenic))
    AN_ID_genic = no_repeats(extract(AN_ID_genic))
    
    # Extracting Gene ID in SnpEff
    SN_ID_inter = no_repeats(extract(SN_ID_intergenic))
    SN_ID_genic = no_repeats(extract(SN_ID_genic))
    
    # Extracting Gene ID in VEP
    VP_ID_inter = no_repeats(extract(VP_ID_intergenic))
    VP_ID_genic = no_repeats(extract(VP_ID_genic))
    
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


    # Names of the variables you want to create
    output_names = ['match_list', 'complete_annotation_list', 'partial_annotation_list', 'no_annotation_list', 'snp_annotation_rate', 'complete_agreement', 'partial_agreement', 'no_agreement']

    AN_ID_inter_check = run_and_store_results(annotation_single_to_master, (AN_ID_inter, united_unique_inter), output_names)
    SN_ID_inter_check = run_and_store_results(annotation_single_to_master, (SN_ID_inter, united_unique_inter), output_names)
    VP_ID_inter_check = run_and_store_results(annotation_single_to_master, (VP_ID_inter, united_unique_inter), output_names)
    
    AN_ID_genic_check = run_and_store_results(annotation_single_to_master, (AN_ID_genic, united_unique_genic), output_names)
    SN_ID_genic_check = run_and_store_results(annotation_single_to_master, (SN_ID_genic, united_unique_genic), output_names)
    VP_ID_genic_check = run_and_store_results(annotation_single_to_master, (VP_ID_genic, united_unique_genic), output_names)

    
    
    """ Giving table """
    # Creating results table
    #column_names = ["Chr", "Position", "Ref", "Alt", "rs ID", "Gene_IDs", "ANNOVAR", "SnpEff", "VEP"]
    #data_values = [chrs, pos, ref, alt, rs, unique_clean, AN_ID_check, SN_ID_check, VP_ID_check]
    
    #results = pd.DataFrame({col: val for col, val in zip(column_names, data_values)})

    # Creating a DataFrame
    df = pd.DataFrame({
        'Chr': chrs,  # Chromosomes or tools
        'Complete Agreement': [AN_ID_inter_check["complete_agreement"], SN_ID_inter_check["complete_agreement"], VP_ID_inter_check["complete_agreement"],
        'Partial Agreement': [AN_ID_inter_check["complete_agreement"], sn_data[1], vp_data[1]],
        'No Agreement': [AN_ID_inter_check["complete_agreement"], sn_data[2], vp_data[2]],
        'Snp_Annotation_rate': [AN_ID_inter_check["complete_agreement"], sn_data[3], vp_data[3]]
    })
    
    print(df)
    
    """ RESULTS FILE CREATION """

    save_processed_file(file_path, results)

    """  HYPOTHESIS TESTING  """

    tool_agreement_intergenic = annotation_agreement_rate(united_unique_inter, AN_ID_inter, SN_ID_inter, VP_ID_inter)
    tool_agreement_genetic = annotation_agreement_rate(united_unique_genic, AN_ID_genic, SN_ID_genic, VP_ID_genic)
    
    
    
    
    
    """ SUMMARY FILE CREATION """
    # Summary results
    l_1 = "For Chromosome: " + chr[1] + "/n"
    l_2 = "There are: %i" % size + " SNPs/n/n"
    l_3 = f"The agreement between all tools is: {tool_agreement_tog}"


    
    # Creating a .txt file with the summary results
    summary_filename = filename_1[0]+"_summary."+filename_1[1]
    f = open("/home1/queme/AnnoQ/processed_hrc_12_2019/"+summary_filename,"w+")
    f.writelines([l_1, l_8, l_9, l_10, l_11, l_12, l_2, l_3, l_4, l_5, l_6, l_7])
    f.close()
