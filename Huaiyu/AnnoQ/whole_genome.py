""" LIBRARIES """
import pandas as pd
import numpy as np
import re
import glob 
import time 

""" FUNCTIONS """
# Precompile the regular expression
pattern = re.compile("ENSG...........")

def load_data(file_path):
    columns_needed = ["chr", "pos", "ref", "alt", "rs_dbSNP151", "ANNOVAR_ensembl_Gene_ID","ANNOVAR_ensembl_Closest_gene(intergenic_only)", "SnpEff_ensembl_Gene_ID", "VEP_ensembl_Gene_ID"]
    df = pd.read_csv(file_path, usecols=columns_needed, delimiter="\t", dtype= object, compression='gzip')
    return df

# Eliminating repeats
def no_repeats(nparray):
    return [np.unique(arr) for arr in nparray]

# Mining gene ID
def extract(lists):
    return[pattern.findall(lst) for lst in lists]

def annotation_agreement_rate(uni_clean, AN, SN, VP):
    # Initialize a dictionary to store various counts.
    counters = {
        'all_agree': 0,
        'two_agree': 0,
        'one_agree': 0,
        'none_agree': 0,
        'SA': 0,
        'SV': 0,
        'AV': 0,
        'A': 0,
        'S': 0,
        'V': 0,
    }
    # Mapping for two-agree scenarios.
    two_agree_mapping = {
        'AS': 'SA',
        'SV': 'SV',
        'AV': 'AV'
    }
    
    # Mapping for one-agree scenarios.
    one_agree_mapping = {
        'A': 'A',
        'S': 'S',
        'V': 'V'
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


def annotation_single_to_master(test, master):
    
    # Get the length of the test and master lists; assume they have the same length
    list_length = len(test)
    snp_annotation_rate = 0
    
    # Initialize lists to store matches and annotation rates.
    match_list = []
    complete_annotation_list = []
    partial_annotation_list = []
    no_annotation_list = []
    snp_annotation_rate = 0  # Initialize the rate counter
    
    # Iterate through both the test and master lists
    for i in range(list_length):
        test_item = test[i]
        master_item = master[i]
        
        zize = len(master_item)  # Number of Gene IDs in the current master item
        rate = len(test_item)   # Number of Gene IDs in the current test item
        snp_annotation_rate += rate  # Increment the annotation rate counter
  
        # Case 1: If Gene ID is empty (this is mostly for VEP in intergenic region)
        if rate == 0: 
            no_annotation_list.append(1)
            partial_annotation_list.append(0)
            complete_annotation_list.append(0)        
        # Case 2: If there is only one Gene ID in the master item (max(test) <= 1)
        elif zize == 1:
            is_match = test_item in master_item  # Check if the test item is in the master item
            
            match_list.append([1] if is_match else [])
            
            # Add to annotation lists, casting boolean to integer (1 if True, 0 if False)
            complete_annotation_list.append(int(is_match))
            partial_annotation_list.append(0)  # Placeholder
            no_annotation_list.append(int(not is_match))
        # Case 3: If there is more than one Gene ID in the master item
        elif zize > 1:
            match_indices = [x + 1 for x in range(len(master_item)) if master_item[x] in test_item]
            match_list.append(match_indices if match_indices else [])
            
            complete_match = False  # Initialize both as False at the start of the loop
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

def data_summary(tool_agreement, S, A, V, chrs, sizes):
    dictionary = {"Chr": chrs, 
                "SAP" :tool_agreement['all_agree'], 
                "SA": tool_agreement['SA'], 
                "SV":tool_agreement['SV'], 
                "AV":tool_agreement['AV'], 
                "S":tool_agreement['S'],
                "A":tool_agreement['A'], 
                "V":tool_agreement['V'], 
                "None":tool_agreement['none_agree'], 
                "Total SNPs":sizes,
                "S_partial": S['partial_agreement'],
                "A_partial": A['partial_agreement'],
                "V_partial": V['partial_agreement'],
                "S_no_agreement": S['no_agreement'],
                "A_no_agreement": A['no_agreement'],
                "V_no_agreement": V['no_agreement'],
                "S_rate": S["snp_annotation_rate"],
                "A_rate": A["snp_annotation_rate"],
                "V_rate": V["snp_annotation_rate"],
                }
    return dictionary 


def data_process(file):

    # Reading file
   # file_path = "https://github.com/quemeb/USC_research/raw/main/Huaiyu/AnnoQ/Test_data.txt.gz"
    cdata = load_data(file)
#    cdata = load_data("C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Mi_lab\\AnnoQ\\Code\\Test_data\\Test_data_set2.txt.gz")
    
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
    
    # # Extracting info from file
    chrs = cdata["chr"] # chromosome number
    # pos = cdata["pos"] # SNP position
    # ref = cdata["ref"] # reference
    # alt = cdata["alt"] # mutation
    # rs = cdata["rs_dbSNP151"] # rs ID


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
    
    # Getting rid of repeats for each SNP
    united_unique_inter = no_repeats(united_inter)
    united_unique_genic = no_repeats(united_genic)
    
    # Converting lists to arrays
    AN_ID_inter = np.array(AN_ID_inter, dtype=object)
    SN_ID_inter = np.array(SN_ID_inter, dtype=object)
    VP_ID_inter = np.array(VP_ID_inter, dtype=object)
    
    AN_ID_genic = np.array(AN_ID_genic, dtype=object)
    SN_ID_genic = np.array(SN_ID_genic, dtype=object)
    VP_ID_genic = np.array(VP_ID_genic, dtype=object)
    
    united_inter = np.array(united_inter, dtype=object)
    united_genic = np.array(united_genic, dtype=object)
    

    # Names of the variables you want to create
    output_names = ['match_list', 'complete_annotation_list', 'partial_annotation_list', 'no_annotation_list', 'snp_annotation_rate', 'complete_agreement', 'partial_agreement', 'no_agreement']

    AN_ID_inter_check = run_and_store_results(annotation_single_to_master, (AN_ID_inter, united_unique_inter), output_names)
    SN_ID_inter_check = run_and_store_results(annotation_single_to_master, (SN_ID_inter, united_unique_inter), output_names)
    VP_ID_inter_check = run_and_store_results(annotation_single_to_master, (VP_ID_inter, united_unique_inter), output_names)
    AN_ID_genic_check = run_and_store_results(annotation_single_to_master, (AN_ID_genic, united_unique_genic), output_names)
    SN_ID_genic_check = run_and_store_results(annotation_single_to_master, (SN_ID_genic, united_unique_genic), output_names)
    VP_ID_genic_check = run_and_store_results(annotation_single_to_master, (VP_ID_genic, united_unique_genic), output_names)


    """  HYPOTHESIS TESTING  """

    tool_agreement_intergenic = annotation_agreement_rate(united_unique_inter, AN_ID_inter, SN_ID_inter, VP_ID_inter)
    tool_agreement_genetic = annotation_agreement_rate(united_unique_genic, AN_ID_genic, SN_ID_genic, VP_ID_genic)
    
    # Summary returns
    inter_summary = data_summary(tool_agreement_intergenic, SN_ID_inter_check, AN_ID_inter_check, VP_ID_inter_check, chrs[1], size_inter)
    genic_summary = data_summary(tool_agreement_genetic, SN_ID_genic_check, AN_ID_genic_check, VP_ID_genic_check, chrs[1], size_genic)
    
    return inter_summary, genic_summary



def process_all_files():
    path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Mi_lab\\AnnoQ\\AnnoQ_data\\*.gz"
    
    # Step 1: Create empty dictionaries for accumulation
    all_data1 = {}
    all_data2 = {}
    
    for filename in glob.glob(path):
        print(f"Processing file: {filename}")  # Print the current filename
        inter, genic = data_process(filename)
        
        # Assuming the chromosome is the key and the data_summary dictionary is the value
        # Update the accumulated dictionaries with new data
        all_data1[inter["Chr"]] = inter
        all_data2[genic["Chr"]] = genic
        
    # Step 3: Convert dictionaries to pandas DataFrame
    df_inter = pd.DataFrame(all_data1).T  # Transpose because keys are chromosomes and values are another dictionary
    df_genic = pd.DataFrame(all_data2).T  # Same reason for Transpose
    
    # Assuming the chromosome numbers are stored in a column named 'Chr', sort by this column
    df_inter_sorted = df_inter.sort_values(by='Chr')
    df_genic_sorted = df_genic.sort_values(by='Chr')
    
    return df_inter_sorted, df_genic_sorted


def main():
    # Start the timer
    start_time = time.time()

    # Process all files and get the sorted dataframes
    df_inter_sorted, df_genic_sorted = process_all_files()
    
    # Save them to CSV
    df_inter_sorted.to_csv('C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Mi_lab\\AnnoQ\\AnnoQ_data\\inter_results.csv', index=False)
    df_genic_sorted.to_csv('C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Mi_lab\\AnnoQ\\AnnoQ_data\\genic_results.csv', index=False)

    # End the timer and calculate elapsed time
    elapsed_time = time.time() - start_time

    print("Files saved successfully!")
    print(f"Total time elapsed: {elapsed_time:.2f} seconds")

# This block ensures that the main function is called only when the script is run directly, 
# and not if it's imported as a module in another script.
if __name__ == "__main__":
    main()


