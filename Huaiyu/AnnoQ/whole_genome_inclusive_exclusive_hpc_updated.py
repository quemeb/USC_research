""" LIBRARIES """
import pandas as pd
import numpy as np
import re
import glob 
import time 


""" FUNCTIONS """
# Precompile the regular expression
pattern = re.compile(r"ENSG\d+",re.IGNORECASE)

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

def complete_annotation_agreement(uni_clean, AN, SN, VP):
    """
    This function compares full aggrement of the unique clean list vs. all
    other individual lists

    Parameters
    ----------
    uni_clean : TYPE
        unique clean list for each SNP.
    AN : TYPE
        Annovar list.
    SN : TYPE
        SnpEff list.
    VP : TYPE
        VEP list.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    
    # Initialize a dictionary to store various counts.
    counters = {
        'all_agree': 0,
        'two_agree': 0,
        'one_agree': 0,
        'none_agree': 0,
        'no_annotation_all': 0,  # Counter for no annotation in all of AN, SN, and VP simultaneously.
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
            
        # Check if all AN, SN, and VP have no annotation simultaneously, regardless of matching to zize
        if len(AN[i]) == 0 and len(SN[i]) == 0 and len(VP[i]) == 0:
            counters['no_annotation_all'] += 1

    return counters  # Return the counters dictionary.


def partial_annotation_agreement(test, master):
    """
    Compares elements of two lists (test and master) to determine their subset (proper or improper), superset, 
    disjoint, and intersection relationships.

    Parameters:
        test (list): A list of sets/items to be compared against the master list. Each element of the list is a set.
        master (list): A list of sets/items considered as the reference for comparison. Each element of the list is a set.

    Returns:
        tuple: Contains the count of occurrences for each relationship type 
        (proper subset, improper subset, disjoint, proper superset, partial overlap).
    """
    
    proper_subset_count = 0
    improper_subset_count = 0
    disjoint_count = 0
    proper_superset_count = 0
    partial_overlap_count = 0
    shouldnt_have_any = 0
    empty_test_count = 0
    empty_master_count = 0
    empty_both = 0

    for test_set, master_set in zip(test, master):
        test_set = set(test_set)
        master_set = set(master_set)
        
        if len(test_set) == 0:
            empty_test_count += 1
        if len(master_set) == 0:
            empty_master_count += 1
        if len(test_set) == 0 and len(master_set) == 0:
            empty_both += 1
            continue  # Skip to the next pair if any is empty


        # Only compare test with master if both are not empty
        if len(test_set) > 0 and len(master_set) > 0:
            # Check if test set is a proper subset of master set
            if test_set < master_set:
                proper_subset_count += 1
            # Check if test set is equal to (improper subset of) master set
            elif test_set == master_set:
                improper_subset_count += 1
            # Check if sets are disjoint (no elements in common)
            elif test_set.isdisjoint(master_set):
                disjoint_count += 1
            # Check if test set is a proper superset of master set
            elif test_set > master_set:
                proper_superset_count += 1
            # Check for partial overlap
            elif test_set & master_set and test_set != master_set:
                partial_overlap_count += 1
            # Check nothing funky is happening
            else:
                shouldnt_have_any += 1

    # Return the counts for each type of relationship and empty counts
    return (proper_subset_count, improper_subset_count, disjoint_count, 
            proper_superset_count, partial_overlap_count, shouldnt_have_any, 
            empty_test_count, empty_master_count, empty_both)


def snp_annotation_rate(test):
    """
    Calculates the SNP annotation rate based on the number of elements in each test item.

    Parameters:
        test (list): A list of sets/items whose annotation rate is to be calculated.

    Returns:
        int: The total SNP annotation rate for the test list.
    """
    snp_annotation_rate = 0

    for item in test:
        snp_annotation_rate += len(set(item))

    return snp_annotation_rate

def looking_for_total_missing(array1, array2, array3):
    """
    Counts the number of empty arrays in each of the three provided numpy arrays.
    
    Parameters:
        array1 (numpy.ndarray): The first numpy array.
        array2 (numpy.ndarray): The second numpy array.
        array3 (numpy.ndarray): The third numpy array.
    
    Returns:
        tuple: A tuple containing the counts of empty arrays in array1, array2, and array3 respectively.
    """
    # Define a function to count empty arrays in a given array
    def count_empty_arrays(arr):
        return sum(1 for item in arr if isinstance(item, np.ndarray) and len(item) == 0)
    
    # Count empty arrays in each array
    count1 = count_empty_arrays(array1)
    count2 = count_empty_arrays(array2)
    count3 = count_empty_arrays(array3)
    
    return (count1, count2, count3)
    

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

def data_summary(tool_agreement, rates, A_S ,V_S, V_A, chrs, sizes, empty):
    dictionary = {"Chr": chrs, 
                "SAV" :tool_agreement['all_agree'], 
                "SA": tool_agreement['SA'], 
                "SV":tool_agreement['SV'], 
                "AV":tool_agreement['AV'], 
                "S":tool_agreement['S'],
                "A":tool_agreement['A'], 
                "V":tool_agreement['V'], 
                "None":tool_agreement['none_agree'], 
                "No_annotations_at_all": tool_agreement['no_annotation_all'],
                "Total_SNPs": sizes,
                # No annotations
                "A_empty": empty[0],
                "S_empty": empty[1],
                "V_empty": empty[2],
                # Annotation rates
                "S_rate": rates['S_rate'],
                "A_rate": rates['A_rate'],
                "V_rate": rates['V_rate'],
                # Partial Aggrements Annovar vs SnpEff
                "A_S_proper": A_S['proper'],
                "A_S_improper": A_S['improper'],
                "A_S_disjoint": A_S['disjoint'],
                "A_S_superset": A_S['superset'],
                "A_S_partial": A_S['partial'],
                "A_S_shouldnt": A_S['shouldnt'],
                "A_S_left_empty": A_S['empty_left'],
                "A_S_right_empty": A_S['empty_right'],
                "A_S_both_empty": A_S['empty_both'],
                # Partial aggrement VEP vs SnpEff
                "V_S_proper": V_S['proper'],
                "V_S_improper": V_S['improper'],
                "V_S_disjoint": V_S['disjoint'],
                "V_S_superset": V_S['superset'],
                "V_S_partial": V_S['partial'],
                "V_S_shouldnt": V_S['shouldnt'],
                "V_S_left_empty": V_S['empty_left'],
                "V_S_right_empty": V_S['empty_right'],
                "V_S_both_empty": V_S['empty_both'],
                # Partial Aggrement VEP vs Annovar
                "V_A_proper": V_A['proper'],
                "V_A_improper": V_A['improper'],
                "V_A_disjoint": V_A['disjoint'],
                "V_A_superset": V_A['superset'],
                "V_A_partial": V_A['partial'],
                "V_A_shouldnt": V_A['shouldnt'],
                "V_A_left_empty": V_A['empty_left'],
                "V_A_right_empty": V_A['empty_right'],
                "V_A_both_empty": V_A['empty_both']
                }
    return dictionary 


def data_process(file):

    # Reading file
    # file_path = "https://github.com/quemeb/USC_research/raw/main/Huaiyu/AnnoQ/Test_data.txt.gz"
    #cdata = load_data("/Users/queme/Desktop/USC_research/Huaiyu/AnnoQ/sample_annotations_ch18.txt.gz")
    cdata = load_data(file)
    
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


    #Rational Counts
    #non_empty_count = [item for item in SN_ID_tog if item != "."]
    #print("Number of non-empty cells:", len(non_empty_count))


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
    output_names = ['proper', 'improper', 'disjoint', 
                    'superset', 'partial', 'shouldnt',
                    'empty_left', 'empty_right', 'empty_both']

    A_S_inter_check = run_and_store_results(partial_annotation_agreement, (AN_ID_inter, SN_ID_inter), output_names)
    V_S_inter_check = run_and_store_results(partial_annotation_agreement, (VP_ID_inter, SN_ID_inter), output_names)
    V_A_inter_check = run_and_store_results(partial_annotation_agreement, (VP_ID_inter, AN_ID_inter), output_names)
    
    A_S_genic_check = run_and_store_results(partial_annotation_agreement, (AN_ID_genic, SN_ID_genic), output_names)
    V_S_genic_check = run_and_store_results(partial_annotation_agreement, (VP_ID_genic, SN_ID_genic), output_names)
    V_A_genic_check = run_and_store_results(partial_annotation_agreement, (VP_ID_genic, AN_ID_genic), output_names)

    
    # Calculating rates and storing them in a dictionary
    rates_inter = {
        "S_rate": snp_annotation_rate(SN_ID_inter),
        "A_rate": snp_annotation_rate(AN_ID_inter),
        "V_rate": snp_annotation_rate(VP_ID_inter)
    }
    rates_genic = {
        "S_rate": snp_annotation_rate(SN_ID_genic),
        "A_rate": snp_annotation_rate(AN_ID_genic),
        "V_rate": snp_annotation_rate(VP_ID_genic)
    }
    
    # total empty arrays
    empty_inter = looking_for_total_missing(AN_ID_inter, SN_ID_inter, VP_ID_inter)
    empty_genic = looking_for_total_missing(AN_ID_genic, SN_ID_genic, VP_ID_genic)

    #Total Annotaiton aggrement 
    tool_agreement_intergenic = complete_annotation_agreement(united_unique_inter, AN_ID_inter, SN_ID_inter, VP_ID_inter)
    tool_agreement_genetic = complete_annotation_agreement(united_unique_genic, AN_ID_genic, SN_ID_genic, VP_ID_genic)
    
    # Summary returns
    inter_summary = data_summary(tool_agreement_intergenic, rates_inter , A_S_inter_check, V_S_inter_check, V_A_inter_check, chrs[1], size_inter,empty_inter)
    genic_summary = data_summary(tool_agreement_genetic, rates_genic, A_S_genic_check, V_S_genic_check, V_A_genic_check, chrs[1], size_genic,empty_genic)
    
    return inter_summary, genic_summary


def process_all_files():
    #WINDOWS
    #path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\AnnoQ_data\\*.gz"
    
    #MAC
    #path = "/Users/queme/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/AnnoQ_data/*.gz"
    
    #HPC
    path = "/home1/queme/AnnoQ/hrc_12_2019/*.gz"
    
    
    # Step 1: Create empty dictionaries for accumulation
    all_data1 = {}
    all_data2 = {}
    
    for filename in glob.glob(path):
        #print(f"Processing file: {filename}")  # Print the current filename
        inter, genic = data_process(filename)
        
        # Assuming the chromosome is the key and the data_summary dictionary is the value
        # Update the accumulated dictionaries with new data
        all_data1[inter['Chr']] = inter
        all_data2[genic['Chr']] = genic
        
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
    # windows platform
    #df_inter_sorted.to_csv('C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\inter_results.csv', index=False)
    #df_genic_sorted.to_csv('C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\genic_results.csv', index=False)

    # Mac Platform
    #df_inter_sorted.to_csv('/Users/queme/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/inter_results_new.csv', index=False)
    #df_genic_sorted.to_csv('/Users/queme/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/genic_results_new.csv', index=False)

    # HPC 
    df_inter_sorted.to_csv("/home1/queme/AnnoQ/hrc_12_2019_subsets_counts/inter_results_new_updated.csv", index=False)
    df_genic_sorted.to_csv("/home1/queme/AnnoQ/hrc_12_2019_subsets_counts/genic_results_new_updated.csv", index=False)



    # End the timer and calculate elapsed time
    elapsed_time = time.time() - start_time

    #print("Files saved successfully!")
    print(f"Total time elapsed: {elapsed_time:.2f} seconds")

# This block ensures that the main function is called only when the script is run directly, 
# and not if it's imported as a module in another script.
if __name__ == "__main__":
    main()


