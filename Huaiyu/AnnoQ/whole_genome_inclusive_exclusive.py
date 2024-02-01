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

def total_annotation_agreement(uni_clean, AN, SN, VP):
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


def partial_annotation_aggrement_and_rates(test, master):
    """
    Compares elements of two lists (test and master) to determine their subset and superset relationships.

    Parameters:
        test (list): A list of sets/items to be compared against the master list.
        master (list): A list of sets/items considered as the reference for comparison.

    Returns:
        tuple: Contains the count of occurrences for each relationship type 
        (proper subset, improper subset, not a subset, proper superset, not a superset).
    """
    
    list_length = len(test)

    proper_subset_list = []
    improper_subset_list = []
    not_subset_list = []
    proper_superset_list = []
    not_superset_list = []

    for i in range(list_length):
        test_set = set(test[i])
        master_set = set(master[i])

          # Checking for subset relations
        if test_set < master_set:  # Proper subset
            proper_subset_list.append(1)
        elif test_set == master_set:  # Equal sets (Improper subset)
            improper_subset_list.append(1)
        else:  # Not a subset
            not_subset_list.append(1)

        # Checking for superset relations
        if master_set < test_set:  # Proper superset
            proper_superset_list.append(1)
        else:  # Not a superset
            not_superset_list.append(1)

    # Summing the occurrences for each condition
    proper_subset_count = sum(proper_subset_list)
    improper_subset_count = sum(improper_subset_list)
    not_subset_count = sum(not_subset_list)
    proper_superset_count = sum(proper_superset_list)
    not_superset_count = sum(not_superset_list)

    return (proper_subset_count, improper_subset_count, not_subset_count, 
            proper_superset_count, not_superset_count)


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

def data_summary(tool_agreement, rates, A_S ,V_S, V_A, chrs, sizes):
    dictionary = {"Chr": chrs, 
                "SAV" :tool_agreement['all_agree'], 
                "SA": tool_agreement['SA'], 
                "SV":tool_agreement['SV'], 
                "AV":tool_agreement['AV'], 
                "S":tool_agreement['S'],
                "A":tool_agreement['A'], 
                "V":tool_agreement['V'], 
                "None":tool_agreement['none_agree'], 
                "Total SNPs": sizes,
                # Annotation rates
                "S_rate": rates['S_rate'],
                "A_rate": rates['A_rate'],
                "V_rate": rates['V_rate'],
                # Partial Aggrements Annovar vs SnpEff
                "A_S_proper": A_S['proper'],
                "A_S_improper": A_S['improper'],
                "A_S_disjoint": A_S['disjoint'],
                "A_S_superset": A_S['superset'],
                # Partial aggrement VEP vs SnpEff
                "V_S_proper": V_S['proper'],
                "V_S_improper": V_S['improper'],
                "V_S_disjoint": V_S['disjoint'],
                "V_S_superset": V_S['superset'],
                # Partial Aggrement VEP vs Annovar
                "V_A_proper": V_A['proper'],
                "V_A_improper": V_A['improper'],
                "V_A_disjoint": V_A['disjoint'],
                "V_A_superset": V_A['superset']
                }
    return dictionary 


def data_process(file):

    # Reading file
    # file_path = "https://github.com/quemeb/USC_research/raw/main/Huaiyu/AnnoQ/Test_data.txt.gz"
    # cdata = load_data("/Users/queme/Desktop/USC_research/Huaiyu/AnnoQ/sample_annotations_ch18.txt.gz")
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
                    'superset', 'no_superset', 'snp_annotation_rate']

    A_S_inter_check = run_and_store_results(partial_annotation_aggrement_and_rates, (AN_ID_inter, SN_ID_inter), output_names)
    V_S_inter_check = run_and_store_results(partial_annotation_aggrement_and_rates, (VP_ID_inter, SN_ID_inter), output_names)
    V_A_inter_check = run_and_store_results(partial_annotation_aggrement_and_rates, (VP_ID_inter, AN_ID_inter), output_names)
    
    A_S_genic_check = run_and_store_results(partial_annotation_aggrement_and_rates, (AN_ID_genic, SN_ID_genic), output_names)
    V_S_genic_check = run_and_store_results(partial_annotation_aggrement_and_rates, (VP_ID_genic, SN_ID_genic), output_names)
    V_A_genic_check = run_and_store_results(partial_annotation_aggrement_and_rates, (VP_ID_genic, AN_ID_genic), output_names)

    
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

    #Total Annotaiton aggrement 
    tool_agreement_intergenic = total_annotation_agreement(united_unique_inter, AN_ID_inter, SN_ID_inter, VP_ID_inter)
    tool_agreement_genetic = total_annotation_agreement(united_unique_genic, AN_ID_genic, SN_ID_genic, VP_ID_genic)
    
    # Summary returns
    inter_summary = data_summary(tool_agreement_intergenic, rates_inter , A_S_inter_check, V_S_inter_check, V_A_inter_check, int(chrs[1]), size_inter)
    genic_summary = data_summary(tool_agreement_genetic, rates_genic, A_S_genic_check, V_S_genic_check, V_A_genic_check, int(chrs[1]), size_genic)
    
    return inter_summary, genic_summary


def process_all_files():
    #WINDOWS
    #path = "C:\\Users\\bryan\\OneDrive - University of Southern California\\Research\\Huaiyu Mi\\AnnoQ\\AnnoQ_data\\*.gz"
    
    #MAC
    path = "/Users/queme/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/AnnoQ_data/*.gz"
    
    
    # Step 1: Create empty dictionaries for accumulation
    all_data1 = {}
    all_data2 = {}
    
    for filename in glob.glob(path):
        print(f"Processing file: {filename}")  # Print the current filename
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
    
    #return df_inter, df_genic
    


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
    df_inter_sorted.to_csv('/Users/queme/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/inter_results_new.csv', index=False)
    df_genic_sorted.to_csv('/Users/queme/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/genic_results_new.csv', index=False)


    # End the timer and calculate elapsed time
    elapsed_time = time.time() - start_time

    print("Files saved successfully!")
    print(f"Total time elapsed: {elapsed_time:.2f} seconds")

# This block ensures that the main function is called only when the script is run directly, 
# and not if it's imported as a module in another script.
if __name__ == "__main__":
    main()


