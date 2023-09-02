""" LIBRARIES """
import pandas as pd
import numpy as np
import re

""" FUNCTIONS """
# Precompile the regular expression
pattern = re.compile("ENSG...........")

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

""" DATA PROCESSING """

""" File Input & Data Mining """
# Reading file
url = "https://github.com/quemeb/USC_research/raw/main/Huaiyu/AnnoQ/Test_data.txt.gz"
cdata = pd.read_csv(url, delimiter="\t", dtype= object, compression='gzip')

# Selecting data
AN_ID_1 = cdata["ANNOVAR_ensembl_Gene_ID"]  #Gene ID should be here
AN_ID_intergenic = cdata["ANNOVAR_ensembl_Closest_gene(intergenic_only)"]    #alternative place for Gene ID
AN_ID_tog = [AN_ID_intergenic[i] if gene_id == "." else gene_id for i, gene_id in enumerate(AN_ID_1)]

SN_ID_1 = cdata["SnpEff_ensembl_Gene_ID"]   #Gene ID for SnpEff
SN_ID_intergenic = [SN_ID_1[i] if gene_id == '.' else '.' for i, gene_id in enumerate(AN_ID_1)]

VP_ID_1 = cdata["VEP_ensembl_Gene_ID"]      #Gene ID for VEP
VP_ID_intergenic = [VP_ID_1[i] if gene_id == "." else "." for i, gene_id in enumerate(AN_ID_1)]

chrs = cdata["chr"] # chromosome number
pos = cdata["pos"] # SNP position
ref = cdata["ref"] # reference
alt = cdata["alt"] # mutation
rs = cdata["rs_dbSNP151"] # rs ID

# Extracting Gene ID in ANNOVAR
AN_ID = extract(AN_ID_tog)

# Extracting Gene ID in SnpEff
SN_ID = extract(SN_ID_1)

# Extracting Gene ID in VEP
VP_ID = extract(VP_ID_1)

# Converting lists to arrays
AN_ID = np.array(AN_ID, dtype=object)
SN_ID = np.array(SN_ID, dtype=object)
VP_ID = np.array(VP_ID, dtype=object)

# Getting rid of repeated Gene IDs
AN_ID_clean = no_repeats(AN_ID)
SN_ID_clean = no_repeats(SN_ID)
VP_ID_clean = no_repeats(VP_ID)

size = len(AN_ID)


"""  Summing total annotations in each tool  """
#This part will help us see the SNPs annotation rate by tool overall

# Getting rid of repeated Gene IDs
AN_ID_clean = no_repeats(AN_ID)
SN_ID_clean = no_repeats(SN_ID)
VP_ID_clean = no_repeats(VP_ID)

# Counting 
AN_total = sum(len(AN_ID_entry) for AN_ID_entry in AN_ID_clean)
SN_total = sum(len(SN_ID_entry) for SN_ID_entry in SN_ID_clean)
VP_total = sum(len(VP_ID_entry) for VP_ID_entry in VP_ID_clean)


""" Comparing Unique Gene IDs  """

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


""" COMPARISON RESULTS ARRAY TABLE """
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
filename_1 = filename.split("/")
filename_1 = filename_1[-1]
filename_1 = filename_1.split(".")
result_filename = filename_1[0]+"_processed."+filename_1[1] #adding the file extension

# Creating a .txt file with results inside
#results.to_csv(path_or_buf=("/home1/queme/AnnoQ/processed_hrc_12_2019/"+result_filename), sep="\t", index=False)
#/home1/queme/AnnoQ/processed_hrc_12_2019/  for personal_directory in my_HPC

""" QUANTIFYING DIFFERENCES """
# Comparing SN vs AN
dif_SN_AN = difference(SN_ID_clean, AN_ID_clean)
dif_SN_AN_sum = counter(dif_SN_AN)

int_SN_AN = inter(SN_ID_clean, AN_ID_clean)
int_SN_AN_sum = counter(int_SN_AN)

int_SN_VP = inter(SN_ID_clean, VP_ID_clean)
int_SN_VP_sum = counter(int_SN_VP)

int_AN_VP = inter(AN_ID_clean, VP_ID_clean)
int_AN_VP_sum = counter(int_AN_VP)



# Comparing SN vs VP
dif_SN_VP = []
dif_SN_VP = difference(SN_ID_clean, VP_ID_clean)
dif_SN_VP_sum = counter(dif_SN_VP)

# Comparing AN vs SN
dif_AN_SN = []
dif_AN_SN = difference(AN_ID_clean, SN_ID_clean)
dif_AN_SN_sum = counter(dif_AN_SN)

# Comparing AN vs VP
dif_AN_VP = []
dif_AN_VP = difference(AN_ID_clean, VP_ID_clean)
dif_AN_VP_sum = counter(dif_AN_VP)

# Comparing VP vs SN
dif_VP_SN = []
dif_VP_SN = difference(VP_ID_clean, SN_ID_clean)
dif_VP_SN_sum = counter(dif_VP_SN)

# Comparing VP vs AN
dif_VP_AN = []
dif_VP_AN = difference(VP_ID_clean, AN_ID_clean)
dif_VP_AN_sum = counter(dif_VP_AN)


""" COMPARING TOOLS TO ALL GENEIDS """
all_agree = 0 #with GeneID union
two_agree = 0 #with GeneID union
one_agree = 0 #with GeneID union
none_agree= 0 #with GeneID union
AN_SN_agree_2 = 0 #missing from VP
SN_VP_agree_2 = 0 #missing form AN
AN_VP_agree_2 = 0 #missing from SN
AN_agree_1 = 0 #only AN
SN_agree_1 = 0 #only SN
VP_agree_1 = 0 #only VP

for i in range(0,size):
  zize = len(unique_clean[i])
  temp = 0
  des = []
  if zize == len(AN_ID_check[i]):
    temp += 1
    des.append("A")
  if zize == len(SN_ID_check[i]):
    temp += 1
    des.append("S")
  if zize == len(VP_ID_check[i]):
    temp += 1
    des.append("V")

# assigning tool aggrements... 
  if temp == 3:             # this means all should agree
    all_agree += 1

  elif temp == 2:           # exactly two agree
    two_agree += 1
    if des == ["A","S"]:
      AN_SN_agree_2 += 1
    if des == ["S","V"]:
      SN_VP_agree_2 += 1
    if des == ["A","V"]:
      AN_VP_agree_2 += 1

  elif temp == 1:           # exactly one agrees
    one_agree += 1
    if des == ["A"]:
      AN_agree_1 += 1
    if des == ["B"]:
      SN_agree_1 += 1
    if des == ["C"]:
      VP_agree_1 += 1

  elif temp == 0:           # no tool has all GeneIDs
    none_agree+= 1          # can hardly see it happening

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
l_12= "No tool agrees with all GeneIDs:  %i (%.2f%%)\n\n" % (none_agree, none_agree/size*100)
l_2 = "%i \t in ANNOVAR \t but not in \t SnpEff  \n" % dif_AN_SN_sum
l_3 = "%i \t in ANNOVAR \t but not in \t VEP \n" % dif_AN_VP_sum
l_4 = "%i \t in SnpEff \t but not in \t ANNOVAR \n" % dif_SN_AN_sum
l_5 = "%i \t in SnpEff \t but not in \t VEP \n" % dif_SN_VP_sum
l_6 = "%i \t in VEP    \t but not in \t ANNOVAR \n" % dif_VP_AN_sum
l_7 = "%i \t in VEP    \t but not in \t SnpEff \n" % dif_VP_SN_sum

# Creating a .txt file with the summary results
summary_filename = filename_1[0]+"_summary."+filename_1[1]
f = open("/home1/queme/AnnoQ/processed_hrc_12_2019/"+summary_filename,"w+")
f.writelines([l_1, l_8, l_9, l_10, l_11, l_12, l_2, l_3, l_4, l_5, l_6, l_7])
f.close()
