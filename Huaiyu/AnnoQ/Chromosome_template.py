import pandas as pd
import numpy as np
import re
import sys

########################################
# Functions
pattern = re.complie("ENSG...........")

# Eliminating repeats
def no_repeats(nparray):
    place = []
    for i in range(0,size):
        temp = []
        temp = np.unique(nparray[i])
        place.append(temp)
    return(place)

# Mining gene ID
def extract(lists):
    place =[]
    for i in range(0,size):
        temp = []
        temp = re.findall("ENSG...........", lists[i])
        place.append(temp)
    return(place)

# Looking for differences
def difference(array1,array2):
    place =[]
    for i in range(0,size):
        temp = []
        temp = np.setdiff1d(array1[i], array2[i], assume_unique = True)
        place.append(temp)
    return(place)

# Counting
def counter(arrays):
    place = 0
    for i in range(0,size):
        if(arrays[i].size > 0):
            place = place + 1
    return place
########################################

# file location
filename = sys.argv[1]

# reading file
cdata = pd.read_csv(filename, delimiter="\t", dtype= object, compression="gzip")

# Selecting data
AN_ID_1 = cdata["ANNOVAR_ensembl_Gene_ID"]  #Gene ID should be here
AN_ID_2 = cdata["ANNOVAR_ensembl_Closest_gene(intergenic_only)"]    #alternative place for Gene ID

SN_ID_1 = cdata["SnpEff_ensembl_Gene_ID"]   #Gene ID for SnpEff

VP_ID_1 = cdata["VEP_ensembl_Gene_ID"]      #Gene ID for VEP

# Finding length of dataset
size = len(AN_ID_1)

# Freeing memory
del cdata

# Finding the Gene ID in ANNOVAR
AN_ID_tog = []
for i in range(0,size):
    if(AN_ID_1[i] == "."):
        AN_ID_tog.append(AN_ID_2[i])
    else:
        AN_ID_tog.append(AN_ID_1[i])

# extracting Gene ID in ANNOVAR
AN_ID = []
AN_ID = extract(AN_ID_tog)

# extracting Gene ID in SnpEff
SN_ID = []
SN_ID = extract(SN_ID_1)

# extracting Gene ID in VEP
VP_ID = []
VP_ID = extract(VP_ID_1)

# converting lists to arrays
AN_ID = np.array(AN_ID, dtype=object)
SN_ID = np.array(SN_ID, dtype=object)
VP_ID = np.array(VP_ID, dtype=object)

# Getting rid of repeated Gene IDs
AN_ID_clean = no_repeats(AN_ID)
SN_ID_clean = no_repeats(SN_ID)
VP_ID_clean = no_repeats(VP_ID)

# merging all Gene IDs
united = []
for i in range(0,size):
  temp = []
  temp = AN_ID[i] + SN_ID[i] + VP_ID[i]
  united.append(temp)

# getting rid of repeats for each SNP
unique_clean = no_repeats(united)

""" Comparing Unique Gene IDs  """
# Creating empty arrays
AN_ID_check = []
SN_ID_check = []
VP_ID_check = []

for i in range(0,size):
  zize = len(unique_clean[i])
  if zize == 1:
    if unique_clean[i] in AN_ID_clean[i]:
      AN_ID_check.append([zize]) #else AN_ID_check.append("")
    if unique_clean[i] in SN_ID_clean[i]:
      SN_ID_check.append([zize]) #else SN_ID_check.append("")
    if unique_clean[i] in VP_ID_clean[i]:
      VP_ID_check.append([zize]) #else VP_ID_check.append("")
  elif zize > 1:
    place_1 = []
    place_2 = []
    place_3 = []
    for x in range(0,zize):
      if unique_clean[i][x] in AN_ID_clean[i]:
        place_1.append(x+1)
      if unique_clean[i][x] in SN_ID_clean[i]:
        place_2.append(x+1)
      if unique_clean[i][x] in VP_ID_clean[i]:
        place_3.append(x+1)
    AN_ID_check.append(place_1)
    SN_ID_check.append(place_2)
    VP_ID_check.append(place_3)

# creating results table
results = pd.DataFrame(data = [unique_clean, AN_ID_check, SN_ID_check, VP_ID_check],
                       index = ["Gene_IDs", "ANNOVAR","SnpEff", "VEP"])
results = results.swapaxes("index", "columns")

# naming output file
filename_1 = filename.split(".")
result_filename = filename_1[0]+"_processed."+filename_1[-1]

# Creating a .txt file with results inside
results.to_csv(path_or_buf=('/Users/bryanqueme/OneDrive - University of Southern California/Research/Mi Lab/AnnoQ/'+result_filename), sep="\t", index=False)
#/home1/queme/AnnoQ/processed_hrc_12_2019  for personal_directory in my_HPC
