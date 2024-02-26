#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Feb 16 11:08:21 2024

@author: queme
"""

import glob
import matplotlib.pyplot as plt

def create_histograms_from_files(folder_path):
    # Find all text files in the specified folder
    file_pattern = f"{folder_path}/*.txt"
    for filename in glob.glob(file_pattern):
        # Read the numbers from the file
        with open(filename, 'r') as file:
            numbers = [int(line.strip()) for line in file]
        
        # Create a histogram for the numbers
        plt.figure(figsize=(10, 6))  # Optional: Adjust figure size
        plt.hist(numbers, bins='auto', alpha=0.7, color='blue', edgecolor='black')

        # Extract the last part of the filename for the title
        file_title = filename.split("/")[-1]  # This splits the filename by '/' and takes the last element
        plt.title(f'Histogram for {file_title}')
        
        plt.xlabel('Value')
        plt.ylabel('Frequency')
        
        # Save the histogram to a file or display it
        # To save, uncomment the next line and replace 'path/to/save' with your desired save location
        # plt.savefig(f'path/to/save/histogram_{filename.split("/")[-1].replace(".txt", ".png")}')
        
        # To display, uncomment the next line (Note: This will pause the script for each plot if running interactively)
        plt.show()

# Example usage
folder_path = "/Users/queme/OneDrive - University of Southern California/Research/Huaiyu Mi/AnnoQ/Raw Data Comparisons/Intergenic_distances"
create_histograms_from_files(folder_path)