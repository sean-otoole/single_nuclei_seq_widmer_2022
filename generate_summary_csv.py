#!/usr/bin/env python
# coding: utf-8

# In[27]:

import csv
import os

def sample_summary(csv_directory_list, csv_output_path, chemistry_dictionary):
    
    ### constructs a header
    
    first_csv = csv_directory_list[0]

    with open(first_csv) as csvfile:
        reader = csv.reader(csvfile, delimiter=',')
        header = next(reader)
        header.insert(0, 'Sample')
        header.append('Chemistry')
        
    ### create the metrics list and add a header

    summary_list = []
    summary_list.append(header)
    
    for sample in csv_directory_list:
        with open(sample, "rt") as infile:
           reader = csv.reader(infile, delimiter = ',')
           next(reader, None)  # skip the headers
           for row in reader:
                sample_name = sample.split('/')[-3]
                row_list = row
                row_list.insert(0, sample_name)
                row_list.append(chemistry_dictionary[sample_name])
                summary_list.append(row)
    
    ### now write a list of lists to your file
    
    if os.path.isfile(csv_output_path) == True:
        os.remove(csv_output_path)  ### this will clear the previous file

    with open(csv_output_path, "a", newline="") as f:
        writer = csv.writer(f)
        writer.writerows(summary_list)
        f.close()