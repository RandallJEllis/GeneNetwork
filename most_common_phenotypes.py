# -*- coding: utf-8 -*-
"""
Created on Mon Feb 01 11:25:21 2016

@author: ellisrj2
"""
def most_common_phenotypes(wdir):
    import os, os.path
    import pandas 
    import csv
    from collections import Counter
    pathy = os.listdir(wdir)
    phenotypes = []
    records = []
    for f in pathy:
        dataframe = pandas.read_csv(f, skipinitialspace=True, engine='python')
        phenotypes+=dataframe['Phenotype'].tolist()
        records+=dataframe['Record'].tolist()
    record_count = Counter(records)
    most_common_records = record_count.most_common(500)
    most_common_records_records = [x[0] for x in most_common_records] #top 200 records 
    most_common_records_counts = [x[1] for x in most_common_records] #top 200 counts
    most_common_phenotypes = []
   
    for record in most_common_records_records:
        index_record = records.index(record) # index of matched gene in gene_name1 list
        most_common_phenotypes.append(phenotypes[index_record]) #append locations to locations list
    rows = zip(most_common_records_records, most_common_records_counts, most_common_phenotypes)
    
    underscore_index = f.find('_')
    with open(f[:underscore_index+1]+'most_common_phenotypes.csv', 'wb') as thefile:
        writer = csv.writer(thefile)
        writer.writerow(['Record', 'Count', 'Phenotype'])
        for row in rows:
            writer.writerow(row)
            
    #WHEN DOING THIS MANUALLY ACROSS MULTIPLE FOLDERS
    #with open("match_genes_P3_P14_adult.csv", "wb") as thefile:
     #   writer = csv.writer(thefile)
      #  writer.writerow(["P3_P14", "P3_adult", "P14_adult", "All_Three_DBs"])
       # writer.writerows(it.izip_longest(match_genes_P3_P14, match_P3_adult, match_P14_adult,match_genes_all3))