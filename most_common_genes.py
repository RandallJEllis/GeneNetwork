# -*- coding: utf-8 -*-
"""
Created on Thu Feb 04 14:24:30 2016
@author: ellisrj2
"""
#def xls_to_csv(xls1):
#    import xlrd
#    import unicodecsv
#    import csv
#    
#    #CONVERT FROM XLS TO CSV
#    wb1 = xlrd.open_workbook(xls1)
#    sh1 = wb1.sheet_by_name('Sheet1')
#    csv1 = str(xls1[:-4]) + '.csv'
#    csv1_obj = open(csv1, 'wb')
#    wr1 = unicodecsv.writer(csv1_obj, quoting=csv.QUOTE_ALL, encoding='utf-8')
#    for rownum in xrange(sh1.nrows):
#        wr1.writerow(sh1.row_values(rownum))
#    return csv1
    
def remove_nonsig(wdir):
    import os, os.path
    import pandas 
    path = os.listdir(wdir)
    for f,i in zip(path,range(len(path))):
        df = pandas.read_excel(f, header=7, skipfooter=3)
        column_names = list(df.columns.values)
        df = df[df[column_names[11]] < .005] # Sample p(rho) or what GN has
        df = df[abs(df[column_names[9]]) >= 0.60] #Sample rho or what GN has
        df.to_excel('sig_' + f, index=False)

def most_common_genes(wdir):
    import os, os.path
    import pandas 
    import csv
    import numpy as np
    import itertools
    from collections import Counter
    path = os.listdir(wdir)
    genes = []
    locations = []
    column_names = [] #for writing CSV column names
    for f,i in zip(path,range(len(path))):
        dot_index = f.rfind('.')
        column_names.append(f[:dot_index])
        dataframe = pandas.read_csv(f, skipinitialspace=True, engine='python')
        genes.append(dataframe['Symbol'].tolist())
        locations.append(dataframe['Location (Chr: Mb)'].tolist())
    allgenes = itertools.chain(*genes)
    allgenes = list(allgenes)
    alllocations = itertools.chain(*locations)
    alllocations = list(alllocations)
    gene_count = Counter(allgenes)
    most_common_genes = gene_count.most_common(500)
    most_common_genes_genes = [x[0] for x in most_common_genes] #top 200 records 
    most_common_genes_counts = [x[1] for x in most_common_genes] #top 200 counts
    most_common_locations = []
    individual_db_counts = np.zeros((500,len(path)))
        
    for i in range(0, len(most_common_genes)):
        for j in range(0, len(path)):
            if most_common_genes_genes[i] in genes[j]:
                indices = [z for z,x in enumerate(genes[j]) if x == most_common_genes_genes[i]] 
                individual_db_counts[i,j] += len(indices)
    
    #create an array for each gene of its individual db counts
    individual_db_counts_columns = []
    for i in range(len(column_names)):
        individual_db_counts_columns.append(individual_db_counts[0:len(most_common_genes), i])
    
    for gene in most_common_genes_genes:
        index_gene = allgenes.index(gene) # index of matched gene in gene_name1 list
        most_common_locations.append(alllocations[index_gene]) #append locations to locations list
        
    real_occurrences = np.zeros((500,1)) #initialize real_occurrences list   
    for i in range(0,500):
        for x in individual_db_counts[i,]: #cycling through each row with 10 columns (DBs)
            if x > 0:
                real_occurrences[i] += 1 #count +1 for each nonzero value
            
    real_occurrences = [int(x) for x in real_occurrences] #convert real_occurrences elements to intsF
    
    rows = zip(most_common_genes_genes, most_common_locations, most_common_genes_counts, real_occurrences)
    
    rows = [list(t) for t in rows]
    
    #add individual db counts to rows
    for i in range(len(rows)):
        for j in range(len(column_names)):
            rows[i].append(individual_db_counts_columns[j][i])
    
    columns = ['Gene', 'Location', 'Count', 'Real Count']
    for column_name in column_names:
        columns.append(column_name)
        
    underscore_index = f.find('_')
    no_first_underscore = f[underscore_index+1:]
    end_of_name = no_first_underscore.find('_')
    name_file = f[underscore_index+1:underscore_index+end_of_name+1]
    with open(name_file+'_'+'most_common_genes.csv', 'wb') as thefile:
        writer = csv.writer(thefile)
        writer.writerow(columns)
        for row in rows:
            writer.writerow(row)