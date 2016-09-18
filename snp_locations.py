# -*- coding: utf-8 -*-
"""
Created on Fri Apr 08 23:41:29 2016
@author: ellisrj2

The input for this script should be the output of a Marker Regression analysis from GeneNetwork. 

GeneNetwork curates these as txt files; all you have to do is open them in any spreadsheet program and save it as a CSV.

The purpose of this script is to gather all relevant SNPs with likelihood ratio statistics (LRS) greater than the suggestive value (p=0.63, [ADD BASIC STATISTICAL SUMMARY WITH LINK]) and format them to be readable by the genome assembly conversion utility from the UCSC Genome Browser.
"""

import pandas

def snp_locations(genenetwork_marker_regression_map_snp_csv, output_filename):
    
    #create empty list of snp locations
    snp_locations = []
    
    #read snp table
    snp_table = pandas.read_csv(genenetwork_marker_regression_map_snp_csv, skipinitialspace=True, engine='python')
    
    #name columns appropriately
    snp_table = snp_table.rename(columns={'Unnamed: 1': 'chromosome', 'Unnamed: 2': 'megabase', 'Unnamed: 3': 'locus' })
    
    #calculate Suggestive LRS value
    suggestive_LRS = snp_table.columns[0]
    idx_before_number = suggestive_LRS.rfind(' ')
    suggestive_LRS_value = suggestive_LRS[idx_before_number+1:]
    
    #convert all LRS values to floats to section values greater than suggestive LRS
    for i in range(3,len(snp_table[str(suggestive_LRS)])):
        snp_table[suggestive_LRS][i] = float(snp_table[suggestive_LRS][i])
    
    #retain all entries with LRS values greater than suggestive value
    snp_table_sig = snp_table[snp_table[suggestive_LRS] > float(suggestive_LRS_value)]
        
    #prepare chromosomes to be read by UCSC 
    chromosomes = snp_table_sig['chromosome'][3:]
    chromosomes = [str(chromosome) for chromosome in chromosomes]
    for i in range(len(chromosomes)):
        chromosomes[i] = 'chr' + chromosomes[i] + ':'
    
    #prepare megabases to be read by UCSC
    megabases = snp_table_sig['megabase'][3:]
    megabases = [str(megabase) for megabase in megabases]
    for i in range(len(megabases)):
        megabases[i] = megabases[i].translate(None, '.')
        megabases[i] = megabases[i] + '-' + str(int(megabases[i]) + 1)
        
    #loop through chromosomes and megabases to create entries and append to list
    for i in range(len(chromosomes)):
        entry = chromosomes[i] + megabases[i]
        snp_locations.append(entry)
        
    with open(output_filename, 'w') as f:
        for snp in snp_locations:
            f.write(snp + '\n')
    