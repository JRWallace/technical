#!/usr/bin/env python

import pandas as pd
import glob2
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

#import the mpileup file 
mpileup=pd.concat([pd.read_csv(f, sep='\t') for f in glob2.glob('*_mem_sorted_pileup_coverage_only_for_plotting.txt')])

#import the bedcoverage file, use glob because the exact file name will change depending on quality threshold parameter
bed_coverage=pd.concat([pd.read_csv(f, sep='\t') for f in glob2.glob('*_mem_sorted_bedcoverage_with_header.txt')])

#import the nuc file, use glob because the exact file name will change depending on quality threshold parameter
bed_nuc=pd.concat([pd.read_csv(f, sep='\t') for f in glob2.glob('*_windows_nuc.txt')])

################################################################################

#merge the 2 dataframes based on the values in the first 3 columns - so the chromosome name and window locations
bed_coverage_nuc=bed_nuc.merge(bed_coverage, how='inner', on=['#1_usercol', '2_usercol', '3_usercol'])

#rename the columns because the first column has a hashtag/pound sign, which seems very likely to mess something up later
bed_coverage_nuc.columns=['chromosome','w_start','w_end','pct_at','pct_gc','num_A','num_C', 'num_G', 'num_T', 'num_N', 'num_OTH', 'seq_len', 'coverage']

################################################################################

#kendall correlation 
bed_coverage_nuc_corr_kendall=bed_coverage_nuc.corr(method='kendall')

#spearman correlation 
bed_coverage_nuc_corr_spearman=bed_coverage_nuc.corr(method='spearman')

#pearson correlation 
bed_coverage_nuc_corr_pearson=bed_coverage_nuc.corr(method='pearson')


################################################################################

#make a function to determine if the correlation is weak, strong, or moderate

def coverage_table (row):
    if abs(row['coverage'])<0.3:
        return 'weak'
    if abs(row['coverage'])>0.7:
        return 'strong'
    if abs(row['coverage'])>0.3 and abs(row['coverage'])<0.7:
        return 'moderate'
   
################################################################################
     
#make kendall text output
bed_coverage_nuc_corr_kendall.fillna(0, inplace=True)
bed_coverage_nuc_corr_kendall_table=bed_coverage_nuc_corr_kendall['coverage']
bed_coverage_nuc_corr_kendall_table=bed_coverage_nuc_corr_kendall_table.reset_index()

#apply the function to add a column indicating the strength of the correlation
bed_coverage_nuc_corr_kendall_table['correlation_strength']=bed_coverage_nuc_corr_kendall_table.apply(lambda row: coverage_table(row),axis=1)

#make columns name more accurately reflect their conents
bed_coverage_nuc_corr_kendall_table.columns=['variable','coverage_corr','kendall_corr_strength']

################################################################################

#make spearman text output
bed_coverage_nuc_corr_spearman.fillna(0, inplace=True)
bed_coverage_nuc_corr_spearman_table=bed_coverage_nuc_corr_spearman['coverage']
bed_coverage_nuc_corr_spearman_table=bed_coverage_nuc_corr_spearman_table.reset_index()

#apply the function to add a column indicating the strength of the correlation
bed_coverage_nuc_corr_spearman_table['correlation_strength']=bed_coverage_nuc_corr_spearman_table.apply(lambda row: coverage_table(row),axis=1)

#make columns name more accurately reflect their conents
bed_coverage_nuc_corr_spearman_table.columns=['variable','coverage_corr','spearman_corr_strength']

################################################################################

#make pearson text output
bed_coverage_nuc_corr_pearson.fillna(0, inplace=True)
bed_coverage_nuc_corr_pearson_table=bed_coverage_nuc_corr_pearson['coverage']
bed_coverage_nuc_corr_pearson_table=bed_coverage_nuc_corr_pearson_table.reset_index()

#apply the function to add a column indicating the strength of the correlation
bed_coverage_nuc_corr_pearson_table['correlation_strength']=bed_coverage_nuc_corr_pearson_table.apply(lambda row: coverage_table(row),axis=1)

#make columns name more accurately reflect their conents
bed_coverage_nuc_corr_pearson_table.columns=['variable','coverage_corr','pearson_corr_strength']

################################################################################

with PdfPages('Report.pdf') as pdf:
    mpileup.plot(x='coordinate', y='coverage')
    pdf.savefig()
    plt.close()
    
    plt.table(cellText=bed_coverage_nuc_corr_pearson_table.values, colLabels=bed_coverage_nuc_corr_pearson_table.columns, loc='upper center')
    plt.axis('off')
    pdf.savefig()
    plt.close()
    
    plt.table(cellText=bed_coverage_nuc_corr_spearman_table.values, colLabels=bed_coverage_nuc_corr_spearman_table.columns, loc='upper center')
    plt.axis('off')
    pdf.savefig()
    plt.close()
    
    plt.table(cellText=bed_coverage_nuc_corr_kendall_table.values, colLabels=bed_coverage_nuc_corr_kendall_table.columns, loc='upper center')
    plt.axis('off')
    pdf.savefig()
    plt.close()
    


