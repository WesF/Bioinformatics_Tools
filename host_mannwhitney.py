from scipy.stats import mannwhitneyu, wilcoxon,ranksums
import numpy as np
import pandas as pd
import sys
import os
import argparse



def load_df(filename,sheetname):
	df=pd.io.excel.read_excel(filename,header=0,index_col=None,sheetname=sheetname)
	return df


def main(argv):
	base_dir = os.path.dirname(os.path.realpath(__file__))
	parser = argparse.ArgumentParser(description = 'Parser for Pairwise_stats')
	parser.add_argument('-i','--infile',  nargs=1, help='An output file from pairwise_mutations.py', required=True)
	parser.add_argument('-g','--gene_type',  type=str, nargs=1, help='Gene type',required=True)
	args=vars(parser.parse_args())
	
	df = load_df(str(args['infile'][0]),str(args['gene_type'][0]))
	#non_host_series = df['Total_norm'][df['Association']=='nonhost-associated']
	host_list=['Gut','Vagina','Tongue','Hoatzin','Dog','Insects']
	nonhost_list=['Soil','Aquatic','Waste']
	
	for host in host_list:
		for nonhost in nonhost_list:
			host_series = df['Total_norm'][df['Environment']==str(host)]
			nonhost_series=df['Total_norm'][df['Environment']==str(nonhost)]
			u=mannwhitneyu(host_series.values,nonhost_series.values)
			#print len(host_series.values)
			print host+' vs '+nonhost+'\t'+str(u)
	
	
if __name__=="__main__":
	main(sys.argv)