import os
import sys
import numpy as np
import pandas as pd
import matplotlib
import plotly
import argparse
import xlrd



def load_key(filename):
	try:
		df=pd.io.excel.read_excel(filename,header=0)
	except xlrd.biffh.XLRDError:
		df=pd.read_csv(filename,header=0,sep='\t')
	locus_dictionary=pd.Series(df['Genome'].values,index=df['Locus Tag']).to_dict()
	geneID_dictionary=pd.Series(df['Genome'].values,index=df['Gene ID']).to_dict()
	
	return locus_dictionary,geneID_dictionary
	
def main(argv):
	gyra_dict={'83':['F','L','Y'],'87':['G','N','Y']}
	rpsl_dict={'43':['N'],'88':['R','Q'],'90':['S']}
	rpob_dict={'504':['F'],'511':['F'],'512':['F','Y','P'],'513':['L','H'],'516':['G','Y'],'522':['F','Y'],'526':['Y','L','P'],'529':['S','L','H'],'531':['F'],'564':['L'],'572':['F']}
	base_dir = os.path.dirname(os.path.realpath(__file__))
	parser = argparse.ArgumentParser(description = 'Parser for Pairwise_stats')
	parser.add_argument('-i','--infile',  nargs=1, help='An output file from pairwise_mutations.py', required=True)
	parser.add_argument('-k','--keyfile',  type=str, nargs=1, help='A keyfile for assigning genomes to sequences',required=True)
	parser.add_argument('-p','--positions',  type=int, nargs='+', help='A space seperated integer list of mutation positions',required=True)
	parser.add_argument('-g','--gene_type',  type=str, nargs=1, help='Gene type',required=True)
	args=vars(parser.parse_args())
	tax_list=['Phylum','Class','Order','Family','Genus','Species','Subspecies']
	header=['uniqueID']
	header.extend(args['positions'])
	header.extend(tax_list)
	locus_dictionary,geneID_dictionary=load_key(str(args['keyfile'][0]))
	df=pd.io.parsers.read_csv(str(args['infile'][0]),sep='\t',header=None,names=header,index_col=False,low_memory=True)
	print df
	z= dict(locus_dictionary.items()+geneID_dictionary.items())
	print z
	df["uniqueID"].replace(z,inplace=True)
	group = df.groupby('uniqueID')
	
	
	#Assign the proper mutation dictionary
	if str(args['gene_type'][0]) == 'gyra':
		mut_dict=gyra_dict
	if str(args['gene_type'][0]) == 'rpsl':
		mut_dict=rpsl_dict
	if str(args['gene_type'][0]) == 'rpob':
		mut_dict=rpob_dict
	to_boxplot=pd.DataFrame()
	#Sum Number of total sequences per group
	df['Seq_count']=df['uniqueID'].map(lambda x: 1 )
	#to_boxplot.append(df.groupby(by=['uniqueID'])['Seq_count'].sum())
	columns=[]
	columns.append(df.groupby(by=['uniqueID'])['Seq_count'].sum())
	
	
	#Sum Number of resistant sequences per group
	for pos in args['positions']:
		df[str(pos)+"_resistant"]=df[pos].map(lambda x: 1 if x in mut_dict[str(pos)] else 0)

	
	for pos in args['positions']:
		columns.append(df.groupby(by=['uniqueID'])[str(pos)+"_resistant"].sum())
	final_df=pd.DataFrame(columns)
	final_df=final_df.transpose()
	final_df['Gene']=final_df['Seq_count'].map(lambda x: str(args['gene_type'][0]))
	final_df=final_df[final_df['Seq_count']>9]
	print final_df.to_csv(sep='\t')
	#Generate Box plots with plotly
		
	#print df
	


if __name__=="__main__":
	main(sys.argv)