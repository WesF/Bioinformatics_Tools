import numpy as np
import pandas as pd
import sys
from Bio import SeqIO
import math
import os.path
import argparse
import os
import Bio
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqRecord import SeqRecord
import re
import csv
import ast
from Bio import pairwise2
from Bio.SubsMat import MatrixInfo as matlist
from Bio import Entrez

class CmdLineParser(object):
	parsed=None
	
	def __init__(self):
		parser = argparse.ArgumentParser(description = 'Pairwise aligns sequences with a reference using Emboss NEEDLE and checks for at given positions for mutations. ')
		
		parser.add_argument('-i','--infile',metavar='Infile', type=str, nargs=1, help='FASTA file handle of sequences to align', dest='infile')
		parser.add_argument('-p','--positions', metavar='Positions', type=int, nargs='+', help='A space seperated integer list of mutation positions')
		parser.add_argument('-o','--outfile', metavar='Outfile', type=str, nargs='?', help='File handle for output, defaults to output.txt', default='./output.txt', dest='outfile')
		parser.add_argument('-r','--refseq', metavar='RefSeq', type=str, nargs=1, help='File handle for Reference Sequence')
		parser.add_argument('-d','--db', metavar='DB', type=str, nargs=1, help='File handle for BLAST database')
		
	def parse(self, state):
		if not self.parsed:
			self.parsed = parser.parse_args()

def makedb(dbpath):
	blastDB = "makeblastdb -in "+dbpath+" -input_type fasta -dbtype prot -out blastdb"
	os.system(str(blastDB))
	
def pairwise(refseq,qseq):
	gap_open=-10
	gap_extend=-0.5
	matrix=matlist.blosum62
	new_query=Bio.Seq.Seq(str(qseq).translate(None,"*").translate(None,"-"),Bio.Alphabet.IUPAC.IUPACProtein())
	alignments = pairwise2.align.globalds(refseq.upper(),new_query.upper(),matrix,gap_open,gap_extend)
	return alignments

def get_residues(alignments,res_no):
	
	j=1
	for i,res in enumerate(alignments[0][0]):
		if res!='-':
			if j==int(res_no):
				print i
				return alignments[0][1][i]
			j+=1
	i=0
def get_taxonomy(org_name):
	handle = Entrez.esearch(db="Taxonomy", term=str(org_name))
	record = Entrez.read(handle)
	try:
		taxonomy = Entrez.efetch(db="Taxonomy", id=record["IdList"][0], retmode="xml")
	except IndexError :
		print "There was an error with Entrez efetch!"
		return ['error','error','error','error','error','error','error','error']
	else:
		taxrecords = Entrez.read(taxonomy)
		tax_list = taxrecords[0]["Lineage"].split('; ')
		return tax_list

def blastp(base_dir,fas, db):
	blast_in = (fas)
	blast_path = os.path.join(base_dir, 'blast.xml')
	blastp_cline = NcbiblastpCommandline(query=blast_in, db=db, outfmt=5, out=blast_path, max_target_seqs=1)
	os.system(str(blastp_cline))
	print blastp_cline
	return blast_path


def main(argv):
	base_dir = os.path.dirname(os.path.realpath(__file__))
	
	parser = argparse.ArgumentParser(description = 'Pairwise aligns sequences with a reference using Emboss NEEDLE and checks for at given positions for mutations. ')
	
	parser.add_argument('-i','--infile',  nargs=1, help='FASTA file handle of sequences to align', dest='infile',required=True)
	parser.add_argument('-p','--positions',  type=int, nargs='+', help='A space seperated integer list of mutation positions',required=True)
	parser.add_argument('-o','--outfile',   nargs=1, help='File handle for output, defaults to output.txt', default='./output.txt', dest='outfile',required=True)
	parser.add_argument('-r','--refseq',   nargs=1, help='File handle for Reference Sequence',required=True)
	parser.add_argument('-d','--db',   nargs=1, help='File handle for BLAST database',required=True)
	
	args=vars(parser.parse_args())
	print args
	
	ref_seq= SeqIO.read(str(args['refseq'][0]),'fasta')
	db_path= makedb(args['db'][0])
	#TODO:initialize DataFrame
	df=pd.DataFrame()
	to_read=open(args['infile'][0],"r+") 
	column_labels=['gi']
	print column_labels
	column_labels.extend(args['positions'])
	column_labels.extend(['organism','Kingdom','Phylum','Class','Order','Family','Genus','Species','Subspecies'])
	df=pd.DataFrame(columns=column_labels)
	print "these are the labels"
	print df
	name = 0
	with open('testing_testing.txt','w') as file_:
		for seq in SeqIO.FastaIO.FastaIterator(to_read):
			to_df=[seq.id]
			alignments=pairwise(ref_seq.seq,seq.seq)
			for position in args['positions']:
				to_df.append(get_residues(alignments,position))
			#Run BLAST
			with open(os.path.join(base_dir, 'temp.fasta'), "w+") as temp:
				temp.write(">"+seq.id+"\n")
				temp.write(str(seq.seq)+"\n")
			
		
			blast=blastp(base_dir,temp.name,'blastdb') 
			temp.close()
		
			with open(blast) as to_read:
				blast_record = NCBIXML.read(to_read)
		
			#Get Top BLAST HIT organism name for taxonomy
			Entrez.email = 'weslfield@gmail.com'
			try:
				org_name=re.split(r'[\[\]]+', blast_record.alignments[0].title)
			except IndexError:
				continue
		
			tax_list=get_taxonomy(org_name[1])
			print org_name[1]+'\t'+tax_list[-1]
			del tax_list[0:2]
			to_df.extend(tax_list)
			file_.write('\t'.join(map(str,to_df)))
			file_.write("\n")
		
			#add to_df to DataFrame as a row
			df_series=pd.Series(data=to_df,name=name)
			name+=1
			df.append(df_series,ignore_index=True)
			os.remove(os.path.join(base_dir, 'temp.fasta'))
		file_.close()
		df.to_csv(path_or_buf='./results.csv',na_rep='0')
		#Send DataFrame to csv

if __name__ == "__main__":
	main(sys.argv)