from collections import defaultdict
import brewer2mpl
import pandas as pd
import sys
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import scipy.spatial.distance as distance
import scipy.cluster.hierarchy as sch
from hcluster import pdist, linkage, dendrogram
import Levenshtein as L
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist, squareform
import csv
from Bio import Phylo
from cStringIO import StringIO



def tree(): return defaultdict(tree)

def tree_add(t, path):
  for node in path:
	node = node.replace(" ","_")
	t = t[node]

def pprint_tree(tree_instance):
	def dicts(t): return {k: dicts(t[k]) for k in t}
	pprint(dicts(tree_instance))

def csv_to_tree(input):
	t = tree()
	for row in csv.reader(input, quotechar='\''):
		tree_add(t, row)
	return t

def tree_to_newick(root):
	items = []
	for k in root.iterkeys():
		s = ''
		if len(root[k].keys()) > 0:
			sub_tree = tree_to_newick(root[k])
			if sub_tree != '':
				s += '(' + sub_tree + ')'
		s += k
		items.append(s)
	return ','.join(items)

def csv_to_weightless_newick(input):
	t = csv_to_tree(input)
	#pprint_tree(t)
	return tree_to_newick(t)



def heatmap(dm):
	
	
	D1 = squareform(pdist(dm, metric='euclidean'))
	D2 = squareform(pdist(dm.T, metric='euclidean'))
	
	f = plt.figure(figsize=(8, 8))

	# add first dendrogram
	ax1 = f.add_axes([0.09, 0.1, 0.2, 0.6])
	Y = linkage(D1, method='complete')
	Z1 = dendrogram(Y, orientation='right')
	ax1.set_xticks([])
	ax1.set_yticks([])

	# add second dendrogram
	ax2 = f.add_axes([0.3, 0.71, 0.6, 0.2])
	Y = linkage(D2, method='complete')
	Z2 = dendrogram(Y)
	ax2.set_xticks([])
	ax2.set_yticks([])

	# add matrix plot
	axmatrix = f.add_axes([0.3, 0.1, 0.6, 0.6])
	idx1 = Z1['leaves']
	idx2 = Z2['leaves']
	D = D1[idx1, :]
	D = D[:, idx2]
	im = axmatrix.matshow(D, aspect='auto', origin='lower', cmap='hot')
	axmatrix.set_xticks([])
	axmatrix.set_yticks([])
	
	return {'ordered' : D, 'rorder' : Z1['leaves'], 'corder' : Z2['leaves']}


def form_df(filepath):
	df = pd.read_table(filepath,sep='\t',names=['gi','83','87','Phylum','Class','Order','Family','Genus','Species','Subspecies'])
	df['83res']=df['83'].map(lambda x: 1 if (x=='L' or x=='F' or x=='Y') else 0)
	df['87res']=df['87'].map(lambda x: 1 if (x=='N' or x=='Y' or x=='G') else 0)
	df['Totres']=df['87res']+df['83res']
	
	score_dict={}
	for tup in df[['Phylum','Class','Order','Family','Genus','Species','Subspecies']].itertuples():
		
		for num in range(1,len(tup)):
			#print df['83res'][tup[0]]
			if tup[num] != 'nan':
				if str(tup[num])+"83" in score_dict:
					score_dict[str(tup[num])+"83"]+=int(df['83res'][tup[0]])
				else:
					score_dict[str(tup[num])+"83"]=int(df['83res'][tup[0]])
				if str(tup[num])+"87" in score_dict:
					score_dict[str(tup[num])+"87"]+=int(df['87res'][tup[0]])
				else:
					score_dict[str(tup[num])+"87"]=int(df['87res'][tup[0]])
				if str(tup[num])+"Tot" in score_dict:
					score_dict[str(tup[num])+"Tot"]+=int(df['Totres'][tup[0]])
				else:
					score_dict[str(tup[num])+"Tot"]=int(df['Totres'][tup[0]])
		
	for tup in df[['Phylum','Class','Order','Family','Genus','Species','Subspecies']]:	
		
		#print tup
		df[tup+"83"]=df[tup].map(lambda x: 0 if pd.isnull(x) else score_dict[str(x)+"83"])
		df[tup+"87"]=df[tup].map(lambda x:0 if pd.isnull(x) else score_dict[str(x)+"87"])
		df[tup+"Tot"]=df[tup].map(lambda x:0 if pd.isnull(x) else score_dict[str(x)+"Tot"])
	
	return df
	
def dendogram(dataframe,taxlevel):
	#TODO
	return 0
def clean_axis(ax):
	"""Remove ticks, tick labels, and frame from axis"""
	ax.get_xaxis().set_ticks([])
	ax.get_yaxis().set_ticks([])
	for sp in ax.spines.values():
		sp.set_visible(False)

if __name__ == '__main__':
	df=form_df(sys.argv[1])
	
	try:
		sys.argv[2]
	except NameError:
		tax_level='Phylum'
	else:
		tax_level= sys.argv[2]
	
	tax_dict={'Phylum':['Phylum'],'Class':['Phylum','Class'],'Order':['Phylum','Class','Order'],'Family':['Phylum','Class','Order','Family'],'Genus':['Phylum','Class','Order','Family','Genus'],'Species':['Phylum','Class','Order','Family','Genus','Species']}


# font size for figures
#rcParams.update({'font.size': 16})
# Arial font
#rc('font',**{'family':'sans-serif','sans-serif':['Arial']})

#amino_df=df[['83res','87res','Totres']]
#X = df[['Phylum','Class','Order','Family','Genus','Species','Subspecies']].T.values
#Y=pdist(X)
#Z=linkage(Y)
#dendogram(Z)


#df['FullTax']=pd.concat([df.ix[:,'Phylum':'Species']])
#dataMatrix<-as.matrix(df[['Phylum83','Phylum87','PhylumTot']])
	#print df.as_matrix(['Phylum83','Phylum87','PhylumTot'])
#heatmap(df.as_matrix(['Phylum83','Phylum87','PhylumTot']))
#print df.ix[:,'Phylum':'Species']
	#print tax_dict[tax_level]
	df[tax_dict[tax_level]].to_csv("gyrA_histo_results.csv",header=False,index=False)
	df.to_csv("gyrA_results.csv")
	tree_list=[]
	with open('gyrA_histo_results.csv','rb') as csvfile:
		spamreader = csv.reader(csvfile, delimiter='\t', quotechar='|')
		for row in spamreader:
			tree_list.append(', '.join(row))
	
	newick= csv_to_weightless_newick(tree_list)
	print str(newick)
	tree=Phylo.read(StringIO(str(newick)),"newick")
	Phylo.draw(tree)
	leaves_list=tree.get_terminals()

	order_list=[]
	for clade in leaves_list:
		order_list.append(clade.name)
	
		
	
	#Heatmap stuff
	taxcount_dict={}
	for index, row in df.iterrows():
		if row[str(tax_level)] in taxcount_dict:
			taxcount_dict[row[str(tax_level)]]+=1
		else:
			taxcount_dict[row[str(tax_level)]]=1
	
	
	
	newdf=df[[str(tax_level),str(tax_level)+'83',str(tax_level)+'87',str(tax_level)+'Tot']].drop_duplicates(subset=[str(tax_level)])
	newdf[str(tax_level)+'Totseq']=newdf[str(tax_level)].apply(lambda x:taxcount_dict[x])
	newdf[str(tax_level)+'83norm']=newdf[str(tax_level)+'83']/newdf[str(tax_level)+'Totseq']*100
	newdf[str(tax_level)+'87norm']=newdf[str(tax_level)+'87']/newdf[str(tax_level)+'Totseq']*100
	newdf[str(tax_level)+'Totnorm']=newdf[str(tax_level)+'Tot']/newdf[str(tax_level)+'Totseq']*100
	
	newdf=newdf.set_index(str(tax_level))
	print order_list
	newdf=newdf.loc[order_list]
	phylo_labels= newdf.index.tolist()
	print newdf
	data=newdf.fillna(0).as_matrix([str(tax_level)+'83norm',str(tax_level)+'87norm',str(tax_level)+'Totnorm'])
	
	
	
	#print phylo_labels.shape
	#plt.figure(1)
	#plt.subplot(121)
	plt.imshow(data,interpolation='none', aspect=3./data.shape[0])
	plt.xticks(range(3),['83 Resistant','87 Resistant', 'Total Resistant'])
	plt.yticks(range(data.shape[0]),phylo_labels)
	plt.jet()
	plt.colorbar()
	
	
	#TODO:this is the list order of the terminal nodes of the tree, need to sort dataframe using this!
	
	
	plt.show()
	
