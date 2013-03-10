"""
.. module:: qMS 
   :platform: Any
   :synopsis: A collection of functions for qMS data processing and analysis

.. moduleauthor:: Josh Silverman <josh.silverman@gmail.com>


"""

import os
import sys
import csv
import urllib2
import re
import string
import math
import fileinput
import numpy
import itertools

#qMS utilities

path = os.getcwd()
sys.setrecursionlimit(10000000)

[startindex, pep_seq_index, exp_mr_index, csvopen] = [0, 0, 0, ''];

ribogenes = ['rpsA', 'rpsB','rpsC','rpsD','rpsE','rpsF','rpsG','rpsH','rpsI','rpsJ','rpsK','rpsL','rpsM','rpsN','rpsO','rpsP','rpsQ','rpsR','rpsS','rpsT','rpsU','sra','rplA','rplB','rplC','rplD','rplE','rplF','rplJ','rplL','rplI','rplK','rplM','rplN','rplO','rplP','rplQ','rplR','rplS','rplT','rplU','rplV','rplW','rplX','rplY','rpmA','rpmB','rpmC','rpmD','rpmE','rpmF','rpmG','rpmH','rpmI','rpmJ']

def seqlength(seq):
#returns the theoretical trypsin digest of a protein
    return re.sub(r'(?<=[RK])(?=[^P])','-', seq).split("-")
    
def importFormatTH(filename):
    # Open the data, separate on tabs and store as a list of lists
    file = list(csv.reader(open(filename, 'rU'), delimiter = "\t"))
    # Take the second line which contains the growth rate information
    header = map(float,  file[1][1:])
    # Store the data series as values with the gene names as keys
    gene_data_dict = {item[0]: map(float, item[1:]) for item in file[2:]}
    return gene_data_dict


#fetches protein sequence from uniprot ID
def getsequence(uniprot):
	try:
		urlbase = "http://www.uniprot.org/uniprot/"
		req = urllib2.Request(urlbase + uniprot + ".fasta")
		response = urllib2.urlopen(req)
		lines = [x for i, x in enumerate(response) if i > 0]
		return "".join(map(lambda x: x.rstrip("\n"), lines))
	except urllib2.URLError:
		print "No internet connection or bad Uniprot id"
		
#Returns the MAD of a list
def MAD(list):
	temp_median = numpy.median(list)
	temp_absdevs = map(lambda item: abs(item - temp_median), list)
	return numpy.median(temp_absdevs)
		
#import a.a. mass information
aa_mass_list = list(csv.reader(open('/Users/joshsilverman/Desktop/qMSmodule/aa_masses.csv', 'rU'), dialect = 'excel'))
aa_mass_list = map(lambda x: (x[0], float(x[1])), aa_mass_list)
aa_mass_dict = dict(aa_mass_list)

#imports a list of ribosomal genes
ribogene_file = list(csv.reader(open('/Users/joshsilverman/Desktop/qMSmodule/ribogenes_alt.csv', 'rU'), dialect = 'excel'))
ribogene_list = map(lambda x: x[0], ribogene_file)

#returns the mass of a peptide, taking into account the water weight of the N,C terminal residues
def sequencemass(str, count = 0, mass = 0, waterweight = 18.01528):
	if count == len(str):
		return mass + waterweight
	elif count < len(str):
		if str[count] != ' ':	
			mass += aa_mass_dict[str[count]]
		return sequencemass(str, count + 1, mass)
		
#returns the entropy of a peptide sequence
def entropy(str):
	seq_dict = {}
	for letter in str:
		if letter != ' ':
			seq_dict.setdefault(letter, 0)
			seq_dict[letter] += 1
		else:
			pass
	seqlen = sum(seq_dict.values())
	return [-1 * sum( value/float(seqlen) * math.log(value/float(seqlen)) for value in seq_dict.values()), len(seq_dict.keys()), seqlen]
			
# spectralcount() returns a dictionary with spectral counts for all genes in a given csvfile
#---------------------------------------------------------#

#open csvfile and calculate indices
def initiate(csvfile):
	csvopen = list(csv.reader(open(path + "/" + csvfile, 'rU')))
	startindex = [i for i, x in enumerate(csvopen) for y in x if y.find("prot_hit_num") > -1][0]
	exp_mr_index = [i for i, x in enumerate(csvopen[startindex]) if x.find("pep_exp_mr") > -1][0]
	pep_seq_index = [i for i, x in enumerate(csvopen[startindex]) if x.find("pep_seq") > -1][0]
	return {'exp_mr_index': exp_mr_index, 'pep_seq_index' : pep_seq_index, 'csvopen' : csvopen[startindex+1:-1]}
	
# Intializes and counts in one step, given the filename of the Mascot CSV
def spectralcount(csvfile, cutoff = 4):
	temp_dict = initiate(csvfile)
	
	csvopen = temp_dict['csvopen']
	exp_mr_index = temp_dict['exp_mr_index']
	pep_seq_index = temp_dict['pep_seq_index']
	
	temp_count = 0
	temp_protein_index = 0
	lastgene = ''
	lastacc = ''
	dict = {}
	
	for index, line in enumerate(csvopen):
		tempseq = line[pep_seq_index]
		tempseqmass = sequencemass(tempseq)
		tempgap = abs(tempseqmass - float(line[exp_mr_index]))
		current_protein_index = float(line[0])
		if tempgap < cutoff: temp_count += 1
		if current_protein_index != temp_protein_index:
			temp_protein_index += 1
			tempgene = line[2].split("=")[2].split(" ")[0]
			tempacc = line[1]
			
			if lastgene != '':
				dict[lastgene] = [temp_count, lastacc]
			lastgene = tempgene
			lastacc = tempacc
			temp_count = 0	
	return dict
	
	
# Returns the spectral fraction of the ribosomal genes in the given Mascot CSV
def ribocount(csvfile, cutoff = 4):
	peptides = spectralcount(csvfile, cutoff)
	total = sum(value[0] for value in peptides.values())
	ribototal = sum(peptides[key][0] for key in peptides.keys() if key in ribogene_list)
	return ribototal/float(total)
	
def ribos(csvfiles, cutoff = 4):
	return sum(map(lambda file: ribocount(file, cutoff), csvfiles))/len(csvfiles)
	
def laccount(csvfile, cutoff = 4):
	peptides = spectralcount(csvfile, cutoff)
	total = sum(value[0] for value in peptides.values())
	lactotal = sum(peptides[key][0] for key in peptides.keys() if key == 'lacZ')
	return lactotal/float(total)	
	
def lacos(csvfiles, cutoff = 4):
	return sum(map(lambda file: laccount(file, cutoff), csvfiles))/len(csvfiles)
	
#returns the GO functional classification given a gene name
#---------------------------------------------------------#
#GO num to GO class dictionary
classesfile = open('/Users/joshsilverman/Desktop/qMSmodule/Classes.txt', 'r').readlines()[7:]
go_num2class = {}

for line in classesfile:
    strippedline = string.strip(line)
    if strippedline == "": continue
    splitted = strippedline.partition(" ")
    go_num, go_class = splitted[0], splitted[2]
    go_num2class[go_num] = go_class
    
# Blattner number to GO tag
assignfile = open('/Users/joshsilverman/Desktop/qMSmodule/multifunassignments.txt', 'r').readlines()[4:]
bnum2go_num = {}
for line in assignfile:
    splitted = line.partition("\t")
    bnum =splitted[0]
    if bnum not in bnum2go_num.keys():bnum2go_num[bnum] = []
    bnum2go_num[bnum].append( string.strip(splitted[2]) )

# Blattner number to gene dictionary and the opposite
productfile = open('/Users/joshsilverman/Desktop/qMSmodule/geneproducts.txt','r')
items = map(lambda x: x.split("\t"), productfile)[4:]
gene2bnum = dict(map(lambda list: (list[3],list[1]), items))
bnum2gene = {value: key for (key, value) in gene2bnum.items()}

def sequential(a):
    return map(lambda x: '.'.join(a.split('.')[0:x]), range(1,len(a.split('.')) + 1))
def tags(gene):
    return bnum2go_num[gene2bnum[gene]]
def labels(gene):
    try: 
    	return map(lambda n : map((lambda x: go_num2class[x]), sequential(tags(gene)[n]) ), range(0,len(tags(gene))) )
    except KeyError:
    	print "gene " + gene + " not found in database"
#---------------------------------------------------------#

# Returns two dictionaries, one which maps operons to their genes and a second which maps genes to their operon
def getOperons():
    # Open file
    datapath = '/Users/joshsilverman/Desktop/qMSmodule/OperonSet.txt'
    file = open(datapath, 'rU')
    preLines = list(file)
    
    # Remove preamble lines and those with speculated operons
    lines = [item for item in preLines if item[0] is not '#']
    # Take the operon column from each line
    preOperons = [line.split('\t')[0] for line in lines]
    # Remove operons with numbers in name
    operons = [item for item in preOperons if not re.search(r'[0-9]', item)]
    
    operon_genes, genes_operon = {}, {}
    for item in operons:
        operon_label = item
        # Split on dashes
        first_pass = item.split('-')
        # Expand all genes of form preABCD to preA, preB, ...
        
        # Strange problem in line below in which 'item' is modified to the first piece of splitting at dashes...
        
        second_pass = [[item[0:3] + end for end in item[3:]] if len(item) >= 4 else [item] for item in first_pass]
        # Make flattened list of genes
        genes = list(itertools.chain(*second_pass))
        # Enter gene list into operons_genes with operon label as key
        operon_genes[operon_label] = genes
        # Enter operon label into genes_operons with gene as key
        for gene in genes:
            genes_operon[gene] = operon_label
    return operon_genes, genes_operon

operon_genes, genes_operon = getOperons()


		