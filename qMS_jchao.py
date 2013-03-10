import os
import sys
import csv
import urllib2
import re
import string

#qMS utilities

path = os.getcwd()
sys.setrecursionlimit(100000)

global startindex, pep_seq_index, exp_mr_index, csvopen
[startindex, pep_seq_index, exp_mr_index, csvopen] = [0, 0, 0, ''];

#returns the theoretical trypsin digest of a protein
def seqlength(str):
	return re.sub(r'(?<=[RK])(?=[^P])','-', str).split("-")

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
		
#import a.a. mass information
aa_mass_list = list(csv.reader(open('/Network/Servers/jchao.scripps.edu/Volumes/Lab5/joshuas5/dropbox/aa_masses.csv', 'rU'), dialect = 'excel'))
aa_mass_list = map(lambda x: (x[0], float(x[1])), aa_mass_list)
aa_mass_dict = dict(aa_mass_list)

#returns the mass of a peptide, taking into account the water weight of the N,C terminal residues
def sequencemass(str, count = 0, mass = 0, waterweight = 18.01528):
	if count == len(str):
		return mass + waterweight
	elif count < len(str):
		mass += aa_mass_dict[str[count]]
		return sequencemass(str, count + 1, mass)

#open csvfile and calculate indices
def initializecsv(csvfile):
	csvopen = list(csv.reader(open(path + "/" + csvfile, 'r')))
	startindex = [i for i, x in enumerate(csvopen) for y in x if y.find("prot_hit_num") > -1][0]
	exp_mr_index = [i for i, x in enumerate(csvopen[startindex]) if x.find("pep_exp_mr") > -1][0]
	pep_seq_index = [i for i, x in enumerate(csvopen[startindex]) if x.find("pep_seq") > -1][0]
	return [startindex, exp_mr_index, pep_seq_index, csvopen]
	
#given a Mascot CSV, returns dict with [spectral counts, Uniprot ids] w/ genes as keys
def countspectrals(csvopen, count = 0, proteinindex = 1, dict = {}, line = startindex + 1, lastgene = '', lastacc = '', pep_seq_index = 0, exp_mr_index = 0):
	try:
		temparray = csvopen[line]
		tempseq = temparray[pep_seq_index]
		tempseqmass = sequencemass(tempseq)
		tempgap = abs(tempseqmass - float(temparray[exp_mr_index]))
		if float(temparray[0]) != proteinindex:
			tempid = temparray[2].split("=")[2].split(" ")[0]
			tempacc = temparray[1]
			dict[lastgene] = [count, lastacc]
			if tempgap < 1:
				return countspectrals(csvopen, 1, proteinindex + 1, dict, line+1, tempid, tempacc, pep_seq_index, exp_mr_index)
			else:
				return countspectrals(csvopen, 0, proteinindex + 1, dict, line+1, tempid, tempacc, pep_seq_index, exp_mr_index)
		elif float(temparray[0]) == proteinindex:
			if tempgap < 1:
				return countspectrals(csvopen, count + 1, proteinindex , dict, line+1, lastgene, lastacc, pep_seq_index, exp_mr_index)
			else:
				return countspectrals(csvopen, count, proteinindex , dict, line+1, lastgene, lastacc, pep_seq_index, exp_mr_index)
	except IndexError:
		return dict

#intializes and counts in one step, given the filename of the Mascot csvfile
def spectralcount(csvfile):
	thing = initializecsv(csvfile)
	return countspectrals(thing[3], 0, 1, {}, thing[0]+1, '', '', thing[2], thing[1])
	
#returns the GO functional classification given a gene name
#---------------------------------------------------------#
#GO num to GO class dictionary
classesfile = open('/Network/Servers/jchao.scripps.edu/Volumes/Lab5/joshuas5/dropbox/Classes.txt', 'r').readlines()[7:]
go_num2class = {}

for line in classesfile:
    strippedline = string.strip(line)
    if strippedline == "": continue
    splitted = strippedline.partition(" ")
    go_num, go_class = splitted[0], splitted[2]
    go_num2class[go_num] = go_class
    
#b-number to GO tag
assignfile = open('/Network/Servers/jchao.scripps.edu/Volumes/Lab5/joshuas5/dropbox/multifunassignments.txt', 'r').readlines()[4:]
bnum2go_num = {}
for line in assignfile:
    splitted = line.partition("\t")
    bnum =splitted[0]
    if bnum not in bnum2go_num.keys():bnum2go_num[bnum] = []
    bnum2go_num[bnum].append( string.strip(splitted[2]) )

#b-number to gene product
productfile = open('/Network/Servers/jchao.scripps.edu/Volumes/Lab5/joshuas5/dropbox/geneproducts.txt','r')
items = map(lambda x: x.split("\t"), productfile)[4:]
gene2bnum = dict(map(lambda list: (list[3],list[1]), items))

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

		