#!/usr/bin/python3

# parse VCF and GFF files and calculate SNP distance matrix within a specified region
# usage: vcf_parse_distance.py gff_file.gff vcf_file.vcf output_file.tsv chromosome start_position end_position
import re
import sys
arguments = sys.argv

gff_file = [1]
vcf = [2]
outfile = open([3], "w")
chromosome = string([4])
start_pos = int([5])
end_pos = int([6])

def getDistMatrix(acclist):
	"""Generates distance matrix of appropriate size (based on no. samples in VCF)
	need as many lists in distmatrix as there are accessions and as many zeros in each list as there are accessions"""
	noAccs = len(acclist)
	distmatrix = [['NA' for x in range(noAccs)] for y in range(noAccs)]

	return distmatrix

def addToDistMatrix(acc1, acc2, distance, acclist, distmatrix):
	"""Enters distance into correct position of distance matrix
	based on position of accs in acclist"""

	# get indices of acc1 and acc2
	acc1i = acclist.index(acc1)
	acc2i = acclist.index(acc2)
	distmatrix[acc2i][acc1i] = distance

	#print(acc1, acc2, acc1i, acc2i)

	return distmatrix

def getAcc2(acc1, acclist):
	
	for acc in acclist:
		if acc != acc1:
			acc2index = acclist.index(acc)
			# check if there is already an entry in distance matrix for this comparison
			# if yes skip over
			# or if empty matrix cells = 0 need to record which comparisons have been done
			acc2 = acc

	return acc2

def getAccList(accdict):
	for item in accdict.items():
		acclist.append(item[0])

	return acclist

def calculateDistance(acc1, acc2, acclist, accdict):
	"""Matches keys of vardict then does a string comparison
	to assess distance between accessions"""
		
	subdict1 = accdict.get(acc1)
	subdict2 = accdict.get(acc2)
	
	distance = 0
	# test if length of both dictionaries the same?
	# need look up table of indicies and accessions
	# think of way to not do comparisons twice

	for key in subdict1.keys():
		allele1 = subdict1.get(key)
		allele2 = subdict2.get(key)
		
		split1=re.split('[/|]', allele1)
		split2=re.split('[/|]', allele2)

		foundEqual=False
		for s1 in split1:
			for s2 in split2:
				if s1 == s2:
					# WHAT IS THE COMMA??? Commas only used in VCFs as separators... some kind of error if they are getting called as alleles?
					# need to check if missign data (.) is because whole gene is missing or just a few positions
					foundEqual=True

		if not foundEqual:
			distance += 1
	return distance

def getAccDict(variation, accdict):
	
	"""Creates a dictionary with acceesion as they key,
	values are a dictionary with chrm, position tuple as key
	and the variant as the value."""
	for acc in list(variation[4].keys()):
		vardict = {}
		if acc in accdict:
			vardict= accdict[acc]
		
		
		vardict[(variation[0], variation[1])] = variation[4][acc]
		
		accdict[acc] = vardict

	
	return accdict

def isexon(genedict, position):
	if (position >= start_pos and position <= end_pos):
		return True
	else:
		return False

def getGenedict(gff_file, seqname):
	"""Create dictionary of gene intervals"""
	genedict = {}

	with open(gff_file) as gff:
		for line in gff:
			if not(line.startswith("#") ): #make sure this is not a comment line
				splitLine = line.split("\t", -1) # split line and store columns as list
				start = int(splitLine[3]) # take the fourth column and convert to a number
				end = int(splitLine[4])    # take the fifth column and convert to a number

				# now, check if the chromosome is correct and that the exon is within our interval. Also check that this is an exon
				if splitLine[0] == seqname and start >= start_pos and end <= end_pos:
					gene = splitLine[8].split("Parent=")[1].split(";")[0]
					gene = gene.rstrip() # rstrip removes potential newline characters.
					interval = (start, end)

					if  not gene in genedict:
						exonset = {interval}
						genedict[gene] = exonset
					else:
						genedict[gene].add(interval)
	return genedict

def convertAlleleToBase(genotypes, refAllele, altAllele):
	"""Convert X/X code for ref/alt allele to the base. 
	Code = index of altAllele +1 (0 = refAlelle)
	Alleles in altAllele comma separated eg GGTCA,GGTCG"""
	splitAlt = altAllele.split(",")
	for gt_key in genotypes.keys():
		gt = genotypes[gt_key]
		for i in re.split("[/|]", gt):
			if i == '0':
				gt = gt.replace(i, refAllele)
			elif i not in ['1', '2', '3', '4']: 
				# account for weird gt values?
				pass
			else:
				gt = gt.replace(i, altAllele[int(i)-1])
		genotypes[gt_key] = gt

	return genotypes

def recordSampleNames(line):
	splitLine = line.split("\t")
	samples = []
	for i in range(9,len(splitLine)):
		samples.append(splitLine[i].strip())
	return samples

def getGT_index(format_column):
	"""check the format colum, at which position the GT field is. Return the index"""
	splitLine = format_column.split(":")
	for i in range(0, len(splitLine)):
		if splitLine[i] == "GT":
			return i

def getGenotypes(splitLine, gt_index, sampleNames):
	genotypes = []
	for i in range(9, len(splitLine)):
		sample_field = splitLine[i]
		split_col = sample_field.split(":")
		genotypes.append(split_col[gt_index])

	genotype_dict = {}
	for i in range(0,len(sampleNames)):
		genotype_dict[sampleNames[i]] = genotypes[i]
	return genotype_dict

def getvars(line, sampleNames, genedict):
	"""Extract required information from each position (line)"""
	splitLine = line.split("\t", -1)
	chromosome = splitLine[0]
	position = int(splitLine[1])
	refAllele = splitLine[3]
	altAllele = splitLine[4]
	gt_index = getGT_index(splitLine[8])

	if isexon(genedict, position):
		genotypes = getGenotypes(splitLine, gt_index, sampleNames)
		convertAlleleToBase(genotypes, refAllele, altAllele)
		return(chromosome, position, refAllele, altAllele, genotypes)

with open(vcf) as infile:
	sampleNames = []
	accdict = {}
	acclist = []
	seqname = chromosome
	genedict = getGenedict(gff_file, seqname)
	for line in infile:
		if line.startswith('#CHROM'):
			sampleNames = recordSampleNames(line)
		elif not line.startswith('#'):
			variation = getvars(line, sampleNames, genedict)
			if variation:
				accdict = getAccDict(variation, accdict)
	
	acclist = getAccList(accdict)
	distmatrix = getDistMatrix(acclist) 
	with outfile as record:
		record.write("Accessions" + "\t" + '\t'.join([str(x) for x in acclist]) + "\n")	
		for acc1 in acclist:
			print("Acc1 " + acc1)
			acc1i = acclist.index(acc1)
			record.write(acc1)
			for acc2 in acclist:
				distance = calculateDistance(acc1, acc2, acclist, accdict)
				record.write("\t" + str(distance))
				addToDistMatrix(acc1, acc2, distance, acclist, distmatrix)
			#record.write(acc1 + "\t" + '\t'.join([str(x) for x in distmatrix[acc1i]]) + "\n")
			record.write("\n")
			