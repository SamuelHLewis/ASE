#!/usr/bin/env python3

import argparse

#########################
# USER ARGUMENT PARSING #
#########################
parser = argparse.ArgumentParser(description="Read arguments")
parser.add_argument("-i","--input",type=str,help="Input exon fasta file")
parser.add_argument("-l","--length",type=str,help="Length of read (this determines the size of the inferred boundary between two exons)")
args=parser.parse_args()
# input file parsing
exonfile = args.input
if exonfile is None:
	print("ERROR: no exon file (-i) specified")
	sys.exit(0)
else:
	print("Exon file = "+exonfile)
# boundary length parsing
ReadLength = args.length
if ReadLength is None:
	print("ERROR: no read length (-l) specified")
	sys.exit(0)
else:
	if int(ReadLength)<2:
		print("ERROR: read length (-l) must be longer than 1")
		sys.exit(0)
	else:
		ReadLength=int(ReadLength)
		print("Read length = "+str(ReadLength))

# function to return widest-spanning boundary (according to read length) between two exons
def BoundaryFinder(exon1,exon2,ReadLength):
	spliced = exon1+exon2
	# define start & end positions
	# NB numbers here are 1 smaller than might expect because python counts from 0
	start = len(exon1)-(ReadLength-1)
	end = len(exon1)+(ReadLength-1)
	boundary = spliced[start:end]
	return(boundary)

def GenomeWideBoundaryFinder(infile):
	# read in gene names, exon names and exon seqs
	GeneNames = []
	ExonNames = []
	ExonSeqs = []
	tempseq = ''
	for line in open(exonfile,'r'):
		if line.startswith('>'):
			ExonNames.append(line.strip('>').strip('\n'))
			name = line.split(':')[0].strip('>')
			if name not in GeneNames:
				GeneNames.append(name)
			# if tempseq has anything in it, add it to ExonSeqs (NB: this will always be one behind the exon name that has been added)
			if tempseq != '':
				ExonSeqs.append(tempseq)
				tempseq = ''
		else:
			# add exon sequence to tempseq (not directly appending to list because of wrapping)
			tempseq += line.strip('\n')
	# add the last exon sequence
	ExonSeqs.append(tempseq)
	print(str(len(ExonNames)) + ' exon names and ' + str(len(ExonSeqs)) + ' sequences read, coming from ' + str(len(GeneNames)) + ' genes')

	# group exon sequences by gene
	GenesExons = {}
	for gene in GeneNames:
		exons = []
		for i in range(len(ExonNames)):
			if gene in ExonNames[i]:
				exons.append(ExonSeqs[i])
		GenesExons[gene] = exons

	# for each gene, generate all exon boundaries
	GenesBoundaries = {}
	for gene in GenesExons:
		boundaries = []
		# for all exons apart from the very last exon... (the first exon of the pair)
		for i in range(0,len(GenesExons[gene])):
			# for all exons apart from the very first exon... (the second exon of the pair)
			for j in range(1,len(GenesExons[gene])):
				# prevent an exon being used with itself (e.g. exon1+exon1), and exons being incorrectly ordered (e.g. exon2+exon1)
				if j > i:
					# define exons to the left and right of junction, and find all boundaries between them (reading into tempbound first to prevent list-of-lists accumulating)
					leftexon = GenesExons[gene][i]
					rightexon = GenesExons[gene][j]
					boundaries.append(BoundaryFinder(leftexon,rightexon,ReadLength))
		# if there are boundaries, add the boundaries to the dict of boundaries for each gene
		if boundaries != []:
			# check that number of boundaries found matches expectation
			expected = 0
			for i in range(1,len(GenesExons[gene])):
				expected += i
			if expected == len(boundaries):
				# if expectation matched, add gene-boundaries pair to dict
				GenesBoundaries[gene]=boundaries
			else:
				# if expectation not matched, throw error
				print('ERROR: ' + gene + ' - ' + str(expected) + ' boundaries expected, ' + str(len(boundaries)) + ' reported')
				raise SystemExit
		print(str(len(boundaries)) + ' boundaries found for ' + gene + ' exons')
	print('All exon boundaries found for all genes')
	return(GenesBoundaries)

GenomeWideBoundaryFinder(infile=exonfile)
