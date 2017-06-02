#!/usr/bin/env python3

# wrapper script to run exonerate one group of (constitutive) exons at a time
# exons run individually (but multiple in parallel) to improve speed

import os
import subprocess
import shutil
import re
import numpy as np
import matplotlib
# explicitly set the DISPLAY environment variable here
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.backends.backend_pdf import PdfPages
import argparse
import sys

#########################
# USER ARGUMENT PARSING #
#########################
parser = argparse.ArgumentParser(description="Read arguments")
parser.add_argument("-e","--exons",type=str,help="Exon fasta file")
parser.add_argument("-g","--genome",type=str,help="Genome fasta file")
args=parser.parse_args()
# exon file parsing
InputExons = args.exons
if InputExons is None:
	print("ERROR: no exon file (-e) specified")
	sys.exit(0)
else:
	if InputExons.endswith('.fasta'):
		InputExons=InputExons.replace('.fasta','.fas')
		print("Exon file = "+InputExons)
Genome = args.genome
if Genome is None:
	print("ERROR: no genome file (-g) specified")
	sys.exit(0)
else:
	print("Genome file = "+Genome)

# function to split each sequence in a fasta file into its own individual fasta file
def FastaSplitter(fastafile):
	######################################
	## read in exon names and sequences ##
	######################################
	ExonNames = []
	ExonSeqs = []
	tempseq = ''
	for line in open(fastafile,"r"):
		if line.startswith('>'):
			# read in exon name (first part of sequence title line)
			if re.search('\:',line) == None:
				print('ERROR - EXON NAME ' + line + ' IS NOT IN THE FORM ">GENE:EXON"')
				sys.exit(0)
			else:
				exonname=re.sub("\s.*","",line.strip(">").strip("\n"))
				ExonNames.append(exonname)			
			# if tempseq has anything in it, add it to ExonSeqs (NB: this will always be one behind the exon name that has been added)
			if tempseq != '':
				ExonSeqs.append(tempseq)
				tempseq = ''
		else:
			# add exon sequence to tempseq (removing gaps, and not directly appending to list because of wrapping)
			tempseq += line.strip('\n').replace('-','')
	# add the last exon sequence
	ExonSeqs.append(tempseq)
	# check that each exon has a corresponding sequence, throw error and crash if not
	if len(ExonNames) == len(ExonSeqs):
		print(str(len(ExonNames)) + ' exon names and ' + str(len(ExonSeqs)) + ' exon sequences read')
	else:
		print('ERROR - EACH EXON NAME DOES NOT HAVE A CORRESPONDING SEQUENCE')
		sys.exit(0)
	# check that directory exists for temporary exon output files (if not, make it)
	if os.path.isdir('./tempexons') is False:
		cmd = 'mkdir ./tempexons'
		subprocess.call(cmd,shell=True)
	#############################################
	## output each exon to separate fasta file ##
	#############################################
	for i in range(len(ExonNames)):
		# construct string of exon name and seq in fasta format
		exon = '>' + ExonNames[i] + '\n' + ExonSeqs[i]
		# write this string to a temp file
		filename = './tempexons/' + ExonNames[i].replace(':','_').strip('\n').strip('>') + '_exon.fas'
		tempout = open(filename,'wt')
		tempout.write(exon)
		tempout.close()
	print('Exons written to individual files')
	return('./tempexons/')

def ExonerateCaller(querydir,targetgenome,outfile):
	################################
	## run Exonerate on all exons ##
	################################
	# if outfile already exists, delete it
	if os.path.exists(outfile) is True:
		os.remove(outfile)
	# run Exonerate on 10 exons in parallel
	cmd = 'ls ' + querydir + '*_exon.fas | parallel -j 10 \'exonerate --model est2genome --softmasktarget yes --bestn 1 --minintron 20 --maxintron 20000 --showvulgar no --showalignment no --showsugar yes --query {} --target ' + targetgenome + '>> ' + outfile + '\''
	subprocess.call(cmd,shell=True)
	print('Exonerate finished')
	# remove temporary individual exon fasta files
	shutil.rmtree(querydir)
	return(outfile)

def ExonerateParser(exonfasta,exonerateoutput):
	#######################################
	## read query exon names and lengths ##
	#######################################
	QueryExonLengths = {}
	for line in open(exonfasta,"r"):
		if line.startswith('>'):
			# read in exon name (first part of sequence title line)
			if re.search('\:',line) == None:
				print('ERROR - EXON NAME ' + line + ' IS NOT IN THE FORM ">GENE:EXON"')
				sys.exit(0)
			else:
				tempname=line.strip('>').strip('\n')		
		else:
			# read in seq length
			templen = len(line.strip('\n'))
			QueryExonLengths[tempname] = templen
	##########################################
	## read target exon names and sequences ##
	##########################################
	TargetExonNames = []
	TargetExonLengths = []
	for line in open(exonerateoutput,'r'):
		if line.startswith('sugar:'):	
			# read in exon name (keeping same as Query name to allow 1:1 matching)
			TargetExonNames.append(line.split(' ')[1])
			# read in exon match length
			TargetExonLengths.append(int(line.split(' ')[3])-int(line.split(' ')[2]))
	##############################################################################################
	## find exons & genes with length differences between query (old) and target (new) versions ##
	##############################################################################################
	# find exons with different lengths
	ExonLengthDiffs = {}
	for i in range(len(TargetExonNames)):
		TargetExonLen = TargetExonLengths[i]
		for j in QueryExonLengths:
			if TargetExonNames[i] == j:
				QueryExonLen = QueryExonLength[j]
				if TargetExonLen != QueryExonLen:
					ExonLengthDiffs[TargetExonNames[i]] = abs(QueryExonLen-TargetExonLen)
	# get total length differences across exons for each gene
	GeneLengthDiffs = {}
	for exonname in ExonLengthDiffs:
		genename = exonname.split(':')[0]
		if genename not in GeneLengthDiffs:
			GeneLengthDiffs[genename] = ExonLengthDiffs[exonname]
		else:
			GeneLengthDiffs[genename] += ExonLengthDiffs[exonname]
	# get total number of exons with different lengths for each gene
	GeneExonDiffs = {}
	for exonname in ExonLengthDiffs:
		genename = exonname.split(':')[0]
		if genename not in GeneExonDiffs:
			GeneExonDiffs[genename] = 1
		else:
			GeneExonDiffs[genename] += 1
	# calculate overall % of exons and genes with length differences
	DiffExonsProp = (len(ExonLengthDiffs)/len(QueryExonLengths))*100
	AllGenes = []
	for exon in TargetExonNames:
		genename = exon.split(':')[0]
		if genename not in AllGenes:
			AllGenes.append(genename)
	DiffGenesProp = (len(GeneLengthDiffs)/len(AllGenes))*100
	####################
	## output results ##
	####################
	# output summary to file
	if len(ExonLengthDiffs) > 0:
		print(str(len(ExonLengthDiffs)) + ' exons (' + str(DiffExonsProp) + '%) and ' + str(len(GeneLengthDiffs)) + ' genes ('+ str(DiffGenesProp) + '%) are different in length')
		output = str(len(ExonLengthDiffs)) + ' exons (' + str(DiffExonsProp) + '%) and ' + str(len(GeneLengthDiffs)) + ' genes ('+ str(DiffGenesProp) + '%) are different in length\n\nGene\tExons different\tTotal length difference\n'
		for gene in GeneLengthDiffs:
			output += gene + '\t' + str(GeneExonDiffs[gene]) + '\t' + str(GeneLengthDiffs[gene]) + '\n'
		report = open('report.txt','wt')
		report.write(output)
		report.close()
	else:
		print('No genes are different in length')
	# output histogram of length differences per gene
	GeneLengthDiffList = []
	for gene in GeneLengthDiffs:
		GeneLengthDiffList.append(GeneLengthDiffs[gene])
	n, bins, patches = plt.hist(np.asarray(GeneLengthDiffList), 50, facecolor='blue', alpha=0.75)
	plt.xlabel('Length Difference')
	plt.ylabel('Count')
	plt.title('Length differences per gene')
	outname = exonfasta.replace('.fas','_GeneLengthDiffs.pdf')
	outpdf = PdfPages(outname)
	plt.savefig(outpdf, format='pdf')
	outpdf.close()
	# output histogram of length differences per exon
	ExonLengthDiffList = []
	for exon in ExonLengthDiffs:
		ExonLengthDiffList.append(ExonLengthDiffs[exon])
	n, bins, patches = plt.hist(np.asarray(ExonLengthDiffList), 50, facecolor='blue', alpha=0.75)
	plt.xlabel('Length Difference')
	plt.ylabel('Count')
	plt.title('Length differences per exon')
	outname = exonfasta.replace('.fas','_ExonLengthDiffs.pdf')
	outpdf = PdfPages(outname)
	plt.savefig(outpdf, format='pdf')
	outpdf.close()

# function to find conserved exons in one genome based on exons in another genome
def ConservedExonFinder(fastafile):
	ExonDir = FastaSplitter(fastafile=InputExons)
	ExonerateOutput = ExonerateCaller(querydir=ExonDir,targetgenome=Genome,outfile=InputExons.replace('.fas','.exonerate'))
	ExonerateParser(exonfasta=InputExons,exonerateoutput=ExonerateOutput)

ConservedExonFinder(fastafile=InputExons)
