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

# user inputs go here
InputFasta = 'test.fas'
if InputFasta.endswith('.fasta'):
	InputFasta.replace('.fasta','.fas')
Genome = 'dsim-all-chromosome-r2.02.fasta'


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
				raise SystemExit
			else:
				ExonNames.append(line.strip('>').strip('\n'))			
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
		raise SystemExit
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

def ExonerateParser(querydir):
	#######################################
	## read query exon names and lengths ##
	#######################################
	QueryExonLengths = {}
	for file in os.listdir(querydir):
		for line in open(file,"r"):
			if line.startswith('>'):
				# read in exon name (first part of sequence title line)
				if re.search('\:',line) == None:
					print('ERROR - EXON NAME ' + line + ' IS NOT IN THE FORM ">GENE:EXON"')
					raise SystemExit
				else:
					tempname=line.strip('>').strip('\n')		
			else:
				# read in seq length
				templen = len(line.strip('\n'))
		QueryExonLengths(tempname) = templen
	##########################################
	## read target exon names and sequences ##
	##########################################
	TargetExonNames = []
	TargetExonLengths = []
	for line in open(outfile,'r'):
		if line.startswith('sugar:'):
			# read in exon name (keeping same as Query name to allow 1:1 matching)
			TargetExonNames.append(line.split(' ')[1])
			# read in exon match length
			TargetExonLengths.append()
	#############################################################################
	## find exons & genes with length differences between old and new versions ##
	#############################################################################
	# find exons with different lengths
	ExonLengthDiffs = {}
	for i in range(len(NewExonNames)):
		NewExonLen = len(NewExonSeqs[i])
		for j in range(len(ExonNames)):
			if NewExonNames[i] == ExonNames[j]:
				OldExonLen = len(ExonSeqs[j])
				if NewExonLen != OldExonLen:
					ExonLengthDiffs[NewExonNames[i]] = abs(OldExonLen-NewExonLen)
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
	DiffExonsProp = (len(ExonLengthDiffs)/len(ExonNames))*100
	AllGenes = []
	for exon in ExonNames:
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
	outname = outfile.replace('.fas','_GeneLengthDiffs.pdf')
	outpdf = PdfPages(outname)
	plt.savefig(outpdf, format='pdf')
	outpdf.close()
	# output histogram of length differences per gene
	ExonLengthDiffList = []
	for exon in ExonLengthDiffs:
		ExonLengthDiffList.append(ExonLengthDiffs[exon])
	n, bins, patches = plt.hist(np.asarray(ExonLengthDiffList), 50, facecolor='blue', alpha=0.75)
	plt.xlabel('Length Difference')
	plt.ylabel('Count')
	plt.title('Length differences per exon')
	outname = outfile.replace('.fas','_ExonLengthDiffs.pdf')
	outpdf = PdfPages(outname)
	plt.savefig(outpdf, format='pdf')
	outpdf.close()

# function to find conserved exons in one genome based on exons in another genome
def ConservedExonFinder(fastafile):
	ExonDir = FastaSplitter(fastafile=InputFasta)
	ExonerateOutput = ExonerateCaller(querydir=ExonDir,targetgenome=Genome,outfile=InputFasta.replace('.fas','.exonerate'))

ConservedExonFinder(fastafile=InputFasta)
