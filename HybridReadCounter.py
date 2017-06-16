#!/usr/bin/env python3

import os
import subprocess
import argparse
import sys
import re

#########################
# USER ARGUMENT PARSING #
#########################
parser = argparse.ArgumentParser(description="Read arguments")
parser.add_argument("-r","--reads",type=str,help="Reads fastq file (gzipped)")
parser.add_argument("-g1","--genome1",type=str,help="Bowtie1 index for genome of species 1")
parser.add_argument("-e1","--exons1",type=str,help="Bowtie1 index for exons of species 1")
parser.add_argument("-n1","--name1",type=str,help="Name of species 1")
parser.add_argument("-g2","--genome2",type=str,help="Bowtie1 index for genome of species 2")
parser.add_argument("-e2","--exons2",type=str,help="Bowtie1 index for exons of species 2")
parser.add_argument("-n2","--name2",type=str,help="Name of species 2")
parser.add_argument("-c","--cores",type=int,help="Number of cores")
args=parser.parse_args()
# reads file parsing
Reads = args.reads
if Reads is None:
	print("ERROR: no reads file (-r) specified")
	sys.exit(0)
else:
	print("Reads file = "+Reads)
# genome 1 file parsing
Genome1 = args.genome1
if Genome1 is None:
	print("ERROR: no genome 1 file (-g1) specified")
	sys.exit(0)
else:
	print("Genome 1 file = "+Genome1)
# exons 1 file parsing
Exons1 = args.exons1
if Exons1 is None:
	print("ERROR: no exons 1 file (-e1) specified")
	sys.exit(0)
else:
	print("Exons 1 file = "+Exons1)
# name 1 parsing
Name1 = args.name1
if Name1 is None:
	print("ERROR: no name 1 file (-n1) specified")
	sys.exit(0)
else:
	print("Name 1 = "+Name1)
# genome 2 file parsing
Genome2 = args.genome2
if Genome2 is None:
	print("ERROR: no genome 2 file (-g2) specified")
	sys.exit(0)
else:
	print("Genome 2 file = "+Genome2)
# exons 2 file parsing
Exons2 = args.exons2
if Exons2 is None:
	print("ERROR: no exons 2 file (-e2) specified")
	sys.exit(0)
else:
	print("Exons 2 file = "+Exons2)
# name 2 parsing
Name2 = args.name2
if Name2 is None:
	print("ERROR: no name 2 file (-n2) specified")
	sys.exit(0)
else:
	print("Name 2 = "+Name2)
# cores parsing
Cores = args.cores
if Cores is None:
	print("Using default number of cores (1)")
	Cores = 1
else:
	if Cores>1:
		print("Using "+str(Cores)+" cores")
	else:
		print("ERROR: cores (-c) must be 1 or more")
		sys.exit(0)

def read_ID_bowtie1(file,genome1,genome2,name1,name2,cores):
	####################################################
	## Map read file to Species1 and Species2 genomes ##
	####################################################
	# unzip read file
	cmd = 'gunzip ' + file
	subprocess.call(cmd,shell=True)
	print(file + ' readfile unzipped')
	# count total number of reads (to be used as scale for summary stats below)
	linecount = 0
	for i in open(file.replace('.fq.gz','.fq')):
		linecount+=1
	TotalReads = linecount/4
	# Species1: map reads to genome
	OutReads1 = file.replace('.fq.gz','_'+name1+'.fq')
	OutAlign1 = file.replace('.fq.gz','_'+name1+'.sam')
	cmd = 'bowtie --sam -p '+ str(Cores) +' -v 0 -y --al '+ OutReads1 + ' ' + genome1 + ' ' + file.replace('.fq.gz','.fq') + ' > ' + OutAlign1
	subprocess.call(cmd,shell=True)
	print(file + ' mapped to ' + name1)
	# Dsim: map reads to genome (NB -m 1 = only keep uniquely-mapping reads)
	OutReads2 = file.replace('.fq.gz','_'+name2+'.fq')
	OutAlign2 = file.replace('.fq.gz','_'+name2+'.sam')
	cmd = 'bowtie --sam -p '+ str(Cores) +' -m 1 -v 0 -y --al '+ OutReads2 + ' ' + genome2 + ' ' + file.replace('.fq.gz','.fq') + ' > ' + OutAlign2
	subprocess.call(cmd,shell=True)
	print(file + ' mapped to '+name2)
	# zip read file
	cmd = 'gzip ' + file.replace('.fq.gz','.fq')
	subprocess.call(cmd,shell=True)
	print(file.replace('.fq.gz','.fq') + ' readfile zipped')
	########################################
	## Find reads mapping to both genomes ##
	########################################
	Species1Reads = []
	Species2Reads = []
	for line in open(OutReads1,'r'):
		if line.startswith('@'):
			Species1Reads.append(line.strip('\n'))
	for line in open(OutReads2,'r'):
		if line.startswith('@'):
			Species2Reads.append(line.strip('\n'))
	print(name1+' total reads = ' + str(len(Species1Reads)))
	print(name2+' total reads = ' + str(len(Species2Reads)))	
	# find shared reads using set().intersection() method (this is only used for validation of numbers of Species1-only & Species2-only reads)
	if len(Species1Reads) < len(Species2Reads):
		SharedReads = list(set(Species2Reads).intersection(set(Species1Reads)))
	else:
		SharedReads = list(set(Species1Reads).intersection(set(Species2Reads)))
	################################################
	## Find Species1-only and Species2-only reads ##
	################################################
	# find reads in Species1Reads that are not in Species2Reads using set comparison
	Species1ReadsOnly = list(set(Species1Reads) - set(Species2Reads))
	print(name1+'-specific reads = ' + str(len(Species1ReadsOnly)) + ' (should be ' + str(len(Species1Reads)-len(SharedReads)) + ')')
	# use set difference method to find reads in Species2Reads that are not in SharedReads list
	Species2ReadsOnly = list(set(Species2Reads) - set(Species1Reads))
	print(name2+'-specific reads = ' + str(len(Species2ReadsOnly)) + ' (should be ' + str(len(Species2Reads)-len(SharedReads)) + ')')
	#######################
	## Output read stats ##
	#######################
	# calculate proportions of shared, mel-specific and sim-specific reads in total library
	SharedReadsProp = (len(SharedReads)/TotalReads)*100
	Species1ReadsOnlyProp = (len(Species1ReadsOnly)/TotalReads)*100
	Species2ReadsOnlyProp = (len(Species2ReadsOnly)/TotalReads)*100
	# construct output line
	output = file.replace('.fq.gz','') + '\t' + str(int(TotalReads)) + '\t' + str(round(SharedReadsProp,2)) + '\t\t' + str(round(Species1ReadsOnlyProp,2)) + '\t\t\t' + str(round(Species2ReadsOnlyProp,2)) + '\n'
	# if the Mapping Summary file exists already, write only the output line
	if os.path.exists('MappingSummary.txt') is True:
		MappingSummary = open('MappingSummary.txt','a')
		MappingSummary.write(output)
		MappingSummary.close()
	# if the Mapping Summary file doesn't exist already, create it and write a header line followed by the output line
	elif os.path.exists('MappingSummary.txt') is False:
		header = 'Library\t\tTotal Reads\tShared (%)\t'+name1+'-specific (%)\t'+name2+'-specific (%)\n'
		MappingSummary = open('MappingSummary.txt','a')
		MappingSummary.write(header)
		MappingSummary.write(output)
		MappingSummary.close()
	###############################################################################################
	# Create fastq files for Species1 (excluding Species2) & Species2 (excluding Species1) reads ##
	###############################################################################################
	# output each readname to a file
	Species1ReadsOnlyNames = open(file.replace('.fq.gz','_'+name1+'OnlyNames.txt'),'a')
	for readname in Species1ReadsOnly:
		Species1ReadsOnlyNames.write(readname + '\n')
	Species2ReadsOnlyNames = open(file.replace('.fq.gz','_'+name2+'OnlyNames.txt'),'a')
	for readname in Species2ReadsOnly:
		Species2ReadsOnlyNames.write(readname + '\n')
	SharedReadsOnlyNames = open(file.replace('.fq.gz','_SharedNames.txt'),'a')
	for readname in SharedReads:
		SharedReadsOnlyNames.write(readname + '\n')
	# grep each Species1 readname (from file) and the 3 trailing lines to get read name, sequence, read name and quality score
	cmd = 'grep --no-group-separator -A 3 -F -x -f ' + file.replace('.fq.gz','_'+name1+'OnlyNames.txt') + ' ' + file.replace('.fq.gz','_'+name1+'.fq') + ' > ' + file.replace('.fq.gz','_'+name1+'Only.fq')
	subprocess.call(cmd,shell=True)
	# grep each Species2 readname (from file) and the 3 trailing lines to get read name, sequence, read name and quality score
	cmd = 'grep --no-group-separator -A 3 -F -x -f ' + file.replace('.fq.gz','_'+name2+'OnlyNames.txt') + ' ' + file.replace('.fq.gz','_'+name2+'.fq') + ' > ' + file.replace('.fq.gz','_'+name2+'Only.fq')
	subprocess.call(cmd,shell=True)
	# grep each shared readname (from file) and the 3 trailing lines to get read name, sequence, read name and quality score
	cmd = 'grep --no-group-separator -A 3 -F -x -f ' + file.replace('.fq.gz','_SharedNames.txt') + ' ' + file.replace('.fq.gz','_'+name1+'.fq') + ' > ' + file.replace('.fq.gz','_Shared.fq')
	subprocess.call(cmd,shell=True)
	#####################################################
	## Remove temporary files and compress final files ##
	#####################################################
	os.remove(OutReads1)
	os.remove(OutAlign1)
	os.remove(file.replace('.fq.gz','_'+name1+'OnlyNames.txt'))
	cmd = 'gzip -f ' + file.replace('.fq.gz','_'+name1+'Only.fq')
	subprocess.call(cmd,shell=True)
	os.remove(OutReads2)
	os.remove(OutAlign2)
	os.remove(file.replace('.fq.gz','_'+name2+'OnlyNames.txt'))
	cmd = 'gzip -f ' + file.replace('.fq.gz','_'+name2+'Only.fq')
	subprocess.call(cmd,shell=True)
	os.remove(file.replace('.fq.gz','_SharedNames.txt'))
	cmd = 'gzip -f ' + file.replace('.fq.gz','_Shared.fq')
	subprocess.call(cmd,shell=True)
	# return the names of the species1, species2 & shared read files as a tuple
	return(file.replace('.fq.gz','_'+name1+'Only.fq.gz'),file.replace('.fq.gz','_'+name2+'Only.fq.gz'),file.replace('.fq.gz','_Shared.fq.gz'))

def GFFmaker(fastafile):
	outfile = open(fastafile.replace('.fas','.gff'),'wt')
	source = 'unknown'
	feature = 'exon'
	seqname = ''
	seqstart = 1
	seq = ''
	for l in open(fastafile,'r'):
		if l.startswith('>'):
			# if a sequence has already been read in (i.e. every sequence except the first
			if seq !='':
				seqlength = len(seq)
				seqend = seqstart + seqlength - 1
				entry = seqname + '\t' + source +'\t'+ feature + '\t'+ str(seqstart) +'\t'+ str(seqend) +'\t.\t+\t.\t' + 'ID=' + seqname + '\n'
				outfile.write(entry)
				seqstart = 1
				seq = ''
				seqname = l.strip('\n').strip('>')
			# this deals with the firstname
			else:
				seqname = l.strip('\n').strip('>')
		# add the bases to the seq
		else:
			seq += l.strip('\n')

	seqlength = len(seq)
	seqend = seqstart + seqlength - 1
	entry = seqname + '\t' + source +'\t'+ feature + '\t' + str(seqstart) +'\t'+ str(seqend) +'\t.\t+\t.\t' + 'ID=' + seqname + '\n'
	outfile.write(entry)
	outfile.close()
	return(fastafile.replace('.fas','.gff'))

# function to map species-specific and shared reads to each set of exons, returns 4 bam files
def HybridMapper(species1reads,species2reads,sharedreads,exons1,exons2):
	# uncompress species1reads and species2reads if necessary
	if species1reads.endswith(".gz"):
		cmd="gunzip -f "+species1reads
		subprocess.call(cmd,shell=True)
		species1reads=species1reads.strip(".gz")
	if species2reads.endswith(".gz"):
		cmd="gunzip -f "+species2reads
		subprocess.call(cmd,shell=True)
		species2reads=species2reads.strip(".gz")
	if sharedreads.endswith(".gz"):
		cmd="gunzip -f "+sharedreads
		subprocess.call(cmd,shell=True)
		sharedreads=sharedreads.strip(".gz")
	# map species1reads to exons1
	if re.search("\/",exons1):
		Species1Exons1Out=species1reads.strip(".fq")+"_"+exons1.split("/")[-1]+".sam"
	else:
		Species1Exons1Out=species1reads.strip(".fq")+"_"+exons1+".sam"
	cmd = 'bowtie --sam -p '+ str(Cores) +' -v 0 -y ' + exons1 + ' ' + species1reads + ' > ' + Species1Exons1Out
	subprocess.call(cmd,shell=True)
	cmd = 'samtools view -bS ' + Species1Exons1Out + ' > ' + Species1Exons1Out.replace(".sam",".bam")
	subprocess.call(cmd,shell=True)
	os.remove(Species1Exons1Out)
	# map shared reads to exons1
	if re.search("\/",exons1):
		SharedExons1Out=sharedreads.strip(".fq")+"_"+exons1.split("/")[-1]+".sam"
	else:
		SharedExons1Out=sharedreads.strip(".fq")+"_"+exons1+".sam"
	cmd = 'bowtie --sam -p '+ str(Cores) +' -v 0 -y ' + exons1 + ' ' + sharedreads + ' > ' + SharedExons1Out
	subprocess.call(cmd,shell=True)
	cmd = 'samtools view -bS ' + SharedExons1Out + ' > ' + SharedExons1Out.replace(".sam",".bam")
	subprocess.call(cmd,shell=True)
	os.remove(SharedExons1Out)
	# map species2reads to exons2
	if re.search("\/",exons2):
		Species2Exons2Out=species2reads.strip(".fq")+"_"+exons2.split("/")[-1]+".sam"
	else:
		Species2Exons2Out=species2reads.strip(".fq")+"_"+exons2+".sam"
	cmd = 'bowtie --sam -p '+ str(Cores) +' -v 0 -y ' + exons2 + ' ' + species2reads + ' > ' + Species2Exons2Out
	subprocess.call(cmd,shell=True)
	cmd = 'samtools view -bS ' + Species2Exons2Out + ' > ' + Species2Exons2Out.replace(".sam",".bam")
	subprocess.call(cmd,shell=True)
	os.remove(Species2Exons2Out)
	# map shared reads to exons2
	if re.search("\/",exons2):
		SharedExons2Out=sharedreads.strip(".fq")+"_"+exons2.split("/")[-1]+".sam"
	else:
		SharedExons2Out=sharedreads.strip(".fq")+"_"+exons2+".sam"
	cmd = 'bowtie --sam -p '+ str(Cores) +' -v 0 -y ' + exons2 + ' ' + sharedreads + ' > ' + SharedExons2Out
	subprocess.call(cmd,shell=True)
	cmd = 'samtools view -bS ' + SharedExons2Out + ' > ' + SharedExons2Out.replace(".sam",".bam")
	subprocess.call(cmd,shell=True)
	os.remove(SharedExons2Out)
	return(Species1Exons1Out.replace(".sam",".bam"),SharedExons1Out.replace(".sam",".bam"),Species2Exons2Out.replace(".sam",".bam"),SharedExons2Out.replace(".sam",".bam"))

# function to generate a count file from a bam file and a gff file
def ReadCounter(bamfile,gff):
	print("GFF file = "+gff)
	print("Bam file = "+bamfile)
	cmd="bedtools coverage -s -counts -a "+gff+" -b "+bamfile+" > "+bamfile.replace(".bam",".counts")
	subprocess.call(cmd,shell=True)
	print("Counts generated for "+bamfile+" based on "+gff)
	return(bamfile.replace(".bam",".counts"))

# function to parse counts from different read files into one output file
def ExonCountParser(countfile,exonfile):
	OutFile=open(countfile.replace(".counts",".genecounts"),"wt")
	OutFile.write("Gene\tCount\n")
	# gather gene names from the exonfile
	GeneNames=[]
	for line in open(exonfile):
		if line.startswith(">"):
			gene=line.strip(">").split(":")[0]
			if gene not in GeneNames:
				GeneNames.append(gene)
	# gather counts for each gene
	readcounts={}
	for line in open(countfile):
		splitline=line.split("\t")
		gene=splitline[0].split(":")[0]
		count=int(splitline[-1])
		# either add to the count for this gene, or create a new name:count entry
		if gene in readcounts:
			readcounts[gene]+=count
		else:
			readcounts[gene]=count
	# output gene counts to file
	for gene in GeneNames:
		OutFile.write(gene+"\t"+str(readcounts[gene])+"\n")
	OutFile.close()
	return(countfile.replace(".counts",".genecounts"))

# function to parse multiple gene count files, and output a summary file of counts for each gene in different species
def GeneCountParser(filelist,librarynames,outname):
	# delete file named "outname" if this already exists
	if os.path.isfile(outname):
		os.remove(outname)
	# write blank outfile
	OutFile=open(outname,"wt")
	# get the gene names from the first file
	GeneNames=[]
	for line in open(filelist[0],"r"):
		if not line.startswith("Gene"):
			GeneNames.append(line.split("\t")[0])
	# create dict to hold libraryname:counts pairs for each library
	LibraryCounts={}
	# get name and counts for each file
	for i in range(len(filelist)):
		name=librarynames[i]
		genes=[]
		counts=[]
		for line in open(filelist[i]):
			if not line.startswith("Gene"):
				genes.append(line.split("\t")[0])
				counts.append(line.split("\t")[1].strip("\n"))
		# check that genes are in the same order as GeneNames list before adding to LibraryCounts dict
		for j in range(len(GeneNames)):
			if GeneNames[j]!=genes[j]:
				print("ERROR: gene names do not match (are they in the right order?)")
				sys.exit(0)
		# add name and counts as a pair into the LibraryCounts dict
		LibraryCounts[name]=counts
	# write header line to content string
	OutContent="Gene"
	for library in librarynames:
		OutContent+="\t"+library
	OutContent+="\n"
	# for each gene, add a new line with the gene name as the first entry, and the count for each library (in the same order) as subsequent entries
	for i in range(len(GeneNames)):
		line=GeneNames[i]
		for library in librarynames:
			line+="\t"+LibraryCounts[library][i]
		OutContent+=line+"\n"
	# write OutContent to file
	OutFile.write(OutContent)
	OutFile.close()
	return(outname)
	

# generate species-specific and shared read file
SeparatedReads=read_ID_bowtie1(file=Reads,genome1=Genome1,genome2=Genome2,name1=Name1,name2=Name2,cores=Cores)
print("The species 1 reads can be found at "+SeparatedReads[0])
# make bowtie1 databases for exon1 and exon2 files
cmd="bowtie-build "+Exons1+" "+Exons1.strip(".fas")
subprocess.call(cmd,shell=True)
cmd="bowtie-build "+Exons2+" "+Exons2.strip(".fas")
subprocess.call(cmd,shell=True)
# make gff files for exon1 and exon2 files
Exons1GFF=GFFmaker(fastafile=Exons1)
Exons2GFF=GFFmaker(fastafile=Exons2)
# map reads to exons
MappedReads=HybridMapper(species1reads=SeparatedReads[0],species2reads=SeparatedReads[1],sharedreads=SeparatedReads[2],exons1=Exons1.strip(".fas"),exons2=Exons2.strip(".fas"))
# count reads mapping to each exon
Species1Exons1Counts=ReadCounter(bamfile=MappedReads[0],gff=Exons1GFF)
SharedExons1Counts=ReadCounter(bamfile=MappedReads[1],gff=Exons1GFF)
Species2Exons2Counts=ReadCounter(bamfile=MappedReads[2],gff=Exons2GFF)
SharedExons2Counts=ReadCounter(bamfile=MappedReads[3],gff=Exons2GFF)
# parse individual exon counts into gene counts
Species1Genes1Counts=ExonCountParser(countfile=Species1Exons1Counts,exonfile=Exons1)
SharedGenes1Counts=ExonCountParser(countfile=SharedExons1Counts,exonfile=Exons1)
Species2Genes2Counts=ExonCountParser(countfile=Species2Exons2Counts,exonfile=Exons2)
SharedGenes2Counts=ExonCountParser(countfile=SharedExons2Counts,exonfile=Exons2)
# parse individual gene count files into one summary file
GeneCountParser(filelist=[Species1Genes1Counts,SharedGenes1Counts],librarynames=[Species1Genes1Counts.split("_")[0]+"_"+Species1Genes1Counts.split("_")[1]+"_"+Name1,SharedGenes1Counts.split("_")[0]+"_"+SharedGenes1Counts.split("_")[1]+"_Shared"],outname=Name1+"Genes.counts")
GeneCountParser(filelist=[Species2Genes2Counts,SharedGenes2Counts],librarynames=[Species2Genes2Counts.split("_")[0]+"_"+Species2Genes2Counts.split("_")[1]+"_"+Name2,SharedGenes2Counts.split("_")[0]+"_"+SharedGenes2Counts.split("_")[1]+"_Shared"],outname=Name2+"Genes.counts")

