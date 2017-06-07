#!/usr/bin/env python3

# script to extract constitutive exons from exon fasta file in FlyBase format
# NB input MUST be in FlyBase format ">GENE:EXON type=TYPE; loc=LOCATION; parent=TRANSCRIPTS; MD5=ID; release=RELEASE; species=SPECIES; length=LENGTH"

import argparse

#########################
# USER ARGUMENT PARSING #
#########################
parser = argparse.ArgumentParser(description="Read arguments")
parser.add_argument("-e","--exons",type=str,help="Exon fasta file")
args=parser.parse_args()
# exon file parsing
exonfile = args.exons
if exonfile is None:
	print("ERROR: no exon file (-e) specified")
	sys.exit(0)
else:
	print("Exon file = "+exonfile)

# read in exon names and sequences
ExonNames = []
ExonSeqs = []
GeneNames = []
TranscriptSets = []
tempseq = ''
for line in open(exonfile,'r'):
	if line.startswith('>'):
		# read in exon names (first part of sequence title line)
		ExonNames.append(line.split(' ')[0].strip('>'))
		# read in gene names (first part of fourth part of sequence title line)
		GeneNames.append(line.split(' ')[3].split(',')[0].strip('parent='))
		# read in transcript names (every part (apart from first) of fourth part of sequence title line)
		TranscriptNameField = line.split(' ')[3].split(',')
		TranscriptNames = []
		for i in TranscriptNameField:
			if i.startswith('FBtr'):
				TranscriptNames.append(i.strip(';'))
		TranscriptSets.append(TranscriptNames)
		# if tempseq has anything in it, add it to ExonSeqs (NB: this will always be one behind the exon name that has been added)
		if tempseq != '':
			ExonSeqs.append(tempseq)
			tempseq = ''
	else:
		# add exon sequence to tempseq (not directly appending to list because of wrapping)
		tempseq += line.strip('\n')
# add the last exon sequence
ExonSeqs.append(tempseq)
print(str(len(ExonNames)) + ' exon names and ' + str(len(ExonSeqs)) + ' sequences read')

# generate complete set of transcripts for each gene
CompleteTranscriptSets = {}
genename = ''
transcriptlist = []
for i in range(len(GeneNames)):
	# if the genename is empty, set it as the current genename (this deals with the first name) and add the transcriptset for that genename to the transcriptlist
	if genename == '':
		genename = GeneNames[i]
		for transcript in TranscriptSets[i]:
			transcriptlist.append(transcript)
	# if the genename is the same as the existing genename, add any previously unseen transcripts to transcriptlist
	elif genename == GeneNames[i]:
		for transcript in TranscriptSets[i]:
			if transcript not in transcriptlist:
				transcriptlist.append(transcript)
	# if a new genename is encountered...
	elif genename != GeneNames[i]:
		# write a new genename:transcriptlist pair to dict
		CompleteTranscriptSets[genename] = transcriptlist
		# set the genename as the new genename
		genename = GeneNames[i]
		# empty the transcriptlist
		transcriptlist = []
		# add the transcriptset for that genename to the transcriptlist
		for transcript in TranscriptSets[i]:
			transcriptlist.append(transcript)		
# write the last gene:transcriptlist pair to dict
CompleteTranscriptSets[genename] = transcriptlist
print(str(len(CompleteTranscriptSets)) + ' complete transcript sets read')

# go through exons, keeping only those that are constitutive (found in all transcripts of the gene)
ConstitutiveExons = []
for i in range(len(ExonNames)):
	gene = GeneNames[i]
	exontranscripts = TranscriptSets[i]
	if len(exontranscripts) == len(CompleteTranscriptSets[gene]):
		ConstitutiveExons.append(ExonNames[i])
print(str(len(ConstitutiveExons)) + ' constitutive exons found')

# write constitutive exons to fasta file
output = ''
for i in range(len(ConstitutiveExons)):
	exon = ConstitutiveExons[i]
	output += '>' + exon + '\n' + ExonSeqs[ExonNames.index(exon)] + '\n'
outfile = open(exonfile.replace('.fasta','_constitutive.fasta'),'wt')
outfile.write(output)
outfile.close()
print('Constitutive exons written to file')
