# ASE
## Purpose
A set of tools to analyse allele-specific expression (ASE).
## Requirements
Written in python3.

Requires:

[Exonerate](http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) v2.2.0 or later

[bedtools](http://bedtools.readthedocs.io/en/latest/) v2.25.0 or later

## ConservedExonFinder
Takes a query exon fasta file, a target genome fasta file, and (optionally) a number of cores. Finds coordinates in the target genome that correspond to each query exon, and screens out those that are different in length. Returns a gff file of exon coordinates for the target genome, and a fasta file of exon sequences (all of which have conserved length in query and target) for the query exons.

**Note: input exon (query) and genome (target) files must end ".fas"**

**Note: input exon (query) file must be named in the format**

`GeneName:ExonNumber`

**Example usage:**
```bash
ConservedExonFinder.py -e queryexons.fas -g targetgenome.fas -c 12
```

## ConstitutiveExonFinder
Takes a fasta file of exons, returns a fasta file of exons that are present in every splice form of the parent gene.

**Note: input file name must end ".fas"**

**Note: input MUST be in FlyBase format**

`>GENE:EXON type=TYPE; loc=LOCATION; parent=TRANSCRIPTS; MD5=ID; release=RELEASE; species=SPECIES; length=LENGTH`

**Example usage:**
```bash
ConstitutiveExonFinder.py -e exons.fas
```

## ExonJunctionMaker
Takes a fasta file of exons and a boundary length, returns a fasta file of exon junctions.

**Note: input file must be named in the format**

`GeneName:ExonNumber`

**Example usage:**
```bash
ExonJunctionMaker.py -i exons.fas -l 100
```

## HybridReadCounter
Takes:

* read file
* set of files for species 1 (path to bowtie2 db for genome, exon fasta and species name)
* set of files for species 2 (path to bowtie2 db for genome, exon fasta and species name)

Returns:

* fastq files for species1-specific reads, species2-specific reads and shared reads
* counts of species1 reads and shared reads mapping to species1 exons and genes
* counts of species2 reads and shared reads mapping to species2 exons and genes

**Note: input fasta files must end ".fas"**

**Example usage:**
```bash
~/Scripts/ASE/ReadID.py -c 12 -r Reads.fastq -g1 ./Bowtie2Index/Species1 -n1 Species1 -e1 Species1_exons.fas -g2 Bowtie2Index/Species2 -n2 Species2 -e2 Species2_exons.fas
```

