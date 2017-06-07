# ASE
## Purpose
A set of tools to analyse allele-specific expression (ASE).
## Requirements
Written in python3.

Requires:

[Exonerate](http://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate) v2.2.0 or later

[bedtools](http://bedtools.readthedocs.io/en/latest/) v2.25.0 or later

### ConstitutiveExonFinder
Takes a fasta file of exons, returns a fasta file of exons that are present in every splice form of the parent gene.
**Note: input MUST be in FlyBase format**
e.g. `>GENE:EXON type=TYPE; loc=LOCATION; parent=TRANSCRIPTS; MD5=ID; release=RELEASE; species=SPECIES; length=LENGTH`
Example usage:
```bash
ConstitutiveExonFinder.py -e exons.fasta
```
