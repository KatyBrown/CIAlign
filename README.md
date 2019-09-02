# CIAlign
CIAlign is a command line tool which performs various functions to parse and analyse a multiple sequence alignment (MSA).

The tool is designed to be highly customisable, allowing users to specify exactly which functions to run and which settings to use.  It is also transparent, generating a clear log file and diagram showing exactly how the alignment has changed.

This allows the user to:
* Remove sources of noise from their MSA
  * Crop of poorly aligned sequence ends
  * Remove of insertions which are not present in the majority of sequences
  * Remove of sequences below a threshold number of bases or amino acids
  * Remove columns containing only gaps
  * Remove sequences above a threshold level percentage of divergence from the majority

* Generate consensus sequences

* Visualise alignments
  * Generate image files showing the alignment before and after parsing and showing which columns and rows have been removed
  * Draw sequence logos
  * Visualise coverage at each postiion in the alignment

* Analyse alignment statistics
  * Generate a similarity matrix showing the percentage identity between each sequence pair
  
## Requirements
* Python >= 3.6
* matplotlib
* numpy
* scipy

The tool runs in Python 3 

## Usage

## Functions

