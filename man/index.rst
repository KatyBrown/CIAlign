CIAlign Documentation
=======================

CIAlign is a command line tool which performs various functions to clean and analyse a multiple sequence alignment (MSA).

Sign up [here](https://t.co/tTyxFV6LR2) for updates when a new feature is added to CIAlign 

The tool is designed to be highly customisable, allowing users to specify exactly which functions to run and which settings to use. It is also transparent, generating a clear log file and alignment markup showing exactly how the alignment has changed and what has been removed by which function.

This allows the user to:

* Remove sources of noise from their MSA
  * Remove insertions which are not present in the majority of sequences
  * Remove sequences below a threshold number of bases or amino acids
  * Crop poorly aligned sequence ends
  * Remove columns containing only gaps
  * Remove sequences above a threshold level percentage of divergence from the majority
  * Remove either end of an alignment where columns don't meet a minimum identity threshold and coverage level

* Generate consensus sequences

* Visualise alignments
  * Generate image files showing the alignment before and after using CIAlign cleaning functions and showing which columns and rows have been removed
  * Draw sequence logos
  * Visualise coverage and conservation at each position in the alignment
  * Generate position frequency, position probability and position weight matrices based on the alignment and produce output formated to be used as input for the [BLAMM](https://github.com/biointec/blamm) and [MEME](https://meme-suite.org/meme/) motif analysis tools.

* Analyse alignment statistics
  * Generate a similarity matrix showing the percentage identity between each sequence pair

* Make changes to the alignment
	* Extract a section of the alignment
	* Unalign the alignment
	* Replace U with T, or T with U in a nucleotide alignment

## Citation

If you found CIAlign useful, please cite: 

[Tumescheit C, Firth AE, Brown K. 2022. CIAlign: A highly customisable command line tool to clean, interpret and visualise multiple sequence alignments. PeerJ 10:e12983 https://doi.org/10.7717/peerj.12983](https://peerj.com/articles/12983/)

.. toctree::
   :maxdepth: 3

   pages/installation.md
   pages/usage.md
