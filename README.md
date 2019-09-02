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

## Installation
? conda

The current release of CIAlign can be downloaded directly using [this link](https://github.com/KatyBrown/CIAlign/archive/v0.1.0.tar.gz)

Add the CIAlign directory to your PATH environment variable as described [here](https://gist.github.com/nex3/c395b2f8fd4b02068be37c961301caa7)

## Usage
### Basic Usage
`CIAlign --infile INFILE --outfile_stem STEM --inifile cialign.ini`
#### Parameters
| infile | required | path to input alignment | None
| inifile | required | path to ini file | None
| outfile_stem | optional | prefix for output files, including the path to the output directoryy | CIAlign

infile 

### Functions
Specify which functions to run by adding the following arguments to the command



optional arguments:
  -h, --help            show this help message and exit
  --inifile INIFILE     path to input alignment
  --infile INFILE       path to input alignment
  --outfile_stem OUTFILE_STEM
                        stem for output files (including path)
  --remove_insertions   run the removeInsertions function to remove insertions
  --insertion_min_size INSERTION_MIN_SIZE
                        minimum size insertion to remove
  --insertion_max_size INSERTION_MAX_SIZE
                        maximum size insertion to remove
  --insertion_min_flank INSERTION_MIN_FLANK
                        minimum number of bases on either side of deleted
                        insertions
  --remove_short        run the removeShort function to remove sequences with
                        less than n non-gap positions
  --remove_min_length REMOVE_MIN_LENGTH
                        minimum length sequence to remove
  --make_consensus      run the findConsensus function to make a consensus
                        sequence
  --consensus_type CONSENSUS_TYPE
                        type of consensus sequence to make
  --consensus_keep_gaps
                        keep gaps in consensus at positions where a gap is the
                        consensus
  --consensus_name CONSENSUS_NAME
                        name of consensus sequence
  --plot_coverage       plot the coverage as an interpolated function
  --crop_ends           run the cropEnds function to remove badly aligned ends
  --crop_ends_mingap CROP_ENDS_MINGAP
                        minimum gap size to crop from ends
  --remove_badlyaligned
                        run the removeBadlyAligned function to remove badly
                        aligned sequences
  --remove_badlyaligned_minperc REMOVE_BADLYALIGNED_MINPERC
                        minimum percentage identity to majority to not be
                        removed
  --remove_gaponly      run the removeGapOnly function to remove gap only
                        columns from the alignment
  --make_similarity_matrix_input
                        run the calculateSimilarityMatrix function to make a
                        similarity matrix for the input alignment
  --make_similarity_matrix_output
                        run the calculateSimilarityMatrix function to make a
                        similarity matrix for the output alignment
  --make_simmatrix_dp MAKE_SIMMATRIX_DP
                        n decimal places for the similarity matrix (output
                        file only)
  --make_simmatrix_minoverlap MAKE_SIMMATRIX_MINOVERLAP
                        minimum overlap between two sequences to have non-zero
                        similarity in the similarity matrix
  --make_simmatrix_keepgaps MAKE_SIMMATRIX_KEEPGAPS
                        include positions with gaps in either or both
                        sequences in the similarity matrix calculation
  --plot_input          run the drawMiniAlignment function to plot the input
                        alignment
  --plot_output         run the drawMiniAlignment function to plot the output
                        alignment
  --plot_markup         run the drawMiniAlignment function to plot the changes
                        made to the alignment
  --plot_dpi PLOT_DPI   dpi for plots
  --plot_format PLOT_FORMAT
                        plot format (png or svg)
  --plot_width PLOT_WIDTH
                        width for plots (inches)
  --plot_height PLOT_HEIGHT
                        height for plots (inches)
  --make_sequence_logo  draw a sequence logo
  --sequence_logo_type SEQUENCE_LOGO_TYPE
                        type of sequence logo - bar/text/both
  --sequence_logo_dpi SEQUENCE_LOGO_DPI
                        dpi for sequence logo image
  --sequence_logo_font SEQUENCE_LOGO_FONT
                        font for text sequence logo
  --sequence_logo_nt_per_row SEQUENCE_LOGO_NT_PER_ROW
                        number of nucleotides or aas to show per row in the
                        sequence logo
  --sequence_logo_filetype SEQUENCE_LOGO_FILETYPE
                        output file type for sequence logo - png/jpg/svg
  --list_fonts_only     make a swatch showing available fonts
(base) PS C:\Users\Katherine Brown\CIAlign_conda>
