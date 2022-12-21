# Parameters

All parameters should be specified in the `pipeline.ini` configuration file

 The template file can be found in `ERvsearch/templates/pipeline.ini`.
 
 Make a copy of this file in your working directory (the directory where you want to run the program and store the output) and change the values of the parameters according to your system and your needs.


## Required parameters
These parameters should always be set (there are no default options).
1. [genome](#genome)<br>
1.1. [file](#file)<br>
2. [paths](#paths)<br>
2.1 [path_to_usearch](#path-to-usearch)<br>
2.2 [path_to_exonerate](#path-to-exonerate)<br>

### genome
*Input file parameters*


#### file
`string`


Path to a single fasta file containing the genome or other sequence you would like to screen for ORFs


e.g. `/home/katy/genomes/hg38.fasta`
<br><br>

### paths
*Paths to software*

#### path_to_usearch
`string`


Path to usearch executable

e.g. `/usr/bin/usearch11.0.667_u86linux32`

#### path_to_exonerate
'string'


Path to exonerate executable

e.g. `/usr/bin/exonerate`
<br><br>

## Optional Parameters

The following parameters are optional (as they have default values)

1. [genomesplits](#genomesplits)<br>
1.1 [split](#split)<br>
1.2 [split_n](#split-n)<br>
1.3 [force](#force)<br>
2. [output](#output)<br>
2.1 [outfile_stem](#outfile-stem)<br>
3. [database](#database)<br>
3.1 [use_custom_db](#use-custom-db)<br>
3.2 [gag](#gag)<br>
3.3 [pol](#pol)<br>
3.4 [env](#env)<br>
4. [exonerate](#exonerate)<br>
4.1 [overlap](#overlap)<br>
4.2 [min_score](#min-score)<br>
5. [usearch](#usearch)<br>
5.1 [min_id](#min-id)<br>
5.2 [min_hit_length](#id1)<br>
5.3 [min_coverage](#min-coverage)<br>
6. [plots](#plots)<br>
6.1 [dpi](#dpi)<br>
6.2 [format](#format)<br>
6.3 [gag_colour](#gag-colour)<br>
6.4 [pol_colour](#pol-colour)<br>
6.5 [env_colour](#env-colour)<br>
7. [ORFs](#orfs)<br>
7.1 [min_orf_len](#min-orf-len)<br>
7.2 [translation_table](#translation-table)<br>
8. [trees](#trees)<br>
8.1 [use_gene_colour](#use-gene-colour)<br>
8.2 [maincolour](#maincolour)<br>
8.3 [highlightcolour](#highlightcolour)<br>
8.4 [outgroupcolour](#outgroupcolour)<br>
8.5 [dpi](#id2)<br>
8.6 [format](#format)<br>
9. [regions](#regions)<br>
9.1 [maxdist](#maxdist)<br>

### genomesplits
*Parameters for dividing the input genome into batches*

If the genome has more than ~100 contigs or scaffolds, it is recommended to batch these rather than running Exonerate on each contig individually, to avoid creating an excessive number of output files

The default is to split the input into 50 batches, this is a manageable number for most systems.

#### split
`string` `True` or `False`<br>
Default: `True`


If split is True the contigs will be batched, if it is False Exonerate will run once for every gene-contig combination (usually 3x number of contigs). If there are less contigs than batches this will be ignored


#### split_n
`integer`<br>
Default: `50`


Number of batches to split the contigs into

#### force
`string` `True` or `False`<br>
Default: `False`


The pipeline will error if running using these genome split settings will run Exonerate more than 500 times. If you want to force the pipeline to run despite this, change force to True.
<br><br>

### output
*Output file parameters*
#### outfile_stem
`string`<br>
Default: `ERVsearch`


Log files will have this as a prefix (e.g. ERVsearch_log.txt)


### database
*Parameters concerning the database of reference ERV sequences*
`string` `True` or `False`<br>
#### use_custom_db
Default: `False`


When screening using Exonerate, query sequences are used to identify ERV like regions of the genome. It is possible to use the default ERV  database provided with the pipeline (recommended) or to use a custom database.

False - use the default database

True - use a custom database

#### gag
`string`<br>
Default: `None`


Path to a custom fasta file of gag amino acid sequences. To skip this gene use None as this value (only if use_custom_db is True)


#### pol
`string`<br>
Default: `None`


Path to a custom fasta file of pol amino acid sequences. To skip this gene use None as this value (only if use_custom_db is True).


#### env
`string`<br>
Default: `None`


Path to a custom fasta file of env amino acid sequences. To skip this gene use None as this value (only if use_custom_db is True). 


### exonerate
*Exonerate parameters*
#### min_hit_length
`integer`<br>
Default: `100`


Minimum Exonerate hit length on the chromosome. Shorter hits are filtered out.

#### overlap
`integer`<br>
Default: `30`


Maximum distance between chromosome regions identified with Exonerate which are merged into single regions.


#### min_score
`integer` <br>
Default: `100`


Minimum score in the second exonerate pass (with ungapped algorithm).

### usearch
*USEARCH parameters*
#### min_id
`float`<br>
Default: `0.5`

Percentage identity used by the UBLAST algorithm. Used to set the `-id` UBLAST parameter


#### min_hit_length
`integer`<br>
Default: `100`


Minimum hit length for UBLAST. Used to set the `-mincols` UBLAST parameter

#### min_coverage
`float` <br>
Default: `0.5`

Minimum proportion of the query sequence which should be covered using UBLAST. Used to set the `-query_cov` UBLAST parameter


### plots
*Plotting parameters*

#### dpi
`integer`<br>
Default: `300`

Dots per inch for output plots (from the summarise functions).


#### format
`string` `png`, `svg`, `pdf` or `jpg`<br>

File format for summary plot files, can be svg, png, pdf or jpg


#### gag_colour
`string` (hex colour code with #) <br>
Default: `#f38918`


Colour for *gag* gene in summary plots. Default is orange.

#### pol_colour
`string` (hex colour code with #) <br>
Default: `#4876f9`


Colour for *pol* gene in summary plots. Default is blue.

#### env_colour
`string` (hex colour code with #) <br>
Default: `#d61f54`


Colour for *env* gene in summary plots. Default is pink


#### other_colour
`string` (hex colour code with #) <br>
Default: `#33b54d`<br>


Colour for anything which doesn't relate to a specific gene in summary plots. Default is green.


#### match_axes
`string` `True` or `False`<br>
Default: `False`


If True, when gag, pol and env are shown as subplots on the same figure they should all have the same axis limits. If they do some can be very small but they are more comparable.

### ORFs
*Parameters for ORF identification*

#### min_orf_len
`integer`<br>
Default: `100`


Minimum length of ORFs to characterise.


#### translation_table
`integer` <br>
Default `1`


Translation table to use when identifying ORFs - listed here https://www.ncbi.nlm.nih.gov/Taxonomy/Utils/wprintgc.cgi. Usually table 1 is fine - default plus alternative initiation codons


### trees
*Phylogenetic tree parameters*

#### use_gene_colour
`string` `True` or `False`<br>
Default: `True`


Use the colours specified in the plots section for each gene to highlight newly identified ERVs. If this is `True`, `highlightcolour` will be ignored and the gene colours will be used instead.If it is `False`, `highlightcolour` is used to highlight newly identified ERVs.

#### maincolour
`string` (hex colour code with #) <br>
Default: `#000000`


Main colour for text in tree images. Default is black.


#### highlightcolour
`string` (hex colour code with #) <br>
Default : `#382bfe`


Text colour for leaves highlighted in tree images - newly identified ERV-like regions. This is ignored if `use_gene_colour` above is `True`. Default is blue.


#### outgroupcolour
`string` (hex colour code with #) <br>
Default: `#0e954b`


Text colour for outgroup in tree images. Default is green.


#### dpi
`integer`<br>
Default: `300`
Dots per inch for phylogenetic tree images.


#### format
`string` `png`, `svg`, `pdf` or `jpg`<br>
Default: `png`


File format for phylogenetic tree images, can be svg, png, pdf or jpg


### regions
*Parameters for defining regions with ORFs from more than one retroviral gene*

#### maxdist
`integer`<br>
Default: 3000


Maximum distance (in nucleotides) between ORFs to be defined as part of the same ERV region.

##### maxoverlap
`float`<br>
Default: 0.5


Maximum proportion of an ORF which can overlap with an ORF from another gene before they are filtered out.