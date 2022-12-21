# Functions
1.[Pipeline Schematic](#pipeline-schematic)<br>
2.[Screen](#screen)<br>
3.[Classify](#classify)<br>
4.[ERVRegions](#ervregions)<br>

## Pipeline Schematic
[large version](../_images/pipeline.png)
![../images/pipeline.png](../images/pipeline.png) 
NB: Where `GENE` is specified in a file path one file will be created for each gene - *gag*, *pol* and *env* (unless otherwise specified by the user).

## Screen
* Screens the genome for ERV like regions by comparing the genome to a set of known retroviral ORFs using Exonerate.
* Confirms the Exonerate regions using UBLAST
* Finds and confirms ORFs within these regions
* Finds the most similar known retroviral ORF in the database to each of the newly identified ORFs

1. [initiate](#initiate)<br>
2. [genomeToChroms](#genometochroms)<br>
3. [prepDBs](#prepdbs)<br>
4. [runExonerate](#runexonerate)<br>
5. [cleanExonerate](#cleanexonerate)<br>
6. [mergeOverlaps](#mergeoverlaps)<br>
7. [makeFastas](#makefastas)<br>
8. [renameFastas](#renamefastas)<br>
9. [makeUBLASTDb](#makeublastdb)<br>
10. [runUBLASTCheck](#runublastcheck)<br>
11. [classifyWithExonerate](#classifywithexonerate)<br>
12. [getORFs](#getorfs)<br>
13. [checkORFsUBLAST](#checkorfsublast)<br>
14. [assignGroups](#assigngroups)<br>
15. [summariseScreen](#summarisescreen)<br>
16. [Screen](#id1)<br>

### initiate

**Input Files** <br>
[`pipeline.ini`](introduction.html#required-input-files)<br>

**Output Files**<br>
[`init.txt`](outputs.html#initiate)<br>

**Parameters**<br>
[`[genome] file`](parameters.html#file)<br>
[`[paths] path_to_usearch`](parameters.html#path-to-usearch)<br>
[`[paths] path_to_exonerate`](parameters.html#path-to-exonerate)<br>

Initialises the pipeline and checks that the required parameters in the pipeline.ini are set and valid and that the required software is in your $PATH.

Checks that:
* The input genome file exists.
* The correct path to ERVsearch is provided.
* samtools, bedtools, FastTree and mafft are in the $PATH
* The correct paths to usearch and exonerate are provided.
 
`init.txt` is a placeholder to show that this step has been completed.


### genomeToChroms

**Input Files**<br>
[`[genome] file`](introduction.html#required-input-files)<br>
[`keep_chroms.txt`](introduction.html#optional-input-files)<br>

**Output Files**<br>
[`host_chromosomes.dir/*fasta`](outputs.html#genometochroms)<br>

**Parameters**<br>
[`[genome] file`](parameters.html#file)<br>
[`[genomesplits] split`](parameters.html#split)<br>
[`[genomesplits] split_n`](parameters.html#split_n)<br>
[`[genomesplits] force`](parameters.html#force)<br>

Splits the host genome provided by the user into FASTA files of a suitable size to run Exonerate efficiently.

If `genomesplits_split` in the pipeline.ini is False, the genome is split into one fasta file for each sequence - each chromosome, scaffold or contig.

If `genomesplits_split` in the pipeline.ini is True, the genome is split into the number of batches specified by the `genomesplits_splitn` parameter, unless the total number of sequences in the input file is less than this number.

The pipeline will fail if the number of sequences which would result from the genomesplits settings would result in >500 Exonerate runs, however it is possible to force the pipeline to run despite this by setting `genomesplits_force` to True.

If the file keep_chroms.txt exists in the working directory only chromosomes listed in this file will be kept.

An unzipped copy of zipped and gzipped fasta files will be created or a link to the file if it is already unzipped. this will be named `genome.fa` and be in the working directory.

This function generates a series of fasta files which are stored in the host_chromosomes.dir directory.

### prepDBs

**Input Files**<br>
None<br>

**Output Files** <br>
[`gene_databases.dir/GENE.fasta`](outputs.html#prepdbs)<br>

**Parameters**<br>
[`[database] use_custom_db`](parameters.html#use-custom-db)<br>
[`[database] gag`](parameters.html#gag)<br>
[`[database] pol`](parameters.html#pol)<br>
[`[database] env`](parameters.html#env)<br>

 Retrieves the gag, pol and env amino acid sequence database fasta files and puts a copy of each gene_databases.dir directory.
    
 If custom databases are used they are retrieved and named as gag.fasta  pol.fasta, env.fasta so the path doesn't need to be changed every time.

### runExonerate

**Input Files** <br>
[`gene_databases.dir/GENE.fasta`](outputs.html#prepdbs)<br>
[`host_chromosomes.dir/*fasta`](outputs.html#genometochroms)<br>

**Output Files** <br>
[`raw_exonerate_output.dir/GENE_*.tsv`](outputs.html#runexonerate)<br>

**Parameters**<br>
[`[paths] path_to_exonerate`](parameters.html#path-to-exonerate)<br>


Runs the `protein2dna` algorithm in the [Exonerate software package](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual) with the host chromosomes (or other regions) in `host_chromosomes.dir` as target sequences and the FASTA files from prepDBs as the query sequences.

The raw output of Exonerate is stored in the raw_exonerate_output directory, one file is created for each combination of query and target sequences.

This step is carried out with low stringency as results are later filtered using UBLAST and Exonerate.
 
### cleanExonerate

**Input Files** <br>
[`raw_exonerate_output.dir/GENE_*.tsv`](outputs.html#runexonerate)<br>

**Output_Files**<br>
[`clean_exonerate_output.dir/GENE_*_unfiltered.tsv`](outputs.html#cleanexonerate)<br>
[`clean_exonerate_output.dir/GENE_*_filtered.tsv`](outputs.html#cleanexonerate)<br>
[`clean_exonerate_output.dir/GENE_*.bed`](outputs.html#cleanexonerate)<br>

**Parameters**<br>
[`[exonerate] min_hit_length`](parameters.html#min-hit-length)<br>

Filters and cleans up the Exonerate output.
* Converts the raw Exonerate output files into dataframes - GENE_unfiltered.tsv
* Filters out any regions containing introns (as defined by Exonerate)
* Filters out regions less than `exonerate_min_hit_length` on the host sequence (in nucleotides).
* Outputs the filtered regions to GENE_filtered.tsv
* Converts this to bed format and outputs this to GENE.bed

### mergeOverlaps

**Input Files** <br>
[`clean_exonerate_output.dir/GENE_*.bed`](outputs.html#cleanexonerate)<br>

**Output_Files** <br>
[`gene_bed_files.dir/GENE_all.bed`](outputs.html#mergeoverlaps)<br>
[`gene_bed_files.dir/GENE_merged.bed`](outputs.html#mergeoverlaps)<br>

**Parameters**<br>
[`[exonerate] overlap`](parameters.html#overlap)<br>

Merges the output bed files for individual sections of the input genome into a single bed file.

Overlapping regions or very close together regions of the genome detected by Exonerate with similarity to the same retroviral gene are then merged into single regions.  This is performed using [bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html) on the bed files output by cleanExonerate.

If there is a gap of less than `exonerate_overlap` between the regions they will be merged.

### makeFastas

**Input Files**<br>
[`gene_bed_files.dir/GENE_merged.bed`](outputs.html#mergeoverlaps)<br>
[`genome.fa`](introduction.html#required-input-files)<br>

**Output Files** <br>
[`gene_fasta_files.dir/GENE_merged.fasta`](outputs.html#makefastas)<br>

**Parameters**<br>
None<br>

Fasta files are generated containing the sequences of the merged regions of the genome identified using mergeOverlaps.
These are extracted from the host chromosomes using [bedtools getfasta](https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html).

### renameFastas

**Input Files** <br>
[`gene_fasta_files.dir/GENE_merged.fasta`](outputs.html#makefastas)<br>

**Output Files** <br>
[`gene_fasta_files.dir/GENE_merged_renamed.fasta`](outputs.html#renamefastas)<br>

**Parameters**<br>
None<br>

 Renames the sequences in the fasta files of ERV-like regions identified with Exonerate so each record has a numbered unique ID (gag1, gag2 etc). Also removes ":" from sequence names as this causes problems later.
 
### makeUBLASTDb

**Input Files**<br>
[`gene_databases.dir/GENE.fasta`](outputs.html#prepdbs)<br>

**Output Files** <br>
[`UBLAST_db.dir/GENE_db.udb`](outputs.html#makeublastdb)<br>

**Parameters**<br>
[`[paths] path_to_usearch`](parameters.html#path-to-usearch)<br>

[USEARCH](https://www.drive5.com/usearch/) requires an indexed database of query sequences to run. This function generates this database for the three gene amino acid fasta files used to screen the genome.

### runUBLASTCheck

**Input Files**<br>
[`UBLAST_db.dir/GENE_db.udb`](outputs.html#makeublastdb)<br>
[`gene_fasta_files.dir/GENE_merged_renamed.fasta`](outputs.html#renamefastas)<br>

**Output Files**<br>
[`ublast.dir/GENE_UBLAST_alignments.txt`](outputs.html#runublastcheck)<br>
[`ublast.dir/GENE_UBLAST.tsv`](outputs.html#runublastcheck)<br>
[`ublast.dir/GENE_filtered_UBLAST.fasta`](outputs.html#runublastcheck)<br>

**Parameters**<br>
[`[paths] path_to_usearch`](parameters.html#path-to-usearch)<br>
[`[usearch] min_id`](parameters.html#min-id)<br>
[`[usearch] min_hit_length`](parameters.html#min-hit-length)<br>
[`[usearch] min_coverage`](parameters.html#min-coverage)<br>

 ERV regions in the fasta files generated by makeFasta are compared to the ERV amino acid database files for a second time, this time using USEARCH (https://www.drive5.com/usearch/). Using both of these tools reduces the number of false positives.

This allows sequences with low similarity to known ERVs to be filtered out.  Similarity thresholds can be set in the pipeline.ini file (`usearch_min_id,` - minimum identity between query and target -  `usearch_min_hit_length` - minimum length of hit on target sequence -  and `usearch_min_coverage` - minimum proportion of the query sequence the hit should cover).

The raw output of running UBLAST against the target sequences is saved in GENE_UBLAST_alignments.txt (equivalent to the BLAST default output) and GENE_UBLAST.tsv (equivalent to the BLAST -outfmt 6 tabular output) this is already filtered by passing the appropriate parameters to UBLAST. The regions which passed the filtering and are therefore in these output files are then output to a FASTA file GENE_filtered_UBLAST.fasta.

### classifyWithExonerate

**Input Files**<br>
[`ublast.dir/GENE_filtered_UBLAST.fasta`](outputs.html#runublastcheck)<br>
[`ERVsearch/ERV_db/all_ERVs_nt.fasta`](introduction.html#databases)<br>

**Output Files**<br>
[`exonerate_classification.dir/GENE_all_matches_exonerate.tsv`](outputs.html#classifywithexonerate)<br>
[`exonerate_classification.dir/GENE_best_matches_exonerate.tsv`](outputs.html#classifywithexonerate)<br>
[`exonerate_classification.dir/GENE_refiltered_matches_exonerate.fasta`](outputs.html#classifywithexonerate)<br>

**Parameters**<br>
[`[paths] path_to_exonerate`](parameters.html#path-to-exonerate)<br>
[`[exonerate] min_score`](parameters.html#min-score)<br>

Runs the Exonerate ungapped algorithm with each ERV region in the fasta files generated by makeFasta as queries and the all_ERVs_nt.fasta fasta file as a target, to detect which known retrovirus is most similar to each newly identified ERV region. Regions which don't meet a minimum score threshold (`exonerate_min_score`) are filtered out.

all_ERVs_nt.fasta contains nucleic acid sequences for many known endogenous and exogenous retroviruses with known classifications.

First all seqeunces are compared to the database and the raw output is saved as exonerate_classification.dirGENE_all_matches_exonerate.tsv. Results need a score greater than `exonerate_min_score`
against one of the genes of the same type (*gag*, *pol* or *env*) in the database. The highest scoring result which meets these critera for each sequence is then identified and output to exonerate_classification.dir/GENE_best_matches_exonerate.tsv. The sequences which meet these critera are also output to a FASTA file exonerate_classification.dir/GENE_refiltered_exonerate.fasta.


### getORFs

**Input Files**<br>
[`exonerate_classification.dir/GENE_refiltered_matches_exonerate.fasta`](outputs.html#classifywithexonerate)<br>

**Output Files**<br>
[`ORFs.dir/GENE_orfs_raw.fasta`](outputs.html#getorfs)<br>
[`ORFs.dir/GENE_orfs_nt.fasta`](outputs.html#getorfs)<br>
[`ORFs.dir/GENE_orfs_aa.fasta`](outputs.html#getorfs)<br>

**Parameters**<br>
[`[orfs] translation_table`](parameters.html#translation-table)<br>
[`[orfs] min_orf_len`](parameters.html#min-orf-len)<br>

Finds the longest open reading frame in each of the ERV regions in the filtered output table.

This analysis is performed using [EMBOSS revseq](http://emboss.open-bio.org/rel/dev/apps/revseq.html) and [EMBOSS transeq](http://emboss.open-bio.org/rel/dev/apps/transeq.html).

The sequence is translated in all six frames using the user specified translation table. The longest ORF is then identified. ORFs shorter than orfs_min_orf_length are filtered out. 

The positions of the ORFs are also convered so that they can be extracted directly from the input sequence file, rather than using the co-ordinates relative to the original Exonerate regions.

The raw transeq output, the nucleotide sequences of the ORFs and the amino acid sequences of the ORFs are written to the output FASTA files.

### checkORFsUBLAST

**Input Files**<br>
[`ORFs.dir/GENE_orfs_nt.fasta`](outputs.html#getorfs)<br>
[`UBLAST_dbs.dir/GENE_db.udb`](outputs.html#makeublastdb)<br>

**Output Files**<br>
[`ublast_orfs.dir/GENE_UBLAST_alignments.txt`](outputs.html#checkorfsublast)<br>
[`ublast_orfs.dir/GENE_UBLAST.tsv`](outputs.html#checkorfsublast)<br>
[`ublast_orfs.dir/GENE_filtered_UBLAST.fasta`](outputs.html#checkorfsublast)<br>

**Parameters**<br>
[`[paths] path_to_usearch`](parameters.html#path-to-usearch)<br>
[`[usearch] min_id`](parameters.html#min-id)<br>
[`[usearch] min_hit_length`](parameters.html#min-hit-length)<br>
[`[usearch] min_coverage`](parameters.html#min-coverage)<br>


ERV ORFs in the fasta files generated by the ORFs function are compared to the original ERV amino acid files using UBLAST. This allows any remaining sequences with poor similarity to known ERVs to be filtered out.

This allows ORFs with low similarity to known ERVs to be filtered out.   Similarity thresholds can be set in the pipeline.ini file (`usearch_min_id,` - minimum identity between query and target -  `usearch_min_hit_length` - minimum length of hit on target sequence -  and `usearch_min_coverage` - minimum proportion of the query sequence the hit should cover).

The raw output of running UBLAST against the target sequences is saved in GENE_UBLAST_alignments.txt (equivalent to the BLAST default output) and GENE_UBLAST.tsv (equivalent to the BLAST -outfmt 6 tabular output) this is already filtered by passing the appropriate parameters to UBLAST. The regions which passed the filtering and are therefore in these output files are then output to a FASTA file GENE_filtered_UBLAST.fasta.

### assignGroups

**Input Files**<br>
[`ublast_orfs.dir/GENE_UBLAST.tsv`](outputs.html#checkorfsublast)<br>
[`ERVsearch/ERV_db/convert.tsv`](introduction.html#databases)<br>

**Output Files**<br>
[`grouped.dir/GENE_groups.tsv`](outputs.html#assigngroups)<br>

**Parameters**<br>
None

Many of the retroviruses in the input database all_ERVs_nt.fasta have been classified into groups based on sequence similarity, prior knowledge and phylogenetic clustering.  Some sequences don't fall into any well defined group, in these cases they are just assigned to a genus, usually based on prior knowledge. The information about these groups is stored in the provided file ERVsearch/ERV_db/convert.tsv.

Each sequence in the filtered fasta file of newly identified ORFs is assigned to one of these groups based on the sequence identified as the most similar in the classifyWithExonerate step. 

The output table is also  tidied up to include the UBLAST output, chromosome, ORF start and end positions, genus and group.

### summariseScreen

**Input Files**<br>
[`gene_bed_files.dir/GENE_merged.bed`](outputs.html#mergeoverlaps)<br>
[`ublast.dir/GENE_UBLAST.tsv`](outputs.html#runublastcheck)<br>
[`ORFs.dir/GENE_orfs_aa.fasta`](outputs.html#getorfs)<br>
[`ublast_orfs.dir/GENE_UBLAST.tsv`](outputs.html#checkorfsublast)<br>

**Output Files**<br>
FMT can be png, svg, pdf, jpg depending on the `plot_format` parameter<br>

[`summary_tables.dir/exonerate_initial_summary.txt`](outputs.html#summarisescreen)<br>
[`summary_tables.dir/ublast_hits_initial_summary.txt`](outputs.html#summarisescreen)<br>
[`summary_tables.dir/orfs_initial_summary.txt`](outputs.html#summarisescreen)<br>
[`summary_tables.dir/ublast_orfs_initial_summary.txt`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/exonerate_initial_lengths.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/exonerate_initial_scores.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/exonerate_initial_strands.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/exonerate_initial_by_sequence.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/exonerate_initial_counts_per_gene.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_hits_alignment_length.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_hits_perc_similarity.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_hits_by_match.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_hits_per_gene.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/orfs_lengths.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/orfs_strands.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/orfs_by_gene.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_orfs_alignment_length.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_orfs_perc_similarity.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_orfs_bit_score.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_orfs_by_match.FMT`](outputs.html#summarisescreen)<br>
[`summary_plots.dir/ublast_orfs_per_gene.FMT`](outputs.html#summarisescreen)<br>
[`screen_results.dir/results.tsv`](introduction.html#id1)<br>
[`screen_results.dir/by_length.FMT`](introduction.html#id1)<br>
[`screen_results.dir/by_genus.FMT`](introduction.html#id1)<br>
[`screen_results.dir/by_group.FMT`](introduction.html#id1)<br>
[`screen_results.dir/by_gene.FMT`](introduction.html#id1)<br>

**Parameters**<br>
[`[plots] dpi`](parameters.html#dpi)<br>
[`[plots] format`](parameters.html#format)<br>
[`[plots] gag_colour`](parameters.html#gag-colour)<br>
[`[plots] pol_colour`](parameters.html#pol-colour)<br>
[`[plots] env_colour`](parameters.html#env-colour)<br>
[`[plots] other_colour`](parameters.html#other-colour)<br>
[`[plots] match_axes`](parameters.html#match-axes)<br>

Generates a series of summary plots and tables showing the results of running the screening functions.

The major outputs of this function are stored in the screen_results.dir directory. Further details of these files are provided in the [Main Output Files](introduction.html#main-output-files) section.

The other files show the output of the intermediate steps.<br>

**Exonerate Initial**<br>
* `summary_tables.dir/exonerate_initial_summary.txt`<br>
  Summary of the output of the initial Exonerate screening step. Note that these are unfiltered and many will not be true ERVs.<br>
* `summary_tables.dir/ublast_hits_initial_summary.txt`<br>
  Summary of the results of running UBLAST on the initial Exonerate output.<br>
* `summary_tables.dir/orfs_initial_summary.txt`<br>
  Summary of the results of the initial ORF identification.<br>
* `summary_tables.dir/ublast_orfs_initial_summary.txt`<br>
  Summary of the results of running UBLAST on these ORFs.<br>
* `summary_plots.dir/exonerate_initial_lengths.FMT`<br>
  Histogram showing the lengths of the initial Exonerate regions for each gene.
* `summary_plots.dir/exonerate_initial_scores.FMT`<br>
  Histogram showing the Exonerate score of the initial Exonerate regions for each gene.<br>
* `summary_plots.dir/exonerate_initial_strands.FMT`<br>
  Bar chart showing the number of regions identified on each strand in the initial Exonerate screen.<br>
* `summary_plots.dir/exonerate_initial_by_sequence.FMT`<br>
  Histogram showing the number of ERV-like regions identified on each sequence in the reference genome being screened.<br>
* `summary_plots.dir/exonerate_initial_counts_per_gene.FMT`<br>
  Bar chart showing the number of ERV regions identified per gene in the initial Exonerate screen.<br>

**UBLAST**<br>
* `summary_plots.dir/ublast_hits_alignment_length.FMT`<br>
  Histogram showing the lengths of the alignments of the UBLAST filtered Exonerate regions and the most similar reference ORF, based on the UBLAST output.<br>
* `summary_plots.dir/ublast_hits_perc_similarity.FMT`<br>
  Histogram showing the percentage identity between the UBLAST filtered Exonerate regions and the most similar reference ORF, based on the UBLAST output.<br>
* `summary_plots.dir/ublast_hits_perc_similarity.FMT`<br>
  Histogram showing the UBLAST bit score between the UBLAST filtered Exonerate regions and the most similar reference ORF, based on the UBLAST output.<br>
* `summary_plots.dir/ublast_hits_by_match.FMT`<br>
  Bar chart showing the number of UBLAST filtered Exonerate regions most similar to each reference ORF in the ERVsearch/ERV_db database.<br>
* `summary_plots.dir/ublast_hits_per_gene.FMT`<br>
  Bar chart showing the number of UBLAST filtered Exonerate regions identified per gene.<br>

**ORFs**<br>
* `summary_plots.dir/orfs_lengths.FMT`<br>
  Histogram of the lengths of ORFs identified in the ERV regions.
* `summary_plots.dir/orfs_strands.FMT`<br>
  Bar chart of the strand (positive (+) or negative (-) sense) of the ORFs identified in the ERV regions.
* `summary_plots.dir/orfs_by_gene.FMT`<br>
  Bar chart of the number of ORFs identified for each gene.<br>

**UBLAST ORFs**<br>
* `summary_plots.dir/ublast_orfs_alignment_length.FMT`<br>
  Histogram showing the lengths of the alignments of the ERV-like ORFs and the most similar reference ORF, based on the UBLAST output.<br>
* `summary_plots.dir/ublast_orfs_perc_similarity.FMT`<br>
  Histogram showing the percentage identity between the ERV-like ORFs and the most similar reference ORF, based on the UBLAST output.<br>
* `summary_plots.dir/ublast_orfs_bit_score.FMT`<br>
  Histogram showing the UBLAST bit score between the ERV-like ORFs and the most similar reference ORF, based on the UBLAST output.<br>
* `summary_plots.dir/ublast_orfs_by_match.FMT`<br>
  Bar chart showing the number of ERV-like ORFs most similar to each reference ORF in the ERVsearch/ERV_db database.<br>
* `summary_plots.dir/ublast_orfs_per_gene.FMT`<br>
  Bar chart showing the number of ERV-like ORFs identified per gene.

### Screen

**Input Files** <br>
None<br>

**Output Files**<br>
 None<br>

**Parameters** <br>
None<br>

Helper function to run all screening functions (all functions prior to this point).

## Classify

* Classifies the newly identified ORFs into groups based on the most similar known ORF
* Aligns the newly identified ORFs with reference sequences within these groups and builds a phylogenetic tree for each group.
* Finds clusters of newly identified ORFs within these trees
* Incorporates representative sequences from these clusters into a summary tree for each retroviral gene and genus  (based on classification into *gamma*, *beta*, *spuma*, *alpha*, *lenti*, *epsilon* and *delta* retroviruses as defined by the [ICTV](https://talk.ictvonline.org/taxonomy).

1. [makeGroupFastas](#makegroupfastas)<br>
2. [makeGroupTrees](#makegrouptrees)<br>
3. [drawGroupTrees](#drawgrouptrees)<br>
4. [makeSummaryFastas](#makesummaryfastas)<br>
5. [makeSummaryTrees](#makesummarytrees)<br>
6. [drawSummaryTrees](#drawsummarytrees)<br>
7. [summariseClassify](#summariseclassify)<br>
8. [Classify](#id2)


### makeGroupFastas

**Input Files**<br>
[`grouped.dir/GENE_groups.tsv`](outputs.html#assigngroups)<br>
[`ERVsearch/phylogenies/group_phylogenies/*fasta`](introduction.html#databases)<br>
[`ERVsearch/phylogenies/summary_phylogenies/*fasta`](introduction.html#databases)<br>
[`ERVsearch/phylogenies/outgroups.tsv`](introduction.html#databases)<br>

**Output Files**<br>
[`group_fastas.dir/GENE_(.*)_GENUS.fasta`](outputs.html#makegroupfastas)<br>
[`group_fastas.dir/GENE_(.*)_GENUS_A.fasta`](outputs.html#makegroupfastas)<br>

**Parameters**<br>
None

Two sets of reference fasta files are available (files are stored in `ERVsearch/phylogenies/group_phylogenies` and `ERVsearch/phylogenies/summary_phylogenies`)

* group_phylogenies - groups of closely related ERVs for fine classification of sequences
* summary_phylogenies - groups of most distant ERVs for broad classification of sequences

Sequences have been assigned to groups based on the most similar sequence in the provided ERV database, based on the score using the Exonerate ungapped algorithm.
Where the most similar sequence is not part of a a well defined group, it has been assigned to a genus.

Fasta files are generated containing all members of the group from the group_phylogenies file (plus an outgroup) where possible and using representative sequences from the same genus, using the summary_phylogenies file, where only a genus has been assigned, plus all the newly identified ERVs in the group. These files are saved as GENE_(group_name_)GENUS.fasta.

A "~" is added to all new sequence names so they can be searched for easily.

The files are aligned using the MAFFT fftns algorithm https://mafft.cbrc.jp/alignment/software/manual/manual.html to generate the GENE_(group_name_)GENUS_A.fasta aligned output files.


### makeGroupTrees

**Input Files**<br>
[`group_fastas.dir/GENE_(.*_)GENUS_A.fasta`](outputs.html#makegroupfastas)<br>

**Output Files**<br>
[`group_trees.dir/GENE_(.*_)GENUS.tre`](outputs.html#makegrouptrees)<br>

**Parameters**<br>
None<br>

Builds a phylogenetic tree, using the FastTree2 algorithm (http://www.microbesonline.org/fasttree) with the default settings plus the GTR model, for the aligned group FASTA files generated by the makeGroupFastas function.

### drawGroupTrees

**Input Files**<br>
[`group_trees.dir/GENE_(.*_)GENUS.tre`](outputs.html#makegrouptrees)<br>

**Output Files**<br>
[`group_trees.dir/GENE_(.*_)GENUS.FMT`](outputs.html#drawgrouptrees) (png, svg, pdf or jpg)<br>

**Parameters**<br>
[`[plots] gag_colour`](parameters.html#gag-colour)<br>
[`[plots] pol_colour`](parameters.html#pol-colour)<br>
[`[plots] env_colour`](parameters.html#env-colour)<br>
[`[trees] use_gene_colour`](parameters.html#use-gene-colour)<br>
[`[trees] maincolour`](parameters.html#maincolour)<br>
[`[trees] highlightcolour`](parameters.html#highlightcolour)<br>
[`[trees] outgroupcolour`](parameters.html#outgroupcolour)<br>
[`[trees] dpi`](parameters.html#dpi)<br>
[`[trees] format`](parameters.html#format)<br>

Generates an image file for each file generated in the makeGroupTrees step, using ete3 (http://etetoolkit.org). Newly identified sequences are labelled as "~" and shown in a different colour.

By default, newly identified sequences are shown in the colours specified in `plots_gag_colour`, `plots_pol_colour` and `plots_env_colour` - to do this then `trees_use_gene_colour` should be set to True in the `pipeline.ini`. Alternatively, a fixed colour can be used by setting `trees_use_gene_colour` to False and settings `trees_highlightcolour`. The text colour of the reference sequences (default black) can be set using `trees_maincolour` and the outgroup using `trees_outgroupcolour`.

The output file DPI can be specified using `trees_dpi` and the format (which can be png, svg, pdf or jpg) using `trees_format`.

### makeSummaryFastas

**Input Files**<br>
[`group_fastas.dir/GENE_(.*_)GENUS.fasta`](outputs.html#makegroupfastas)<br>
[`group_trees.dir/GENE_(*_)GENUS.tre`](outputs.html#makegrouptrees)<br>
[`ERVsearch/phylogenies/summary_phylogenies/GENE_GENUS.fasta`](introduction.html#databases)<br>
[`ERVsearch/phylogenies/group_phylogenies/(.*)_GENUS_GENE.fasta`](introduction.html#databases)<br>

**Output Files**<br>
[`summary_fastas.dir/GENE_GENUS.fasta`](outputs.html#makesummaryfastas)<br>
[`summary_fastas.dir/GENE_GENUS.tre`](outputs.html#makesummaryfastas)<br>

**Parameters**<br>
None

Based on the group phylogenetic trees generated in makeGroupTrees, monophyletic groups of newly idenified ERVs are identified. For each of these groups, a single sequence (the longest) is selected as representative. The representative sequences are combined with the FASTA files in `ERVsearch/phylogenies/summary_phylogenies`, which contain representative sequences for each retroviral gene and genus. These are extended to include further reference sequences from the same small group as the newly identified sequences.

For example, if one MLV-like pol and one HERVF-like pol was identified in the gamma genus, the gamma_pol.fasta summary fasta would contain:
	* The new MLV-like pol sequence
	* The new HERVF-like pol sequence
	* The reference sequences from `ERVsearch/phylogenies/group_phylogenies/MLV-like_gamma_pol.fasta` - highly related sequences from the MLV-like group
	* The reference sequences from `ERVsearch/phylogenies/group_phylogenies/HERVF-like_gamma_pol.fasta` - highly related sequences from the HERVF-like group.
	* The reference sequences from `ERVsearch/phylogenies/summary_phylogenies/gamma_pol.fasta` - a less detailed but more diverse set of gammaretroviral pol ORFs.
	* A epsilonretrovirus outgroup
 
 This ensures sufficient detail in the groups of interest while avoiding excessive detail in groups where nothing new has been identified.
 
 These FASTA files are saved as GENE_GENUS.fasta
 
The files are aligned using the MAFFT fftns algorithm https://mafft.cbrc.jp/alignment/software/manual/manual.html to generate the GENE_GENUS_A.fasta aligned output files.

### makeSummaryTrees

**Input Files**<br>
[`summary_fastas.dir/GENE_GENUS_A.fasta`](outputs.html#makesummaryfastas)<br>

**Output Files**<br>
[`summary_trees.dir/GENE_GENUS.tre`](outputs.html#makesummarytrees)<br>

**Parameters**<br>
None<br>

Builds a phylogenetic tree, using the FastTree2 algorithm (http://www.microbesonline.org/fasttree) with the default settings plus the GTR model, for the aligned group FASTA files generated by the makeSummaryFastas function.

### drawSummaryTrees

**Input Files**<br>
[`summary_trees.dir/GENE_GENUS.tre`](outputs.html#makesummarytrees)<br>

**Output Files**<br>
[`summary_trees.dir/GENE_GENUS.FMT`](outputs.html#drawsummarytrees)(FMT = png, svg, pdf or jpg)<br>

**Parameters**<br>
[`[plots] gag_colour`](parameters.html#gag-colour)<br>
[`[plots] pol_colour`](parameters.html#pol-colour)<br>
[`[plots] env_colour`](parameters.html#env-colour)<br>
[`[trees] use_gene_colour`](parameters.html#use-gene-colour)<br>
[`[trees] maincolour`](parameters.html#maincolour)<br>
[`[trees] highlightcolour`](parameters.html#highlightcolour)<br>
[`[trees] outgroupcolour`](parameters.html#outgroupcolour)<br>
[`[trees] dpi`](parameters.html#dpi)<br>
[`[trees] format`](parameters.html#format)<br>

Generates an image file for each file generated in the makeSummaryTrees step, using ete3 (http://etetoolkit.org). Newly identified sequences are labelled as "~" and shown in a different colour. Monophyletic groups of newly identified ERVs have been collapsed (by choosing a single representative sequence) and the number of sequences in the group is added to the label and represented by the size of the node tip.

By default, newly identified sequences are shown in the colours specified in `plots_gag_colour`, `plots_pol_colour` and `plots_env_colour` - to do this then `trees_use_gene_colour` should be set to True in the `pipeline.ini`. Alternatively, a fixed colour can be used by setting `trees_use_gene_colour` to False and settings `trees_highlightcolour`. The text colour of the reference sequences (default black) can be set using `trees_maincolour` and the outgroup using `trees_outgroupcolour`.

The output file DPI can be specified using `trees_dpi` and the format (which can be png, svg, pdf or jpg) using `trees_format`.

### summariseClassify

**Input Files**<br>
[`summary_fastas.dir/GENE_GENUS.fasta`](outputs.html#makesummaryfastas)<br>
[`summary_trees.dir/GENE_GENUS.tre`](outputs.html#makesummarytrees)<br>

**Output Files**<br>
[`classify_results.dir/results.tsv`](outputs.html#summariseclassify)<br>
[`classify_results.dir/by_gene_genus.FMT`](outputs.html#summariseclassify) (png, svg, pdf or jpg)<br>

**Parameters**<br>
None

Combines the results of the classify steps to generate additional summary files. The `results.tsv` output file lists the number of genes which have been collapsed into each group in the trees in the `summary_trees.dir` directory. The by_gene_genus.FMT plot is a bar chart of the same information, organised by gene and genus.


### Classify

**Input Files** None<br>

**Output Files** None<br>

**Parameters** None<br>

Helper function to run all screening functions and classification functions (all functions prior to this point).



## ERVRegions

Identifies regions of the genome containing ORFs resembling more than one different retroviral gene within a certain distance

1. [makeCleanBeds](#makecleanbeds)
2. [makeCleanFastas](#makecleanfastas)
3. [findERVRegions](#findervregions)
4. [makeRegionTables](#makeregiontables)
5. [plotERVRegions](#plotervregions)
6. [summariseERVRegions](#summariseervregions)
7. [ERVRegions](#id3)
8. [Full](#full)

### makeCleanBeds

**Input Files**<br>
[`grouped.dir/GENE_groups.tsv`](outputs.html#assigngroups)<br>

**Output Files**<br>
[`clean_beds.dir/GENE.bed`](outputs.html#makecleanbeds)<br>

**Parameters**<br>
None<br>

Generates a bed file for each gene which contains the co-ordinates of the ORFs which have passed all filtering criteria in the Screen section.


### makeCleanFastas

**Input Files**<br>
[`clean_beds.dir/GENE.bed`](outputs.html#makecleanbeds)<br>
[`genome.fa`](introduction.html#required-input-files)<br>

**Output Files**<br>
[`clean_fastas.dir/GENE.fasta`](outputs.html#makecleanfastas)<br>

**Parameters**<br>
None<br>
 
Fasta files are generated containing the sequences of the regions listed by makeCleanBeds. These are extracted from the host chromosomes using bedtools getfasta (https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html).

### findERVRegions

**Input Files**<br>
[`clean_fastas.dir/*.fasta`](outputs.html#makecleanfastas)<br>

**Output Files**<br>
[`ERV_regions.dir/all_ORFs.bed`](outputs.html#findervregions)<br>
[`ERV_regions.dir/all_regions.bed`](outputs.html#findervregions)<br>
[`ERV_regions.dir/multi_gene_regions.bed`](outputs.html#findervregions)<br>
[`ERV_regions.dir/regions.fasta`](outputs.html#findervregions)<br>

**Parameters**<br>
[`[regions] maxdist`](parameters.html#maxdist)<br>

Combines the files containng the ORF regions for the different retroviral genes and merges any regions which are within `regions_maxdist` of each other to find larger regions containing multiple genes.
The `all_ORFs.bed` output file is the concatenated and sorted bed files, `all_regions.bed `contains the merged regions with any ORFs within regions_maxdist of each other (end to end) combined, plus all regions with a single ORF, generated from all_regions.bed using bedtools merge (https://bedtools.readthedocs.io/en/latest/content/tools/merge.html). The name, strand and score columns are concatenated for merged regions, delimited with a ",".
`multi_gene_regions.bed` contains only the regions which were found to contain multiple ORFs, `regions.fasta` is the sequence of these regions in FASTA format. At this point this includes regions with multiple ORFs from the same gene (e.g. two *pol* ORFs).

### makeRegionTables

**Input Files**<br>
[`ERV_regions.dir/multi_gene_regions.bed`](outputs.html#findervregions)<br>
[`grouped.dir/*_groups.tsv`](outputs.html#assigngroups)<br>
[`genome.fa`](introduction.html#required-input-files)<br>

**Output Files**
[`ERV_regions.dir/ERV_regions_final.tsv`](outputs.html#makeregiontables)<br>
[`ERV_regions.dir/ERV_regions_final.bed`](outputs.html#makeregiontables)<br>
[`ERV_regions.dir/ERV_regions_final.fasta`](outputs.html#makeregiontables)<br>

**Parameters**<br>
[`[regions] maxoverlap`](parameters.html#maxoverlap)<br>

Takes a merged bed file consisting of regions of the genome identified as having more than one ERV-like ORF, finds the regions within this file
which contain more than one different gene (e.g. gag and pol instead of two gag ORFs) and outputs a formatted table of information about these
regions.

The output table (`ERV_regions_final.tsv`) will usually have 37 columns:
* `name` - the final ID of the ERV region - the genes found plus an integer  e.g. gag_pol_12<br>
* `chrom` - chromosome<br>
* `start` - start position of the ERV region<br>
* `end` - end position of the ERV region<br>
* `strand` - strand of the ERv region<br>
* `genus` - genus of the ERV region, can be multiple genera delimted by "|" if different genes had different genera<br>
* for each gene screened for (usually gag, pol and env)<br>
    * `GENE_name` - the names of the ORFs for this gene in this region<br>
    * `GENE_ID` - the original IDs of the ORFs for this gene in this region<br>
    * `GENE_start` - the start position of this gene in this region (genome co-ordinates)<br>
    * `GENE_relative_start` - the start position of this gene in this region (relative to the start of the region)<br>
    * `GENE_end` - the end position of this gene in this region (genome co-ordinates)<br>
    * `GENE_relative_end` - the end position of this gene in this region (relative to the start of the region)<br>
    * `GENE_strand` - the strand for this gene in this region<br>
    * `GENE_match` - the closest reference retrovirus to this gene in this region<br>
    * `GENE_group` - the group of the closest reference retrovirus to this gene in this region<br>
    * `GENE_genus` - the genus of the closest reference retrovirus to this gene in this region<br>
* `orig_name` - the name of the region in the input table<br>

If not all genes are screened for the table will not have the columns for this gene.

A bed file (`ERV_regions_final.bed`) is generated with the co-ordinates of the identified regions and a FASTA file (`ERV_regions_final.fasta`) containing their sequences.


### plotERVRegions
**Input_Files**<br>
[`ERV_regions.dir/ERV_regions_final.tsv`](outputs.html#makeregiontables)<br>

**Output_Files**<br>
[`ERV_region_plots.dir/*FMT`](outputs.html#plotervregions)<br>

**Parameters**<br>
[`[plots] format`](parameters.html#format)<br>
[`[plots] dpi`](parameters.html#dpi)<br>
[`[plots] gag_colour`](parameters.html#gag-colour)<br>
[`[plots] pol_colour`](parameters.html#pol-colour)<br>
[`[plots] env_colour`](parameters.html#env-colour)<br>

For each region containing ORFs resembling more than one retroviral gene, a plot is generated showing how these ORFs are distributed on the genome relative to each other.

Each gene is shown on a different line on the y axis, the x axis is chromosome co-ordinates.
  
### summariseERVRegions

**Input Files**<br>
[`ERV_regions.dir/ERV_regions_final.tsv`](outputs.html#makeregiontables)<br>


**Output Files**<br>
[`erv_regions_results.dir/results.tsv`](introduction.html#id3)<br>
[`erv_regions_results.dir/erv_regions.FMT`](introduction.html#id3<br>
	
**Parameters**<br>
[`[plots] other_colour`](parameters.html#other-colour)<br>


Combines the results of the ERVregions steps to generate additional summary files.

The `results.tsv` output file is a copy of the output of the makeRegionTables functions.
This will usually have 37 columns:
* `name` - the final ID of the ERV region - the genes found plus an integer  e.g. gag_pol_12<br>
* `chrom` - chromosome<br>
* `start` - start position of the ERV region<br>
* `end` - end position of the ERV region<br>
* `strand` - strand of the ERv region<br>
* `genus` - genus of the ERV region, can be multiple genera delimted by "|" if different genes had different genera<br>
* for each gene screened for (usually gag, pol and env)<br>
    * `GENE_name` - the names of the ORFs for this gene in this region<br>
    * `GENE_ID` - the original IDs of the ORFs for this gene in this region<br>
    * `GENE_start` - the start position of this gene in this region (genome co-ordinates)<br>
    * `GENE_relative_start` - the start position of this gene in this region (relative to the start of the region)<br>
    * `GENE_end` - the end position of this gene in this region (genome co-ordinates)<br>
    * `GENE_relative_end` - the end position of this gene in this region (relative to the start of the region)<br>
    * `GENE_strand` - the strand for this gene in this region<br>
    * `GENE_match` - the closest reference retrovirus to this gene in this region<br>
    * `GENE_group` - the group of the closest reference retrovirus to this gene in this region<br>
    * `GENE_genus` - the genus of the closest reference retrovirus to this gene in this region<br>
* `orig_name` - the name of the region in the input table<br>

A bar chart - `erv_regions_results.dir/erv_regions.FMT` is also generated showing the number of ERV regions found with each combination of genes.


### ERVRegions

**Input Files**
None

**Output Files**
None

**Parameters**
None

Helper function to run all screening functions and ERVRegions functions.


### Full

**Input Files**
None

**Output Files**
None

**Parameters**
None

Helper function to run all functions.