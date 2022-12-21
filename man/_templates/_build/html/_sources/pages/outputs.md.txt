# Minor Output Files

1.[Screen](#screen)<br>
2.[Classify](#classify)<br>
3.[ERVRegions](#ervregions)<br>

## Screen
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
`init.txt`<br>
Placeholder file to show initial checks have been run.<br>

### genomeToChroms
`host_chromosomes.dir/*fasta`<br>
FASTA format files containing the input genome (or other sequence), divided into regions based on the [genomesplits](parameters.html#genomesplits) parameters.<br>

### prepDBs
`gene_databases.dir/GENE.fasta`<br>
Copies of the fasta files of reference retroviral amino acid sequences.<br>

### runExonerate
`raw_exonerate_output.dir/GENE_*.tsv`<br>
Raw "vulgar" output of Exonerate as described [here](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual)<br>

### cleanExonerate
`clean_exonerate_output.dir/GENE_*_unfiltered.tsv`<br>
Raw exonerate output converted into a table.<br>
Columns:<br>
* **query_id**<br>
Reference amino acid sequence ID for retroviral gene<br>
* **query_start**<br>
Start position of match within reference sequence<br>
* **query_end**<br>
End position of match within reference sequence<br>
* **query_strand**<br>
Strand of match relative to reference sequence<br>
* **target_id**<br>
Nucleotide sequence from input which matched the retroviral gene.<br>
* **target_start**<br>
Start position of match within input sequence.<br>
* **target_end**<br>
End position of match within output sequence.<br>
* **target_strand**<br>
Strand of match relative to input sequence<br>
* **score**<br>
Exonerate score of match
* **details**<br>
Additional columns (columns 11+) from [Exonerate vulgar output](https://www.ebi.ac.uk/about/vertebrate-genomics/software/exonerate-manual), delimited by "|".
* **length**<br>
Length of match on input sequence<br>

`clean_exonerate_output.dir/GENE_*_filtered.tsv`<br>
Output of Exonerate filtered to remove regions containing introns and regions which are shorter than [exonerate_min_hit_length](parameters.html#min-hit-length). Columns are the same as in the unfiltered table.

`clean_exonerate_output.dir/GENE_*.bed`<br>
ERV-like regions from the table above in [bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1)<br>

### mergeOverlaps
`gene_bed_files.dir/GENE_all.bed`,<br>
ERV-like regions from `clean_exonerate_output.dir/GENE_*.bed` combined into a single [bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) file.<br>

`gene_bed_files.dir/GENE_merged.bed`<br>
ERV-like regions from the previous file merged using [Bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html)<br>

### makeFastas
`gene_fasta_files.dir/GENE_merged.fasta`<br>
Fasta files of the merged regions in `gene_bed_files.dir/GENE_merged.bed`<br>

### renameFastas
`gene_fasta_files.dir/GENE_merged_renamed.fasta`<br>
Fasta files from the `makeFastas` step with the ERV-like regions renamed with unique IDs.<br>
 
### makeUBLASTDb
`UBLAST_db.dir/GENE_db.udb`<br>
UBLAST databases for the reference retroviral amino acid sequences.

### runUBLASTCheck
`ublast.dir/GENE_UBLAST_alignments.txt`<br>
Raw UBLAST output files for the Exonerate regions vs the retroviral amino acid databases. Equivalent to the BLAST [pairwise](https://www.ncbi.nlm.nih.gov/books/NBK279684/) output.<br>

`ublast.dir/GENE_UBLAST.tsv`<br>
UBLAST tabular output for Exonerate regions vs the retrovrial amino acid databases. Equivalent to the BLAST [tabular](https://www.ncbi.nlm.nih.gov/books/NBK279684/) output.<br>

`ublast.dir/GENE_filtered_UBLAST.fasta`<br>
Fasta file of the regions which passed the UBLAST filter.<br>

### classifyWithExonerate
`exonerate_classification.dir/GENE_all_matches_exonerate.tsv`<br>
Raw output of ungapped Exonerate algorithm for the UBLAST verified regions against the ERV amino acid database.<br>
Columns:
* ID of newly identified ERV-like region
* ID of reference ERV amino acid
* Exonerate score

`exonerate_classification.dir/GENE_best_matches_exonerate.tsv`<br>
The previous table filtered to list only the highest scoring hit for each ERV-like region.<br>

`exonerate_classification.dir/GENE_refiltered_matches_exonerate.fasta`<br>
Highest scoring hits from the previous table in FASTA format.<br>

### getORFs
`ORFs.dir/GENE_orfs_raw.fasta`<br>
Raw output of running [EMBOSS transeq -frame 6](http://emboss.sourceforge.net/apps/release/6.6/emboss/apps/transeq.html) on the regions output from `classifyWithExonerate`.<br>

`ORFs.dir/GENE_orfs_nt.fasta`<br>
ORFs longer than [ORFs_min_orf_length](parameters.html#min-orf-length) as nucleotide sequences, with IDs redefined to include the chromosome, start and end position and strand of the ORF.

`ORFs.dir/GENE_orfs_aa.fasta`<br>
ORFs longer than [ORFs_min_orf_length](parameters.html#min-orf-length) as amino acid sequences, with IDs redefined to include the chromosome, start and end position and strand of the ORF.


### checkORFsUBLAST
`ublast_orfs.dir/GENE_UBLAST_alignments.txt`<br>
Raw UBLAST output files for the newly identified ORFs vs the retroviral amino acid databases. Equivalent to the BLAST [pairwise](https://www.ncbi.nlm.nih.gov/books/NBK279684/) output.<br>

`ublast_orfs.dir/GENE_UBLAST.tsv`<br>
UBLAST tabular output for newly identified ORFs vs the retrovrial amino acid databases. Equivalent to the BLAST [tabular](https://www.ncbi.nlm.nih.gov/books/NBK279684/) output.<br>

`ublast_orfs.dir/GENE_filtered_UBLAST.fasta`<br>
Fasta file of the newly identified ORFs which passed the UBLAST filter.<br>

### assignGroups
`grouped.dir/GENE_groups.tsv`<br>
Summarised output for the previous steps. This is identical to [`screen_results.dir/results.tsv`](introduction.html#id1).

### summariseScreen
FMT can be png, svg, pdf, jpg depending on the `plot_format` parameter

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

## Classify
1. [makeGroupFastas](#makegroupfastas)<br>
2. [makeGroupTrees](#makegrouptrees)<br>
3. [drawGroupTrees](#drawgrouptrees)<br>
4. [makeSummaryFastas](#makesummaryfastas)<br>
5. [makeSummaryTrees](#makesummarytrees)<br>
6. [drawSummaryTrees](#drawsummarytrees)<br>
7. [summariseClassify](#summariseclassify)<br>
8. [Classify](#id2)

### makeGroupFastas
`group_fastas.dir/GENE_(.*)_GENUS.fasta`<br>
Fasta files for each small subgroup of ERV-like ORFs and reference sequences.<br>

`group_fastas.dir/GENE_(.*)_GENUS_A.fasta`<br>
Aligned version of the above Fasta file generated using [MAFFT](https://mafft.cbrc.jp/alignment/software/)<br>

### makeGroupTrees
`group_trees.dir/GENE_(.*)_GENUS.tre<br>
Phylogenetic trees in Newick format for each small subgroup of ERV-like ORFs and reference sequences generated using [FastTree](http://www.microbesonline.org/fasttree/)<br>

### drawGroupTrees
`group_trees.dir/GENE_(.*_)GENUS.FMT` (png, svg, pdf or jpg)<br>
Image files of the phylogenetic trees for each small subgroup of ERV-like ORFs and reference sequences, with newly identified sequences highlighted.<br>

### makeSummaryFastas
`summary_fastas.dir/GENE_GENUS.fasta`<br>
Fasta files for each retroviral gene and genus combining reference and newly identified sequences. Monophyletic groups of newly identified ORfs sequences are represented by a single sequence.<br>

`group_lists.dir/*tsv`
Lists of sequences in each monophyletic group of newly identified ORFs which has been collated to be represented by a single sequnce.<br>

### makeSummaryTrees
`summary_trees.dir/GENE_GENUS.tre`<br>
Tree files in Newick format for each retroviral gene and genus combining reference and newly identified sequences. Monophyletic groups of newly identified ORfs sequences are represented by a single sequence.<br>

### drawSummaryTrees
`summary_trees.dir/GENE_GENUS.FMT` (FMT = png, svg, pdf or jpg)<br>
Images files of the phylogenetic trees for each retroviral gene and genus. Different sized circles are used to show the relative size of collapsed monophyletic groups. Newly identified ERVs are highlighted.<br>

### summariseClassify
`classify_results.dir/results.tsv`<br>
Table listing the number of genes which have been collapsed into each monophyletic group in the trees in the `summary_trees.dir` directory.<br>
Columns:
* **gene** Retroviral gene for this group
* **genus** Retroviral genus for this group
* **group** Group ID
* **count** Number of sequences in this group

`classify_results.dir/by_gene_genus.png`<br>
Bar chart showing the number of genes which have been collapsed into each monophyletic group, organised by gene and genus.

## ERVRegions
1. [makeCleanBeds](#makecleanbeds)
2. [makeCleanFastas](#makecleanfastas)
3. [findERVRegions](#findervregions)
4. [makeRegionTables](#makeregiontables)
5. [summariseERVRegions](#summariseervregions)
6. [plotERVRegions](#plotervregions)
7. [ERVRegions](#id3)
8. [Full](#full)


### makeCleanBeds
`clean_beds.dir/GENE.bed`<br>
[Bed](https://genome.ucsc.edu/FAQ/FAQformat.html#format1) files containing the co-ordinates all the ERV-like ORFs output by the Screen section.

### makeCleanFastas
`clean_fastas.dir/GENE.fasta`<br>
FASTA files of the ERV-like ORFs output by the Screen section.

### findERVRegions
`ERV_regions.dir/all_ORFs.bed`<br>
Concatenated version of the bed files output by the `makeCleanBeds` function with all three genes in a single file.<br>

`ERV_regions.dir/all_regions.bed`<br>
Previous bed file merged using [Bedtools merge](https://bedtools.readthedocs.io/en/latest/content/tools/merge.html), with regions within [regions_maxdist](parameters.html#maxdist) merged.<br>

`ERV_regions.dir/multi_gene_regions.bed`<br>
Filtered version of the previous bed file with only regions which were merged.<br>

`ERV_regions.dir/regions.fasta`<br>
Fasta file of the merged regions.<br>

### makeRegionTables
`ERV_regions.dir/ERV_regions_final.tsv`<br>
Table showing the details of the combined regions containing ORFs resembling more than one retroviral gene. This table is identical to [`erv_region_results.dir/results.tsv`](introduction.html#id3).<br<

`ERV_regions.dir/ERV_regions_final.bed`<br>
Bed file with the co-ordinates of the identified regions.

`ERV_regions.dir/ERV_regions_final.fasta`<br>
FASTA file of the regions in the bed file above.

### plotERVRegions
`ERV_region_plots.dir/*.FMT`<br>
Plots showing the distributions of ORFs resembling each retroviral gene on the genome. Each gene is shown on a different line on the y axis, the x axis is chromosome co-ordinates. One plot is generated for each multi-gene region.

### summariseERVRegions
`erv_region_results.dir/results.tsv`<br>
Table showing the overall results for regions with ORFs resembling multiple retroviral genes. This table is described in full in the [Main Output Files](introduction.html#main-output-files) section.

`erv_region_results.dir/erv_regions.png`<br>
Bar chart showing the number of ERV regions identified with each combination of retroviral genes.