

## ERVRegions
---
### makeCleanBeds

**Input Files**<br>
`grouped.dir/GENE_groups.tsv`<br>

**Output Files**<br>
`clean_beds.dir/GENE.bed`<br>

**Parameters**<br>
None<br>

Generates a bed file for each gene which contains the co-ordinates of the ORFs which have passed all filtering criteria in the Screen section.


### makeCleanFastas

**Input Files**<br>
`clean_beds.dir/GENE.bed`<br>
`genome.fa`<br>

**Output Files**<br>
`clean_fastas.dir/GENE.fasta`<br>

**Parameters**<br>
None<br>
 
Fasta files are generated containing the sequences of the regions listed by makeCleanBeds. These are extracted from the host chromosomes using bedtools getfasta (https://bedtools.readthedocs.io/en/latest/content/tools/getfasta.html).

### findERVRegions

**Input Files**<br>
`clean_fastas.dir/*.fasta`<br>

**Output Files**<br>
`ERV_regions.dir/all_ORFs.bed`<br>
`ERV_regions.dir/all_regions.bed`<br>
`ERV_regions.dir/multi_gene_regions.bed`<br>
`ERV_regions.dir/regions.fasta`<br>

**Parameters**<br>
`[regions] maxdist`<br>

Combines the files containng the ORF regions for the different retroviral genes and merges any regions which are within `regions_maxdist` of each other to find larger regions containing multiple genes.
The `all_ORFs.bed` output file is the concatenated and sorted bed files, `all_regions.bed `contains the merged regions with any ORFs within regions_maxdist of each other (end to end) combined, plus all regions with a single ORF, generated from all_regions.bed using bedtools merge (https://bedtools.readthedocs.io/en/latest/content/tools/merge.html). The name, strand and score columns are concatenated for merged regions, delimited with a ",".
`multi_gene_regions.bed` contains only the regions which were found to contain multiple ORFs, `regions.fasta` is the sequence of these regions in FASTA format. At this point this includes regions with multiple ORFs from the same gene (e.g. two *pol* ORFs).

### makeRegionTables

**Input Files**<br>
`ERV_regions.dir/multi_gene_regions.bed`<br>
`grouped.dir/*_groups.tsv`<br>
`genome.fa`<br>

**Output Files**
`ERV_regions.dir/ERV_regions_final.tsv`<br>
`ERV_regions.dir/ERV_regions_final.bed`<br>
`ERV_regions.dir/ERV_regions_final.fasta`<br>

**Parameters**<br>
None<br>

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


#### summariseERVRegions

**Input Files**

**Output Files**

**Parameters**


### ERVRegions

**Input Files**

**Output Files**

**Parameters**


### Full