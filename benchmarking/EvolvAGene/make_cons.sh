for i in {1..100}
do
   CIAlign --infile sim\_$i/human_gapdh_True_alignment.FASTA --make_consensus --consensus_type majority_nongap --outfile_stem sim\_$i/true\_consensus
done
