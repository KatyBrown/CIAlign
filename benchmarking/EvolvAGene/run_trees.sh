
for i in {1..100};
do
    FastTree sim\_$i/clustal/nucleotide/auto.fasta > sim\_$i/clustal/nucleotide/auto.tre
    FastTree sim\_$i/mafft/nucleotide/local_max100.fasta > sim\_$i/mafft/nucleotide/local_max100.tre
    FastTree sim\_$i/mafft/nucleotide/global_max100.fasta > sim\_$i/mafft/nucleotide/global_max100.tre
    FastTree sim\_$i/muscle/nucleotide/max100.fasta > sim\_$i/muscle/nucleotide/max100.tre

    for stri in highly_stringent med_stringent low_stringent;
    do
	FastTree sim\_$i/clustal/nucleotide/$stri\_auto_cleaned.fasta > sim\_$i/clustal/nucleotide/$stri\_auto.tre
	FastTree sim\_$i/mafft/nucleotide/$stri\_local_max100_cleaned.fasta > sim\_$i/mafft/nucleotide/$stri\_local_max100.tre
	FastTree sim\_$i/mafft/nucleotide/$stri\_global_max100_cleaned.fasta > sim\_$i/mafft/nucleotide/$stri\_global_max100.tre
	FastTree sim\_$i/muscle/nucleotide/$stri\_max100_cleaned.fasta > sim\_$i/muscle/nucleotide/$stri\_max100.tre
    done
done
