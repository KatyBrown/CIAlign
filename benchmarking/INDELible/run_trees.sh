for i in {1..100};
do
    FastTree -gtr -nt nucleotide/sim\_$i/clustal/auto.fasta > nucleotide/sim\_$i/clustal/auto.tre
    FastTree -gtr -nt nucleotide/sim\_$i/muscle/max100.fasta > nucleotide/sim\_$i/muscle/max100.tre
    FastTree -gtr -nt nucleotide/sim\_$i/mafft/local\_max100.fasta > nucleotide/sim\_$i/mafft/auto.tre
    FastTree -gtr -nt nucleotide/sim\_$i/mafft/global_max100.fasta > nucleotide/sim\_$i/mafft/global\_max100.tre
    for stri in highly_stringent med_stringent low_stringent;
    do
	FastTree -gtr -nt nucleotide/sim\_$i/clustal/$stri\_auto\_cleaned.fasta > nucleotide/sim\_$i/clustal/$stri\_auto.tre
	FastTree -gtr -nt nucleotide/sim\_$i/muscle/$stri\_max100\_cleaned.fasta > nucleotide/sim\_$i/muscle/$stri\_max100.tre
	FastTree -gtr -nt nucleotide/sim\_$i/mafft/$stri\_local\_max100\_cleaned.fasta > nucleotide/sim\_$i/mafft/$stri\_local\_max100.tre
	FastTree -gtr -nt nucleotide/sim\_$i/mafft/$stri\_global\_max100\_cleaned.fasta > nucleotide/sim\_$i/mafft/$stri\_global\_max100.tre
    done
done
