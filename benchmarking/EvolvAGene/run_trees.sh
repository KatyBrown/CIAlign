
for i in {1..100};
do
    FastTree sim\_$i/clustal/pep/auto.fasta > sim\_$i/clustal/pep/auto.tre
    FastTree sim\_$i/mafft/pep/local_max100.fasta > sim\_$i/mafft/pep/local_max100.tre
    FastTree sim\_$i/mafft/pep/global_max100.fasta > sim\_$i/mafft/pep/global_max100.tre
    FastTree sim\_$i/muscle/pep/max100.fasta > sim\_$i/muscle/pep/max100.tre

    for stri in highly_stringent med_stringent low_stringent;
    do
	FastTree sim\_$i/clustal/pep/$stri\_auto_cleaned.fasta > sim\_$i/clustal/pep/$stri\_auto.tre
	FastTree sim\_$i/mafft/pep/$stri\_local_max100_cleaned.fasta > sim\_$i/mafft/pep/$stri\_local_max100.tre
	FastTree sim\_$i/mafft/pep/$stri\_global_max100_cleaned.fasta > sim\_$i/mafft/pep/$stri\_global_max100.tre
	FastTree sim\_$i/muscle/pep/$stri\_max100_cleaned.fasta > sim\_$i/muscle/pep/$stri\_max100.tre
    done
done
