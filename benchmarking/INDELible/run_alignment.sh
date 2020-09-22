for n in nucleotide amino_acid codon;
do
    
    for tool in clustal muscle mafft;
    do
	for i in {1..100};
	do
	    mkdir $n/sim\_$i/$tool
	done
    done
done



for i in {1..100};
do
    for n in nucleotide amino_acid codon;
	     do
		 clustalo --auto -i $n/sim\_$i/sequences.fasta -o $n/sim\_$i/clustal/auto.fasta
		 mafft --localpair --maxiterate 100 $n/sim\_$i/sequences.fasta > $n/sim\_$i/mafft/local_max100.fasta
		 mafft --maxiterate 100 $n/sim\_$i/sequences.fasta > $n/sim\_$i/mafft/global_max100.fasta
		 muscle -maxiters 100 -in $n/sim\_$i/sequences.fasta -out $n/sim\_$i/muscle/max100.fasta
    done
done
