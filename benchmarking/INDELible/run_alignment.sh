
for tool in clustal muscle mafft;
do
    for i in {1..100};
    do
	mkdir nucleotide/sim\_$i/$tool
    done
done

for i in {1..100};
do
    clustalo --auto -i nucleotide/sim\_$i/sequences.fasta -o nucleotide/sim\_$i/clustal/auto.fasta
    mafft --localpair --maxiterate 100 nucleotide/sim\_$i/sequences.fasta > nucleotide/sim\_$i/mafft/local_max100.fasta
    mafft --maxiterate 100 nucleotide/sim\_$i/sequences.fasta > nucleotide/sim\_$i/mafft/global_max100.fasta
    muscle -maxiters 100 -in nucleotide/sim\_$i/sequences.fasta -out nucleotide/sim\_$i/muscle/max100.fasta
done
