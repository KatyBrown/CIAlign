


for i in {1..100};
do
    clustalo --auto -i sim\_$i/human_gapdh_Unaligned.FASTA -o sim\_$i/clustal/nucleotide/auto.fasta
    mafft --localpair --maxiterate 100 sim\_$i/human_gapdh_Unaligned.FASTA > sim\_$i/mafft/nucleotide/local_max100.fasta
    mafft --maxiterate 100 sim\_$i/human_gapdh_Unaligned.FASTA > sim\_$i/mafft/nucleotide/global_max100.fasta
    muscle -maxiters 100 -in sim\_$i/human_gapdh_Unaligned.FASTA -out sim\_$i/muscle/nucleotide/max100.fasta
done
