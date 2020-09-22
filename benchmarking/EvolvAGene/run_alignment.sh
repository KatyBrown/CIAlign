for tool in clustal muscle mafft;
do
    for i in {1..100};
    do
	mkdir sim\_$i/$tool
	for typ in nucleotide pep
	do
	    mkdir sim\_$i/$tool/$typ
	done
    done
done


for i in {1..100};
do
    clustalo --auto -i sim\_$i/human_gapdh_Unaligned.FASTA -o sim\_$i/clustal/nucleotide/auto.fasta
    clustalo --auto -i sim\_$i/human_gapdh_pep_Unaligned.FASTA -o sim\_$i/clustal/pep/auto.fasta
    mafft --localpair --maxiterate 100 sim\_$i/human_gapdh_Unaligned.FASTA > sim\_$i/mafft/nucleotide/local_max100.fasta
    mafft --localpair --maxiterate 100 sim\_$i/human_gapdh_pep_Unaligned.FASTA > sim\_$i/mafft/pep/local_max100.fasta
    mafft --maxiterate 100 sim\_$i/human_gapdh_Unaligned.FASTA > sim\_$i/mafft/nucleotide/global_max100.fasta
    mafft --maxiterate 100 sim\_$i/human_gapdh_pep_Unaligned.FASTA > sim\_$i/mafft/pep/global_max100.fasta
    muscle -maxiters 100 -in sim\_$i/human_gapdh_Unaligned.FASTA -out sim\_$i/muscle/nucleotide/max100.fasta
    muscle -maxiters 100 -in sim\_$i/human_gapdh_pep_Unaligned.FASTA -out sim\_$i/muscle/pep/max100.fasta
done
