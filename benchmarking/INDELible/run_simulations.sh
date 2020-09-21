for typ in nucleotide amino_acid codon;
do
    mkdir $typ
    cp control\_$typ.txt control.txt
    for j in {1..100};
    do
	mkdir $typ/sim\_$j
	indelible
	for i in {2..9}; do
	    awk -v i=$i '{if(NR % 11 == i && x == 0)
                            {printf(">%s\n%s", $1, $2); x = x + 1}
                          else if (NR % 11 == i)
                            {printf($2)}}; END{printf("\n")}' currentOut_TRUE.phy;
	done > $typ/sim\_$j/alignment_true.fasta
    sed 's/-//g' $typ/sim\_$j/alignment\_true.fasta > $typ/sim\_$j/sequences.fasta
    rm currentOut*
    mv LOG.txt $typ/sim\_$j/log.txt
    awk 'NR == 7' trees.txt  | cut -f9 > $typ/sim\_$j/tree.nwk
    rm trees.txt
    done
done
