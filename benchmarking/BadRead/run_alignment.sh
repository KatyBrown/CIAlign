for i in {1..100}
do
    mkdir sim\_$i/mafft
    for stem in good_nanopore medium_nanopore bad_nanopore;
    do
	for i in {1..100};
	do
	    mafft --localpair --adjustdirection sim\_$i/$stem.fasta > sim\_$i/mafft/$stem\_mafft_localpair.fasta;
	done;
    done
done

