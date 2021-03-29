for i in {1..100};
do
    for typ in good medium bad;
    do
	ali=sim\_$i/mafft/$typ\_nanopore\_mafft\_localpair.fasta
	CIAlign --infile $ali --outfile_stem sim\_$i/mafft/$typ\_nanopore\_local\_max100 --make_consensus --consensus_type majority_nongap;
    done
done

