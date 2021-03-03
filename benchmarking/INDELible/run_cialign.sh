for i in {1..100};
do
    for typ in nucleotide amino_acid codon;
    do
	
	for stri in highly_stringent med_stringent low_stringent;
	do
	
            CIAlign --inifile $stri\_config\_INDELible.ini --infile $typ/sim\_$i/mafft/local_max100.fasta --outfile_stem $typ/sim\_$i/mafft/$stri\_local_max100
            CIAlign --inifile $stri\_config\_INDELible.ini --infile $typ/sim\_$i/mafft/global_max100.fasta --outfile_stem $typ/sim\_$i/mafft/$stri\_global_max100
            CIAlign --inifile $stri\_config\_INDELible.ini --infile $typ/sim\_$i/muscle/max100.fasta --outfile_stem $typ/sim\_$i/muscle/$stri\_max100
            CIAlign --inifile $stri\_config\_INDELible.ini --infile $typ/sim\_$i/clustal/auto.fasta --outfile_stem $typ/sim\_$i/clustal/$stri\_auto
	done
    done
done

