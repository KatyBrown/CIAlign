for i in {1..100};
do	
    for stri in highly_stringent med_stringent low_stringent;
    do
	
            CIAlign --inifile $stri\_config\_INDELible.ini --infile nucleotide/sim\_$i/mafft/local_max100.fasta --outfile_stem nucleotide/sim\_$i/mafft/$stri\_local_max100
            CIAlign --inifile $stri\_config\_INDELible.ini --infile nucleotide/sim\_$i/mafft/global_max100.fasta --outfile_stem nucleotide/sim\_$i/mafft/$stri\_global_max100
            CIAlign --inifile $stri\_config\_INDELible.ini --infile nucleotide/sim\_$i/muscle/max100.fasta --outfile_stem nucleotide/sim\_$i/muscle/$stri\_max100
            CIAlign --inifile $stri\_config\_INDELible.ini --infile nucleotide/sim\_$i/clustal/auto.fasta --outfile_stem nucleotide/sim\_$i/clustal/$stri\_auto
	done
    done
done

