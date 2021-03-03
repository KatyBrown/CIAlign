for i in {1..100};
do
    for stri in highly_stringent med_stringent low_stringent;
    do
	
        CIAlign --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/mafft/nucleotide/local_max100.fasta --outfile_stem sim\_$i/mafft/nucleotide/$stri\_local_max100
        CIAlign.py --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/mafft/nucleotide/global_max100.fasta --outfile_stem sim\_$i/mafft/nucleotide/$stri\_global_max100
        CIAlign.py --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/muscle/nucleotide/max100.fasta --outfile_stem sim\_$i/muscle/nucleotide/$stri\_max100
        CIAlign.py --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/clustal/nucleotide/auto.fasta --outfile_stem sim\_$i/clustal/nucleotide/$stri\_auto
        CIAlign.py --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/mafft/pep/local_max100.fasta --outfile_stem sim\_$i/mafft/pep/$stri\_local_max100
        CIAlign.py --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/mafft/pep/global_max100.fasta --outfile_stem sim\_$i/mafft/pep/$stri\_global_max100
        CIAlign.py --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/muscle/pep/max100.fasta --outfile_stem sim\_$i/muscle/pep/$stri\_max100
        CIAlign.py --inifile $stri\_config\_EvolvAGene.ini --infile sim\_$i/clustal/pep/auto.fasta --outfile_stem sim\_$i/clustal/pep/$stri\_auto
	
    done
done

