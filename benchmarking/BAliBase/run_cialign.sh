for i in {1..10};
do
    for ref in reference\_set\_$i/aligned/*fasta;
    do
	stem=$(echo $ref | rev | cut -d "/" -f1 | rev | cut -d "." -f1)
        base=reference\_set\_$i
	clustal=$base/clustal/$stem\_auto.fasta
	mafftloc=$base/mafft/$stem\_local\_max100.fasta
	mafftglob=$base/mafft/$stem\_global\_max100.fasta
	muscle=$base/muscle/$stem\_max100.fasta
        for stri in high med low;
	do
	    CIAlign --inifile ~/CIAlign/benchmarking/BAliBase/$stri\_config\_BAliBase.ini --infile $clustal --outfile_stem $base/clustal/$stri\_stringency\_$stem\_auto
	    CIAlign --inifile ~/CIAlign/benchmarking/BAliBase/$stri\_config\_BAliBase.ini --infile $mafftloc --outfile_stem $base/mafft/$stri\_stringency\_$stem\_local\_max100
	    CIAlign --inifile ~/CIAlign/benchmarking/BAliBase/$stri\_config\_BAliBase.ini --infile $mafftglob --outfile_stem $base/mafft/$stri\_stringency\_$stem\_global\_max100
	    CIAlign --inifile ~/CIAlign/benchmarking/BAliBase/$stri\_config\_BAliBase.ini --infile $muscle --outfile_stem $base/muscle/$stri\_stringency\_$stem\_max100
	done
	    
    done
done
