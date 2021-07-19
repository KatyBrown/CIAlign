for n in {1..10}
do
    for file in reference\_set\_$n/input/*fasta
    do
	for stri in highly_stringent med_stringent low_stringent;
	do
	    stem=$(echo $file | cut -d "." -f1 | cut -d "/" -f3)

	    inf=reference\_set\_$n/clustal/$stem\_auto.fasta
	    outf=reference\_set\_$n/clustal/$stri\_$stem\_auto
	    CIAlign --inifile $stri\_config\_BB.ini --infile $inf --outfile_stem $outf

	    inf=reference\_set\_$n/muscle/$stem\_max100.fasta
	    outf=reference\_set\_$n/muscle/$stri\_$stem\_max100
	    CIAlign --inifile $stri\_config\_BB.ini --infile $inf --outfile_stem $outf

	    inf=reference\_set\_$n/mafft/$stem\_local\_max100.fasta
	    outf=reference\_set\_$n/mafft/$stri\_$stem\_local\_max100
	    CIAlign --inifile $stri\_config\_BB.ini --infile $inf --outfile_stem $outf

	    inf=reference\_set\_$n/mafft/$stem\_global\_max100.fasta
	    outf=reference\_set\_$n/mafft/$stri\_$stem\_global\_max100
	    CIAlign --inifile $stri\_config\_BB.ini --infile $inf --outfile_stem $outf	    
	done
    done
done

