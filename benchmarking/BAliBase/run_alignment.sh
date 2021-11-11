for i in {1..10}
do
    for file in reference\_set\_$i/input/*;
    do
	stem=$(echo $file | rev | cut -d "/" -f1 | cut -d "." -f2 | rev)
        clustalo --auto -i reference\_set\_$i/input/$stem.fasta -o reference\_set\_$i/clustal/$stem\_auto.fasta
	mafft --localpair --maxiterate 100 reference\_set\_$i/input/$stem.fasta > reference\_set\_$i/mafft/$stem\_local\_max100.fasta
	mafft --maxiterate 100 reference\_set\_$i/input/$stem.fasta > reference\_set\_$i/mafft/$stem\_global\_max100.fasta
	muscle -maxiters 100 -in reference\_set\_$i/input/$stem.fasta > reference\_set\_$i/muscle/$stem\_max100.fasta	
    done
done
