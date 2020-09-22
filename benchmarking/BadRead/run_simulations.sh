
for i in {1..100};
do
    mkdir sim\_$i
    badread simulate --reference dwv.fasta --quantity 10x --length 15000,13000 --identity 95,100,4 --start_adapter_seq "" --end_adapter_seq "" --chimeras 0 --glitches 0,0,0 --junk_reads 0 --random_reads 0 | paste - - - - | awk -F "\t" '{printf(">%s\n%s\n", $1, $2)}' > sim\_$i/good_nanopore.fasta

    badread simulate --reference dwv.fasta --quantity 10x --length 7500,7500 --identity 95,100,4 --start_adapter_seq "" --end_adapter_seq "" --chimeras 0 --glitches 0,0,0 --junk_reads 0 --random_reads 0 --error_model pacbio --qscore_model pacbio | paste - - - - | awk -F "\t" '{printf(">%s\n%s\n", $1, $2)}'  > sim\_$i/good_pacbio.fasta

    badread simulate --reference dwv.fasta --quantity 10x --length 15000,13000 --identity 85,95,5 --chimeras 1 --glitches 10000,25,25 --junk_reads 1 --random_reads 1 | paste - - - - | awk -F "\t" '{printf(">%s\n%s\n", $1, $2)}'  > sim\_$i/medium_nanopore.fasta

    badread simulate --reference dwv.fasta --quantity 10x --length 7500,7500 --identity 85,95,5 --chimeras 1 --glitches 10000,25,25 --junk_reads 1 --random_reads 1 --error_model pacbio --qscore_model pacbio  | paste - - - - | awk -F "\t" '{printf(">%s\n%s\n", $1, $2)}' > sim\_$i/medium_pacbio.fasta

    badread simulate --reference dwv.fasta --quantity 10x --length 15000,13000 --identity 75,90,8 --chimeras 10 --glitches 1000,100,100 --junk_reads 5 --random_reads 5 | paste - - - - | awk -F "\t" '{printf(">%s\n%s\n", $1, $2)}'  > sim\_$i/bad_nanopore.fasta

    badread simulate --reference dwv.fasta --quantity 10x --length 7500,7500 --identity 75,90,8 --chimeras 10 --glitches 1000,100,100 --junk_reads 5 --random_reads 5 --error_model pacbio --qscore_model pacbio | paste - - - - | awk -F "\t" '{printf(">%s\n%s\n", $1, $2)}'  > sim\_$i/bad_pacbio.fasta
done
