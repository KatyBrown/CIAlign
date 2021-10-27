#!/bin/bash

FILES=$1'*'
# DESTDIR_ALI="/home/charlotte/homfam/aligned/"
DESTDIR_CIALIGN="/Users/lotti/CIAlign/benchmarking/QuanTest2/CIAlign_run/"

for f in $FILES
  do
    NAME=${f##*/}      # remove part before last slash
    NAME=${NAME%.*}    # from the new var remove the part after the last period
    SS="/Users/lotti/CIAlign/benchmarking/QuanTest2/SS/"$NAME".ss"
    aligned=$DESTDIR_CIALIGN$NAME"_muscle.fasta"
    OUTFILE_CIAlign=$DESTDIR_CIALIGN$NAME
    FAKE_MSA=$OUTFILE_CIAlign"_fake.fasta"
    LOG=$DESTDIR_CIALIGN$NAME"_log.txt"

    echo $f
    # echo $NAME >> $LOG
    # muscle -in $f -out $aligned
    # python3 /Users/lotti/CIAlign/CIAlign/CIAlign.py --infile $aligned --outfile $OUTFILE_CIAlign --clean --remove_divergent_minperc 0.1 --visualise
    # python3 /Users/lotti/CIAlign/benchmarking/functions/adjust_msa_to_quantest2.py $aligned $OUTFILE_CIAlign"_removed.txt" $FAKE_MSA
    # echo "with cialign" >> $LOG
    # echo $SS
    cd ..
    # python3 /Users/lotti/CIAlign/benchmarking/QuanTest2/quantest2_cialign.py ct518@cam.ac.uk $OUTFILE_CIAlign"_cleaned.fasta" $SS $FAKE_MSA
    # echo "without cialign" >> $LOG
    python3 /Users/lotti/CIAlign/benchmarking/QuanTest2/quantest2.py ct518@cam.ac.uk $aligned $SS >> $LOG
    exit
    cd CIAlign_run
done
