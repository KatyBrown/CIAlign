#!/bin/bash

FILES=$1'*'
DESTDIR_CIALIGN=$2

for f in $FILES
  do
    NAME=${f##*/}      # remove part before last slash
    NAME=${NAME%.*}    # from the new var remove the part after the last period
    SS="../QuanTest2/SS/"$NAME".ss"
    aligned=$DESTDIR_CIALIGN"MSAs/"$NAME"_muscle.fasta"
    OUTFILE_CIAlign=$DESTDIR_CIALIGN"CIAlign_MSAs/"$NAME
    REPLACED_MSA=$OUTFILE_CIAlign"_replaced.fasta"
    LOG=$DESTDIR_CIALIGN$NAME"_log.txt"

    echo $f
    echo $NAME >> $LOG
    muscle -in $f -out $aligned
    CIAlign --infile $aligned --outfile $OUTFILE_CIAlign --clean --remove_divergent_minperc 0.15
    python3 adjust_msa_to_quantest2.py $aligned $OUTFILE_CIAlign"_removed.txt" $REPLACED_MSA
    echo "with cialign" >> $LOG
    python3 quantest2_cialign.py mail@mail.com $OUTFILE_CIAlign"_cleaned.fasta" $SS $REPLACED_MSA >> $LOG
    echo "without cialign" >> $LOG
    python3 quantest2.py mail@mail.com $aligned $SS >> $LOG
done
