#!/bin/bash

FILES=$1'*'
# DESTDIR_ALI="/home/charlotte/homfam/aligned/"
DESTDIR_CIALIGN="/home/charlotte/homfam/results/"

for f in $FILES
  do
    NAME=${f##*/}      # remove part before last slash
    NAME=${NAME%.*}    # from the new var remove the part after the last period
    OUTFILE_CIAlign=$DESTDIR_CIALIGN$NAME
    echo $f
    # python3 check_length.py $f

    python3 /home/charlotte/CIAlign/CIAlign/CIAlign.py --infile $f --outfile $OUTFILE_CIAlign"crop_ends" --crop_ends
    python3 /home/charlotte/CIAlign/CIAlign/CIAlign.py --infile $f --outfile $OUTFILE_CIAlign"remove_insertions" --remove_insertions
    python3 /home/charlotte/CIAlign/CIAlign/CIAlign.py --infile $f --outfile $OUTFILE_CIAlign"remove_divergent" --remove_divergent --remove_divergent_minperc 0.3
    python3 /home/charlotte/CIAlign/CIAlign/CIAlign.py --infile $f --outfile $OUTFILE_CIAlign"remove_short" --remove_short

done
