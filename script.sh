#!/usr/bin/env bash

#cmd = "-R -x-X"

MDIR="/mnt/data0/carlom/samples_project/"

subs = `ls $MDIR`
mutect = 'mutect'

for dir in $subs; do
mv $MDIR/$dir/$mutect/* $main_dir/$dir/
done