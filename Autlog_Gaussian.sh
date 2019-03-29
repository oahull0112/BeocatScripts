#!/bin/bash

# Take gaussian job submissions and log them automatically in the file "Glog"
# Stores the path to the submission and the submission ID number
# pair with the bashrc alias: alias gsub='bash path/to/Autolog_Gaussian.sh'
# and the aliases: alias push='echo >> ~/Glog.txt && data >> ~/Glog.txt' to make a line of today's date
# and the alias: alias log='tail -15 ~/Glog.txt' to view the last 15 lines of the file --> to see recently submitted jobs.

DESCRIPTION="$1"

SUBLINE=$(sbatch sbatch.in)
echo $SUBLINE >> ltmp0
awk '{print $4}' ltmp0 > ltmp1  # has job id in file only
grep 'gdv' sbatch.in | tail -1 > ftmp0
awk '{print $2}' ftmp0 > ftmp1 # has file name in file only
pwd >> ptmp
#paste ltmp1 ftmp1 > tmp0
paste ltmp1 ftmp1 ptmp > tmp0
awk -v var="${DESCRIPTION}" '{print "    ", $1," : ", $2, " : ", $3, " : " ,var}' tmp0 >> /homes/oliviahull/Glog.txt
rm tmp0 ftmp0 ltmp1 ltmp0 ptmp ftmp1
