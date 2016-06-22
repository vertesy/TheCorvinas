#!/bin/sh
# This script merges single output files to one big table
# It assumes that Copasi outputs a single-line textual outputs, and merges them line by line, adding a number in the first column denoting which .cps file (folder) they correspond to.
# The header assumes that the first column is the objective value, followed by the corresponding parameters set within brackets that are separate columns.

source ~/.bash_profile;


workdir="/hpc/hub_oudenaarden/Abel/ppw/Q6/"
NrOfFits=50

# get to main directory
cd $workdir
	echo "working directory is $workdir"

# Calulate the numeber of columns and parameters based on output #1
NrColumns=$(head "$workdir/1/output.tab" | awk '{print NF}')
echo $NrColumns" Columns"
NrParams=`expr $NrColumns + 3`

# Add header, and copy output lines one by one. Did not work tab-separated, so we convert from space-separated to tab-separated at the end of this script.
echo -e "ID ObjValue x p1 p2 p3 p4 p5 p6 p7 p8 p9 p10 p11 p12 p13 p14 p15 p16 p17 p18 p19 p20 p21 p22 p23 p24 p25 p26 p27 p28 p29 p30 p31 p32 p33 p34 p35 p36 p37 p38 p39 p40 p41 p42 p43 p44 p45 p46 p47 p48 p49 p50" > 'FitResults.tsv'

for ((i="1"; i<=NrOfFits; i++)) ; do
	echo -e $i" "$(cat "$i/output.tab") >> 'FitResults.tsv'
done

# Cut out the 3rd (2nd in new file) and the last columns (brackets)
cut -f1-2,4-$NrColumns -d " " 'FitResults.tsv' > $workdir'_FitResults.tsv'

# Convert from space-separated to tab-separated values
sed "s/ \+/\t/g" _FitResults.tsv > FitResults.tsv
rm -rf "_FitResults.tsv"

# check results and download from the server.
head $workdir'FitResults.tsv'
download $workdir'FitResults.tsv'
