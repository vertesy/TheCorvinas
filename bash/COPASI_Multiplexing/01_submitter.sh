#! /bin/bash
# COPASI multiplex parameter estimation wrapper scripts for unix servers with a queue system
# This script prepares the files/folders and submit the jobs on the cluster

# Replace the following: ---------------------------------
	# - Iteration count
	# - Folder name HAS TO END WITH '/'
	# - Modelfile's name
# In copasi you have to:
	# - Make the model executable
	# - Define a single-line output
	# - Check "Randomize Initial values"

NrOfFits=50

# workdir="/hpc/hub_oudenaarden/Abel/ppw/Q5/"
workdir="/hpc/hub_oudenaarden/Abel/ppw/Q6/"
file="Q5.cps"

ScriptLoc="/home/hub_oudenaarden/avertesy/bin/COPASI/Scripts/04_run_sub.sh"

# get to main directory
cd $workdir
# if echo $workdir | grep \/$ ; then echo ""; else echo "WORKDIR HAS TO END WITH '/"; exit 1; fi
if echo $workdir | grep \/$ ; then echo ""; else echo "WORKDIR HAS TO END WITH '/"; fi

# copy model and data files into subdirectories and submit the job
for ((i="1"; i<=NrOfFits; i++)) ; do
	mkdir -p $i
	cp $file $i
	cd $i
	qsub -N c$i'-'$file -cwd -V -l h_rt=11:55:00 -l h_vmem=3G -M youremail@hubrecht.eu -m a -v file=$workdir$i/$file $ScriptLoc -d $workdir/$i/
	cd ..
done
echo -N c$i'-'$file -cwd -V -l h_rt=11:55:00 -l h_vmem=3G -M youremail@hubrecht.eu -m a -v file=$workdir$i/$file $ScriptLoc -d $workdir/$i/

qstat
