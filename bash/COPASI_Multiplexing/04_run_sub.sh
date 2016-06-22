#!/bin/bash
# runs one optimization task on the cluster
# This script is invited by 01_submitter, therefore do not modify the name of the script.

/home/hub_oudenaarden/avertesy/bin/COPASI/COPASI-4.16.104-Linux-64bit/bin/CopasiSE $file -s $file -c /home/hub_oudenaarden/avertesy/bin/COPASI/COPASI-4.16.104-Linux-64bit/bin
	echo "submitting CopasiSE -s $file $file"
