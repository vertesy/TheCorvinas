######################################################################
# Copy.files.sh
######################################################################




# scratch-cbe ------------------------------------------------------------------------------------------
# This scratch space is located on our very fast, all flash-based BeeGFS cluster filesystem. It is accessible under /scratch-cbe/user/<username> or via the $SCRATCHDIR environment variable.
# Data lifetime is limited to 20 days with automatic cleanup according to the modification date of individual files. Please make sure to move important data elsewhere on time! There are no snapshots or backup.

$SCRATCHDIR
mkdir CONscr/ CONscr/Connectomics.6.to.19 CONscr/Connectomics.6.to.19/Analysis/


scratch_dirCON=$SCRATCHDIR'/CONscr/Connectomics.6.to.19/Analysis/'
ls $scratch_dirCON


cd /groups/knoblich/raw.fastq.files.from.NGS/single.cell.RNA.seq/A.Vertesy/CONN/Connectomics.6.to.19/Analysis/
cp ./combined.obj.all.cells_CON_computed_2023.09.05_14.06.Rds.gz $scratch_dirCON &
cp ./ls.Seurat_CON_b4.GT.annot_2023.09.04_18.51.Rds.gz $scratch_dirCON &
ls $scratch_dirCON

$scratch_dir
