######################################################################
# Merge.CountTables.R
######################################################################
# source ('~/GitHub/TheCorvinas/R/Merge.CountTables.R')
# rm(list=ls(all.names = TRUE)); try(dev.off(), silent = T)

# Combinations.all.vs.all.R  ------------------------


listFiles =list(
  i2 = p0(inDir, 'INT02/Abel-INT2-3Aligned.toTranscriptome.out.coutt.csv.gz'),
  i7 = p0(inDir, 'INT07/Abel-int007Aligned.toTranscriptome.out.coutt.csv.gz'),
  i4 = p0(inDir, 'INT04_6_8_9/int4.TranscriptCounts.tsv'),
  i6 = p0(inDir, 'INT04_6_8_9/int6-9.TranscriptCounts.tsv'),
  i8 = p0(inDir, 'INT04_6_8_9/int8.TranscriptCounts.tsv')
)

LS_tables = lapply(listFiles, read.simple.tsv)

COMBINED = merge_dfs_by_rn(LS_tables)


