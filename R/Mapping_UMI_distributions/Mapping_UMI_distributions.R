

setup_MarkdownReports(OutDir = "~/Github_repos/TheCorvinas/R/Mapping_UMI_distributions/")


UMI_distr_in_mapped_reads = 	read.simple.tsv("~/Github_repos/TheCorvinas/R/Mapping_UMI_distributions/UMI_distr_in_mapped_reads.tsv")
UMI_distr_in_unmapped_reads = 	read.simple.tsv("~/Github_repos/TheCorvinas/R/Mapping_UMI_distributions/UMI_distr_in_unmapped_reads.tsv")

AllMapped_Reads =	sum(UMI_distr_in_mapped_reads)
AllUnmapped_Reads =	sum(UMI_distr_in_unmapped_reads)

ExpectedOccurence = round(AllMapped_Reads* 1/256)
UMI_counts_in_Mapped_Reads = as.named.vector(UMI_distr_in_mapped_reads)
wbarplot(UMI_counts_in_Mapped_Reads, hline = ExpectedOccurence, ylab ="UMI observed", w=20)

UMI_enrichment_Mapped_Reads =	iround(as.named.vector(UMI_distr_in_mapped_reads/ (AllMapped_Reads* 1/256)))
UMI_enrichment_Unmapped_Reads =	iround(as.named.vector(UMI_distr_in_unmapped_reads/ (AllUnmapped_Reads* 1/256)))


wbarplot(UMI_enrichment_Mapped_Reads, hline = 1, ylab ="UMI observed/UMI expected", w=20)


order = sort(intersect(names(UMI_enrichment_Mapped_Reads),names(UMI_enrichment_Unmapped_Reads) ))


UMI_correlations_fold_above_random = cbind( "Mapped_Reads" = UMI_enrichment_Mapped_Reads[order],
                          "Unmapped_Reads" = UMI_enrichment_Unmapped_Reads[order])

wplot(UMI_correlations_fold_above_random)

fc_UMIwithAtLeast3GsandMaybe1T = (sum(UMI_distr_in_mapped_reads[c('GGGG', 'GGGT', 'TGGG', 'GGTG', 'GTGG'),]) / AllMapped_Reads)
"9.71 %"
percentage_formatter(fc_UMIwithAtLeast3GsandMaybe1T -(5*1/256))

