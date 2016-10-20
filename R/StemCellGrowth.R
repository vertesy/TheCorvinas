######################################################################
# Stem cell growth rate
######################################################################

# Functions ------------------------
require(MarkdownReports)
# Setup ------------------------

# Go ------------------------
F2a_result = read.csv("/Users/abelvertesy/Downloads/F2a_result.csv")

MeanCellCount = rowMeans(F2a_result[, 2:4])

# apply(array, margin (1 for row , 2 for col), FUNCTION)
SD_CellCount = apply(X = F2a_result[, 2:4], MARGIN = 1, FUN = sd)

MeanGrowth =  cbind(  "Hours" = F2a_result$X,
                      "Cell Count" = MeanCellCount )

plot(MeanGrowth, ylim= c(0, 30), pch=18, 
     sub="Standard Deviation is depicted")
error_bar(x = F2a_result$X, y = MeanCellCount, upper = SD_CellCount)
wplot_save_this(plotname = "ESC Cell Growth")

