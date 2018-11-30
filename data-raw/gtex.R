set.seed(666)

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
gtex <- gtex$strong.z
row.sample <- sample(1:nrow(gtex), 1000)
gtex <- gtex[row.sample, ]

devtools::use_data(gtex)


missing.tissues <- c(7, 8, 19, 20, 24, 25, 31, 34, 37)
gtex.colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", sep = '\t', comment.char = '')[-missing.tissues, 2]
gtex.colors <- as.character(gtex.colors)

devtools::use_data(gtex.colors)
