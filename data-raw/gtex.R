set.seed(666)

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
gtex <- gtex$strong.z
row.sample <- sample(1:nrow(gtex), 1000)
gtex <- gtex[row.sample, ]

usethis::use_data(gtex, overwrite = TRUE)


missing.tissues <- c(7, 8, 19, 20, 24, 25, 31, 34, 37)
gtex_colors <- read.table("https://github.com/stephenslab/gtexresults/blob/master/data/GTExColors.txt?raw=TRUE", sep = '\t', comment.char = '')[-missing.tissues, 2]
gtex_colors <- as.character(gtex_colors)
names(gtex_colors) <- colnames(gtex)

usethis::use_data(gtex_colors, overwrite = TRUE)
