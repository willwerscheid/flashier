set.seed(666)

gtex <- readRDS(gzcon(url("https://github.com/stephenslab/gtexresults/blob/master/data/MatrixEQTLSumStats.Portable.Z.rds?raw=TRUE")))
gtex <- gtex$strong.z
row.sample <- sample(1:nrow(gtex), 1000)
gtex <- gtex[row.sample, ]

devtools::use_data(gtex)
