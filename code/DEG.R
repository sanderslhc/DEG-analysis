rm(list = ls())
options(stringsAsFactors = F)
#import count files
a <- read.table("D:/R-3.6.2/wrokspace/counts_DP.txt", header=T)
b <- read.table("D:/R-3.6.2/wrokspace/counts_NP.txt", header=T)
#constract count matrix
count=a[,c(1,7:9)]
names(count)=c("Gene","DP1","DP2","DP3")
count$NP1=b[,7]
count$NP2=b[,8]
count$NP3=b[,9]
row.names(count) <- count $Gene
count <- count[,-1]

X = c("DP1","DP2","DP3","NP1","NP2","NP3")
condition = c(1,1,1,2,2,2)
type = c(1,1,1,1,1,1)
coldata = data.frame(X, condition, type)
row.names(coldata) <- coldata$X
coldata <- coldata[,-1]

library("DESeq2")
dds <- DESeqDataSetFromMatrix(countData = count,
                              colData = coldata,
                              design = ~ condition)
dds <- DESeq(dds)
res <- results(dds)
res
summary(res)
write.csv(res, file = 'D:/R-3.6.2/results/DEG.csv')