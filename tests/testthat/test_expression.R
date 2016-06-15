e=lapply(janedoe$expression, Expression)

deseq=e[[1]]@DESeq2

deseq

#edgeR
nrow(e[[1]]@edgeR$counts)


