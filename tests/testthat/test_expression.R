e=lapply(janedoe$expression, Expression)

test_that("number of genes in DEseq2 object equals 397", {
	expect_equal(length(rownames(e[[1]]@DESeq2)), 398)
})

test_that("number of columns in object equals 14", {
	expect_equal(length(colnames(e[[1]]@DESeq2)), 14)
})

#edgeR
test_that("number of genes in edgeR input equal 397", {
	expect_equal(nrow(e[[1]]@edgeR$counts), 398)
})


