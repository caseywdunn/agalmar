e=lapply(janedoe$expression, Expression)

context("DESeq2")

de2 = create_DESeq2( e[[1]], design = ~individuals+treatments )

test_that("number of genes in DEseq2 object equals 397", {
	expect_equal(length(rownames(de2)), 398)
})

test_that("number of columns in object equals 14", {
	expect_equal(length(colnames(de2)), 14)
})

context("EdgeR")
test_that("number of genes in edgeR input equal 397", {
	expect_equal(nrow(e[[1]]@edgeR$counts), 398)
})


