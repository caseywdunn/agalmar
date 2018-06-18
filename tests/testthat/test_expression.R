e=lapply(janedoe$expression, Expression)

context("DESeq2")

de2 = create_DESeq2( e[[1]], design = ~individual+treatment )

test_that("number of genes in DEseq2 object equals 422", {
	expect_equal(length(rownames(de2)), 422)
})

test_that("number of columns in object equals 14", {
	expect_equal(length(colnames(de2)), 14)
})

context("EdgeR")
test_that("number of genes in edgeR input equal 423", {
	expect_equal(nrow(e[[1]]@edgeR$counts), 423)
})
#this number is different than the DEseq object, as one of the libraries has rowSums <1 

context("summarize_libraries")
test_that("can construct a library summary", {
	expect_equal(nrow(summarize_libraries(e[[1]])), 14)
})