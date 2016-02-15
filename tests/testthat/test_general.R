

library(hutan)
library(ape)
context("general")

test_tree = read.tree(text="((Prayidae_D27SS7@1152004:0.310638840467,(Hippopodius_hippopus@867239:0.242845544461,(Craseoa_lathetica@180617:1.21823973672E-6,Prayidae_D27D2@922608:0.013252166117)n38047893-N-100:0.143340517912)n38047889-N-83:0.0697164264107)n38047887-Y-100:0.155909099947,(Kephyes_ovata@1402412:0.0624504448762,Chuniphyes_multidentata@278022:0.247764507969)n38047885-N-100:0.155909099947)n38048123-N;")

test_that("tree is decomposed into orthologs based on notung node names", {
	expect_equal( length( decompose_orthologs(test_tree) ), 3 )
})


decompose_orthologs