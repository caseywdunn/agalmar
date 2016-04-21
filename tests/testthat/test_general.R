

library(hutan)
library(ape)
context("general")

# Example nhx tree from ggtree, https://github.com/GuangchuangYu/ggtree/blob/master/inst/extdata/ADH.nhx
test_nhx_text = "(((ADH2:0.1[&&NHX:S=human], ADH1:0.11[&&NHX:S=human]):0.05[&&NHX:S=primates:D=Y:B=100], ADHY:0.1[&&NHX:S=nematode],ADHX:0.12[&&NHX:S=insect]):0.1[&&NHX:S=metazoa:D=N], (ADH4:0.09[&&NHX:S=yeast],ADH3:0.13[&&NHX:S=yeast], ADH2:0.12[&&NHX:S=yeast],ADH1:0.11[&&NHX:S=yeast]):0.1 [&&NHX:S=Fungi])[&&NHX:D=N];"

test_that("can parse example ggtree nhx tree string", {
	nhx = parse_gene_tree( test_nhx_text )
	ntips = length(nhx@phylo$tip.label)
	expect_equal( ntips , 8 )
})




test_tree = read.tree(text="((Prayidae_D27SS7@1152004:0.310638840467,(Hippopodius_hippopus@867239:0.242845544461,(Craseoa_lathetica@180617:1.21823973672E-6,Prayidae_D27D2@922608:0.013252166117)n38047893-N-100:0.143340517912)n38047889-N-83:0.0697164264107)n38047887-Y-100:0.155909099947,(Kephyes_ovata@1402412:0.0624504448762,Chuniphyes_multidentata@278022:0.247764507969)n38047885-N-100:0.155909099947)n38048123-N;")

#test_that("tree is decomposed into orthologs based on notung node names", {
#	expect_equal( length( decompose_orthologs(test_tree) ), 3 )
#})


decompose_orthologs