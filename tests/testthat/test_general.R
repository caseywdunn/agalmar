

library(hutan)
library(ape)
context("nhx_manipulation")

test_newick_text = "((Prayidae_D27SS7@1152004:0.310638840467,(Hippopodius_hippopus@867239:0.242845544461,(Craseoa_lathetica@180617:1.21823973672E-6,Prayidae_D27D2@922608:0.013252166117)n38047893-N-100:0.143340517912)n38047889-N-83:0.0697164264107)n38047887-Y-100:0.155909099947,(Kephyes_ovata@1402412:0.0624504448762,Chuniphyes_multidentata@278022:0.247764507969)n38047885-N-100:0.155909099947)n38048123-N;"

# Example nhx tree from ggtree, https://github.com/GuangchuangYu/ggtree/blob/master/inst/extdata/ADH.nhx
test_nhx_text = "(((ADH2:0.1[&&NHX:S=human], ADH1:0.11[&&NHX:S=human]):0.05[&&NHX:S=primates:D=Y:B=100], ADHY:0.1[&&NHX:S=nematode],ADHX:0.12[&&NHX:S=insect]):0.1[&&NHX:S=metazoa:D=N], (ADH4:0.09[&&NHX:S=yeast],ADH3:0.13[&&NHX:S=yeast], ADH2:0.12[&&NHX:S=yeast],ADH1:0.11[&&NHX:S=yeast]):0.1 [&&NHX:S=Fungi])[&&NHX:D=N];"

test_phyldog_nhx_text = "(((Prayidae_D27SS7@2825365:0.0682841[&&NHX:Ev=S:S=58:ND=0],(Kephyes_ovata@2606431:0.0193941[&&NHX:Ev=S:S=69:ND=1],Chuniphyes_multidentata@1277217:0.0121378[&&NHX:Ev=S:S=70:ND=2]):0.0217782[&&NHX:Ev=S:S=60:ND=3]):0.0607598[&&NHX:Ev=S:S=36:ND=4],((Apolemia_sp_@1353964:0.11832[&&NHX:Ev=S:S=31:ND=9],(((Bargmannia_amoena@263997:0.0144549[&&NHX:Ev=S:S=37:ND=10],Bargmannia_elongata@946788:0.0149723[&&NHX:Ev=S:S=38:ND=11]):0.0925388[&&NHX:Ev=S:S=33:ND=12],Physonect_sp_@2066767:0.077429[&&NHX:Ev=S:S=61:ND=13]):0.0274637[&&NHX:Ev=S:S=24:ND=14],(Stephalia_dilata@2960089:0.0761163[&&NHX:Ev=S:S=52:ND=15],((Frillagalma_vityazi@1155031:0.0906068[&&NHX:Ev=S:S=53:ND=16],Resomia_ornicephala@3111757:1e-06[&&NHX:Ev=S:S=54:ND=17]):1e-06[&&NHX:Ev=S:S=45:ND=18],((Lychnagalma_utricularia@2253871:0.120851[&&NHX:Ev=S:S=65:ND=19],Nanomia_bijuga@717864:0.133939[&&NHX:Ev=S:S=71:ND=20]):1e-06[&&NHX:Ev=S:S=56:ND=21],Cordagalma_sp_@1525873:0.0693814[&&NHX:Ev=S:S=64:ND=22]):1e-06[&&NHX:Ev=S:S=46:ND=23]):0.0333823[&&NHX:Ev=S:S=40:ND=24]):1e-06[&&NHX:Ev=S:S=35:ND=25]):0.0431861[&&NHX:Ev=D:S=24:ND=26]):1e-06[&&NHX:Ev=S:S=19:ND=27],Rhizophysa_filiformis@3073669:0.22283[&&NHX:Ev=S:S=26:ND=28]):0.0292362[&&NHX:Ev=S:S=17:ND=29]):0.185603[&&NHX:Ev=D:S=17:ND=8],(Hydra_magnipapillata@52244:0.0621782[&&NHX:Ev=S:S=16:ND=5],Ectopleura_larynx@3556167:0.332505[&&NHX:Ev=S:S=15:ND=6]):0.185603[&&NHX:Ev=S:S=12:ND=7])[&&NHX:Ev=S:S=9:ND=30];"

test_phyldog_nhx_text_simple = '(((Prayidae_D27SS7@2825365,(Kephyes_ovata@2606431,Chuniphyes_multidentata@1277217)Ev-S_S-60_ND-3)Ev-S_S-36_ND-4,((Apolemia_sp_@1353964,(((Bargmannia_amoena@263997,Bargmannia_elongata@946788)Ev-S_S-33_ND-12,Physonect_sp_@2066767)Ev-S_S-24_ND-14,(Stephalia_dilata@2960089,((Frillagalma_vityazi@1155031,Resomia_ornicephala@3111757)Ev-S_S-45_ND-18,((Lychnagalma_utricularia@2253871,Nanomia_bijuga@717864)Ev-S_S-56_ND-21,Cordagalma_sp_@1525873)Ev-S_S-46_ND-23)Ev-S_S-40_ND-24)Ev-S_S-35_ND-25)Ev-D_S-24_ND-26)Ev-S_S-19_ND-27,Rhizophysa_filiformis@3073669)Ev-S_S-17_ND-29)Ev-D_S-17_ND-8,(Hydra_magnipapillata@52244,Ectopleura_larynx@3556167)Ev-S_S-12_ND-7)Ev-S_S-9_ND-30;'


test_notung_nhx_text = "((((Rhizophysa_filiformis@2564549:0.09666991738603078[&&NHX:S=Rhizophysa_filiformis],((Marrus_claudanielis@2027078:0.03368582974818837[&&NHX:S=Marrus_claudanielis],((Erenna_richardi@1434201:0.014306889954561298[&&NHX:S=Erenna_richardi],Marrus_claudanielis@2027079:0.010842363778569869[&&NHX:S=Marrus_claudanielis])n5940011:0.01779384958849464[&&NHX:S=n57:D=N],(((Agalma_elegans@88626:0.05872379503260147[&&NHX:S=Agalma_elegans],Lychnagalma_utricularia@1828459:0.04211137470826968[&&NHX:S=Lychnagalma_utricularia])n5940018:0.02375590664436535[&&NHX:S=n47:D=N],(((Bargmannia_amoena@3459111:0.19058396964770352[&&NHX:S=Bargmannia_amoena],Bargmannia_elongata@469437:1.00000050002909E-6[&&NHX:S=Bargmannia_elongata])n5939974:0.11560220708003867[&&NHX:S=n22:D=N],Cordagalma_sp_@1115328:0.04829417133033771[&&NHX:S=Cordagalma_sp_])n5939976:0.011316847557531757[&&NHX:S=n62:D=N],Forskalia_asymmetrica@1220430:0.01667566952752948[&&NHX:S=Forskalia_asymmetrica])n5939978:0.0063213422810751655[&&NHX:S=n62:D=Y])n5940014:0.017792661031819083[&&NHX:S=n62:D=Y],(Resomia_ornicephala@2657185:0.004262563771468986[&&NHX:S=Resomia_ornicephala],Frillagalma_vityazi@663744:0.028441637105547157[&&NHX:S=Frillagalma_vityazi])n5939981:0.006136291467151878[&&NHX:S=n51:D=N])n5940013:0.013546839136761205[&&NHX:S=n62:D=Y])n5940012:0.011839606018978143[&&NHX:S=n62:D=Y])n5940008:0.013840645450221475[&&NHX:S=n62:D=Y],(((Chelophyes_appendiculata@1615707:0.007647023552225329[&&NHX:S=Chelophyes_appendiculata],Clytia_hemisphaerica@756642:0.643907456299178[&&NHX:S=Clytia_hemisphaerica])n5939984:0.08603691877960613[&&NHX:S=n67:D=N],(Chuniphyes_multidentata@930929:0.01248550133310033[&&NHX:S=Chuniphyes_multidentata],Kephyes_ovata@1966030:0.014671165587181996[&&NHX:S=Kephyes_ovata])n5939987:0.013285803501636162[&&NHX:S=n27:D=N])n5939988:0.008000411801689693[&&NHX:S=n67:D=Y],(((Hippopodius_hippopus@1084434:0.0505718831943577[&&NHX:S=Hippopodius_hippopus],Prayidae_D27D2@2878798:0.00905875758406546[&&NHX:S=Prayidae_D27D2])n5939991:0.021772123626769023[&&NHX:S=n38:D=N],Prayidae_D27SS7@2181711:0.029009000260863272[&&NHX:S=Prayidae_D27SS7])n5939993:1.00000050002909E-6[&&NHX:S=n38:D=Y],Prayidae_D27D2@2878801:1.00000050002909E-6[&&NHX:S=Prayidae_D27D2])n5939995:0.00916688375355408[&&NHX:S=n38:D=Y])n5939996:0.05191099091093772[&&NHX:S=n67:D=Y])n5940006:0.03953811931719265[&&NHX:S=n67:D=Y])n5940005:0.10134081070615458[&&NHX:S=n67:D=Y],(Podocoryna_carnea@3033951:0.11270255504816476[&&NHX:S=Podocoryna_carnea],Hydractinia_symbiolongicarpus@1679508:0.030168043235021993[&&NHX:S=Hydractinia_symbiolongicarpus])n5939999:0.17223048099362362[&&NHX:S=n11:D=N])n5940003:0.16233679521228994[&&NHX:S=n67:D=Y],Hydra_magnipapillata@801936:0.585696573276294[&&NHX:S=Hydra_magnipapillata])n5940002:0.4403044529817829[&&NHX:S=n68:D=N],Aegina_citrea@825314:0.4403044529817829[&&NHX:S=Aegina_citrea])n5942419[&&NHX:S=n70:D=N];"

test_that("can parse example ggtree nhx tree string", {
	nhx = parse_gene_tree( test_nhx_text )
	ntips = length(nhx@phylo$tip.label)
	expect_equal( ntips , 8 )
})

test_that("can parse phyldog nhx tree", {
	nhx = parse_gene_tree( test_phyldog_nhx_text )

	# Test that NHX tags correctly parsed into nhx object
	tags = nhx@nhx_tags
	tags$node = as.numeric(tags$node)
	tags = tags[ !is.na(tags$node), ]
	tags = tags[ tags$node > length(nhx@phylo$tip.label), ] # Consider internal nodes o
	tags = tags[order(tags$node),]

	# Compare a node via mrca
	n = ape::mrca(nhx@phylo)["Hydra_magnipapillata@52244", "Ectopleura_larynx@3556167"]
	expect_equal( as.numeric(tags[tags$node==n,]$S), 12 )

	# Compare all nodes, assuming in same order as independently 
	# parsed phylo object of a simplified tree text that can be
	# read by ape
	phylo=read.tree(text=test_phyldog_nhx_text_simple)
	phylo_S=unlist(lapply(strsplit(phylo$node.label, "_"), function(x) x[[2]])) # Get the S field
	phylo_S=unlist(lapply(strsplit(phylo_S, "-"), function(x) x[[2]])) # Get the value
	phylo_S=as.numeric(phylo_S)
	expect_equal( phylo_S, as.numeric(tags$S) )

	# Test that the NHX annotations were correctly added to the phylo node labels
	# This is a feature of how agalmar parses nhx trees, not how ggplot 
	# parses them
	get_nhx_by_tips = function(phy, tip1, tip2){
		n = ape::mrca(phy)[tip1, tip2]
		label = phy$node.label[ n - length(phy$tip.label) ]
		fields = nhx_label_to_list(label)
		return(fields)
	}

	# Test a single node
	p=nhx@phylo
	fields = get_nhx_by_tips(p, "Hydra_magnipapillata@52244", "Ectopleura_larynx@3556167")
	expect_equal( as.numeric(fields$S), 12 )

	# Now check them all, assuming conserved order
	p_nodes = sapply(p$node.label, function(x) as.numeric(nhx_label_to_list(x)$S))
	names(p_nodes) = NULL
	expect_equal( phylo_S, p_nodes )
})


test_that("tree is decomposed into orthologs based on notung nhx annotations", {
	nhx = parse_gene_tree( test_notung_nhx_text )
	decomposed = decompose_orthologs(nhx)
	expect_equal( length(decomposed), 14 ) # haven't checked that it should be 14
	# ggtree(nhx) + geom_tiplab() + geom_point(aes(color=D), size=5, alpha=.5) + theme(legend.position="right")


	# Check that all tips in the original tree are in the final trees
	decomposed_tips = unlist( lapply( decomposed, function(x) x$tip.label ) )
	expect_true( setequal( decomposed_tips, nhx@phylo$tip.label ) )
})


test_that("tree is decomposed into orthologs based on phyldog nhx annotations", {
	nhx = parse_gene_tree( test_phyldog_nhx_text )
	decomposed = decompose_orthologs(nhx)
	expect_equal( length(decomposed), 5 )
	# ggtree(nhx) + geom_tiplab() + geom_point(aes(color=D), size=5, alpha=.5) + theme(legend.position="right")

})

test_that("can summarize edges", {
	nhx = parse_gene_tree( test_phyldog_nhx_text )
	edges = summarize_edges(nhx)
	expect_equal( nrow(edges), 30 )
})


test_that("can summarize nodes", {
	nhx = parse_gene_tree( test_phyldog_nhx_text )
	nodes = summarize_nodes(nhx)
	expect_equal( nrow(nodes), 31 )
})

test_that("can drop nhx tips", {
	nhx = parse_gene_tree( test_phyldog_nhx_text )
	to_drop = c("Physonect_sp_@2066767", "Lychnagalma_utricularia@2253871", "Kephyes_ovata@2606431")
	
	nhx_reduced = drop.tip.nhx(nhx, to_drop, test=FALSE)
	expect_equal( length(nhx_reduced@phylo$tip.label), 13 )

	nhx_reduced = drop.tip.nhx(nhx, to_drop, test=TRUE)
	expect_equal( length(nhx_reduced@phylo$tip.label), 13 )
})
drop.tip.nhx

# test_tree = read.tree( text=test_newick_text )



