# Make S3 style classes available to S4 code
setOldClass(c("DGEList"))
################################################################################
# Classes

#' Summarize experiment.
#' 
#' @param object An Expression object
#' @return A data frame that summarizes each sample in the experiment
#' @export
setGeneric (
	name = "summary_frame", 
	def = function( object ) 
		{ standardGeneric("summary_frame") }
)

#' Get species
#' 
#' @param object An Expression object
#' @return Character string with species
#' @export
setGeneric (
	name = "species", 
	def = function( object ) 
		{ standardGeneric("species") }
)

#' An S4 class to represent gene expression data for multiple samples
#' for a given species. Assumes that all expession data are derived 
#' from mapping to the same reference. Applies to data for g genes
#' across s samples.
#'
#' @slot species  The species
#' @slot dge  an edgeR DGEList object holding data for all samples. Expression
#' matrix of dimension g,s.
#' @slot lengths  Summary of length, in bp, of transctipts for each gene. Length g.
#' @slot individuals  Factors indicating which individual each sample is from. Length s.
#' @slot treatments  Factors indicating which treatment applies to each sample. Length s.
#' @slot id  Factors indicating the unique id of each sample. Length s.
#' @slot samples  Factors indicating the unique library id of each sample. Length s.
#' @slot sample_prep  Sample prep strategy. Length s.
#' @slot genome_type  Character indicating genome type, eg nuclear. See agalma documentation. Length g.
#' @slot molecule_type  Character indicating encoded molecule type, eg protein. See agalma documentation. Length g.
#' @slot blast_hit  Blast hit. Length g.
#' @slot rRNA  Fraction of reads in sample that are rRNA. Length s.
#' @slot protein  Fraction of reads in sample that are protein coding. Length s.
setClass(
	Class = "Expression",
	representation = representation (
		species = "character", 
		edgeR = "DGEList",
		DESeq2 = "DESeqDataSet",
		lengths = "matrix",
		individuals = "vector",
		treatments = "vector",
		id = "vector",
		samples = "vector",
		sample_prep = "vector",
		genome_type = "vector",
		molecule_type = "vector",
		blast_hit =  "vector",
		rRNA = "vector",
		protein = "vector"
	)
)


#' Construct an Expression object from list of experiment data and metadata 
#' provided by Agalma.
#' 
#' @param data_list A list containing the expression data
#' @return An Expression object
#' @export
Expression <- function( data_list ) {
	object <- new("Expression")
	
	object@species <- data_list$species


	
	# Parse column annotations
	object@individuals <- as.factor( data_list$individual )
	object@treatments  <- as.factor( data_list$treatment )
	object@id  <- as.factor( data_list$id )
	object@samples <- as.factor( data_list$library_id )
	object@sample_prep <- data_list$sample_prep
	
	# Simplify sample prep names
	object@sample_prep[grep("Illumina TruSeq", object@sample_prep)] <- "Illumina TruSeq"
	object@sample_prep[grep("Illumina mRNA-Seq", object@sample_prep)] <- "Illumina mRNA-Seq"
	object@sample_prep[grep("NEBNext", object@sample_prep)] <- "NEBNext"
	
	# Parse row annotations
	object@genome_type <- as.factor( data_list$genome_type )
	object@molecule_type <- as.factor( data_list$molecule_type )

	if ( exists( 'blast_hit', where=data_list ) ){
		# old name, retained for compatability
		object@blast_hit <- data_list$blast_hit
	} else {
		object@blast_hit <- data_list$blast_title
	}
	
	# Parse counts matrix
	x <- data_list$count
	if ( is.null( dim(x) ) ){
		# if there is a single samply, then x is a vector rather than an array, 
		# which messes things up for row sampling later
		dim( x ) <- c( length(x), 1 )
	}
	rownames( x ) <- data_list$gene
	colnames( x ) <- object@samples
	
	# Parse the lengths, if present
	if ( exists( 'length', where=data_list ) ){
		object@lengths <- data_list$length
	} else{
		empty_lengths <- rep( NA, length(x) )
		dim( empty_lengths ) <- dim( x )
		object@lengths <- empty_lengths
	}


	# Calculate total counts, including rRNA
	totals <- colSums( x )
	
	# Quantify rRNA and identify non ribosomal RNA
	rrna <- object@molecule_type %in% c('L','S')

	if ( sum(rrna) == 0 ){
		object@rRNA <- rep( 0, ncol(x) )
	} else if ( sum(rrna) == 1 ){
		object@rRNA <- x[rrna,]/totals
	} else {
		object@rRNA <- colSums(x[rrna,])/totals
	}
	
	# Identify protein coding genes
	protein_coding <- ( object@molecule_type == 'P' )
	object@protein <- colSums(x[protein_coding,])/totals

	# Identify rows that have at last 2 libraries with count greater than 0
	passes_sampling_criterion <- rowSums(x > 0) > 2

	# Exclude plastid genomes
	genome_keep <- (object@genome_type != 'P') & (object@genome_type != 'M')
	
	# Create a vector of rows to keep
	keep <- protein_coding & genome_keep & passes_sampling_criterion
	
	# Subsample matrix and row annotations
	x <- x[keep,]
	object@lengths <- object@lengths[keep,]
	object@genome_type <- object@genome_type[keep]
	object@molecule_type <- object@molecule_type[keep]
	object@blast_hit <- object@blast_hit[keep]
	
	# Prepare EdgeR DGE object
	object@edgeR <- edgeR::DGEList( counts=x, group=object@treatments )
	object@edgeR <- edgeR::calcNormFactors( object@edgeR )
	
	# Prepare DESeq2 DGE object
	colData=data.frame(treatment=data_list$treatment,individual=data_list$individual, row.names=data_list$library_id)
	
	dds <- DESeqDataSetFromMatrix(countData = x,
	                              colData = colData,
	                              design = ~individual+treatment)
	
	object@DESeq2 <- dds[ rowSums(counts(dds)) > 1, ]
	
	object
}


setMethod("summary_frame", signature(object = "Expression"),
	function(object) {
		sample_summary <- data.frame(Species=rep(object@species, length(object@samples)), Individuals=object@individuals, Treatments=object@treatments, Samples=object@samples, Preparation=object@sample_prep, rRNA=object@rRNA, Protein=object@protein, Reads=colSums(object@dge$counts), Run=as.factor(sapply(object@id, get_run)), Lane=as.factor(sapply(object@id, get_lane)))
		
		return( sample_summary )
	}
)


setMethod("species", signature(object = "Expression"),
	function(object) {
				
		return( object@species )
	}
)


################################################################################
# Functions

#' Get run from Illumina fastq sequence header
#' 
#' @param header character strings
#' @return run character strings
#' @export
get_run <- function( header ) {
	header <- as.character(header)
	run <- ''
	if (substring(header,1,1) == 'H') {
		fields <- strsplit(header, "-")
		fields <- unlist(fields)
		run <- paste(fields[1],fields[2],fields[3], sep = "-")
	} else {
	  	fields <- strsplit(header, "-")
		fields <- unlist(fields)
		run <- paste(fields[1],fields[2], sep = "-")
    }
    
	return( run )
}



#' Get lane from Illumina fastq sequence header
#' 
#' @param header character strings
#' @return lane character strings
#' @export
get_lane <- function( header ) {
	header <- as.character(header)

	fields <- strsplit(header, "-")
	fields <- unlist(fields)
	lane <- paste(fields[length(fields)-1])

	return( lane )
}


#' Parses newick text to a gene tree
#' 
#' @param tree_text Text representation of a newick tree in newick format.
#' @return phy The tree, as an ape phylo object
#' @export
parse_gene_tree <- function( tree_text ){

	phy <- ape::read.tree( text= tree_text )
	return( phy )
}


#' Parse support from node name
#' 
#' @param node_name A character string of format 'n38047889-N-83', where 'n38047889-N' is 
#' notung annotation and '83' is node support.
#' @return Numeric indicating node support
#' @export
node_support <- function( node_name ) {

	ns <- NA
	fields <- strsplit(node_name, "-")[[1]]
	if (length( fields ) == 3 ){
		ns = as.numeric( fields[3] )
	}

	return( ns )
}


#' Parse species and sequence id from a phy
#' 
#' @param phy The tree, as an ape phylo object. Tip names must have `species@@id` format.
#' @return Dataframe with one row per tip, a species column, and an id column
#' @export
get_tip_info <- function( phy ) {
	
	tips <- phy$tip.label
	tip_info <- as.data.frame(matrix(unlist(strsplit(tips, split='@')), nrow=length(tips), byrow=T))
	names(tip_info) <- c("species", "id")
	tip_info$id <- as.numeric( as.character( tip_info$id ) ) # Change from factors to numeric

	# species names have spaces in expression data, remove underscores to make them consistent
	tip_info[,1] <- sub('_', ' ', tip_info[,1])
	
	return( tip_info )

}


#' Returns true if all species are in the tree tips
#' 
#' @param phy The tree, as an ape phylo object
#' @param species Vector of species names
#' @return logical
#' @export
has_species <- function( phy, species ) {

	psp <- get_tip_info( phy )[,1]
	
	present = species %in% psp
	
	return( all( present ) )  
	
}


#' Takes a DGEList object, and returns a matrix of normalized counts
#' 
#' @param dge edgeR DGEList object, to which normalizations have been applied
#' @return matrix of normalized counts
#' @export
apply_normalizations <- function(dge){
	
	# Calculate the normalization multipliers, scaled to reads per million
	norm <- 1e6 / (dge$samples[,"lib.size"] * dge$samples[,"norm.factors"])
	
	# Apply them to the counts
	norm_counts <- t(t(dge$counts)*norm)
	
	return(norm_counts)
}


#' Plots a matrix
#' 
#' @param m The matrix
#' @param ... additional arguments for image
#' @export
plot_matrix <- function(m, ... ) {

	nr <- nrow(m)
	nc <- ncol(m)
	image(1:nc, 1:nr, t(m[nr:1, ]), axes=F,xlab="", ylab="", ... )
}


#' Decomposes a gene tree into a list of subtrees that have no duplication events. Assumes notung
#' style node annotations.
#' 
#' @param phy The tree, as an ape phylo object
#' @return The subtrees as a list of phylo object
#' @export
decompose_orthologs <- function( phy ){

	# Identify duplicated nodes 
	duplications = grep( "-Y", phy$node.label ) + length( phy$tip.label )

	# Get their immediate descendants, which define the clades we want to excise
	to_prune = phy$edge[,2][ phy$edge[,1] %in% duplications ]

	subtrees = hutan::decompose_tree( phy, to_prune )

	return( subtrees )

}