# Make S3 style classes available to S4 code

#' "DGEList" class
#'
#' @name DGEList-class
setOldClass( c( "DGEList" ) )

################################################################################
# Classes

#' Summarize expression libraries
#'
#' @param object An Expression object
#' @return A data frame that summarizes each library with expression data in the
#' Expression object
#' @export
setGeneric (
	name = "summarize_libraries",
	def = function( object )
		{ standardGeneric( "summarize_libraries" ) }
)


#' Summarize reference sequences
#'
#' @param object An Expression object
#' @return A data frame that summarizes the reference sequences in the
#' Expression object
#' @export
setGeneric (
	name = "summarize_reference",
	def = function( object )
		{ standardGeneric( "summarize_reference" ) }
)



#' Get species
#'
#' @param object An Expression object
#' @return Character string with species
#' @export
setGeneric (
	name = "species",
	def = function( object )
		{ standardGeneric( "species" ) }
)


#' Creates a DESeqDataSet object
#'
#' @param object An Expression object
#' @param design a formula that explains the project design
#' @return A DESeqDataSet object
#' @export
setGeneric (
	name = "create_DESeq2",
	def = function( object, design )
		{ standardGeneric( "create_DESeq2" ) }
)


#' An S4 class to represent gene expression data for multiple samples
#' for a given species. Assumes that all expession data are derived
#' from mapping to the same reference. Applies to data for g genes
#' across s samples (ie, sequenced libraries).
#'
#' The fields that apply to the s samples correspond to those
#'
#' @slot species  The species
#' @slot edgeR  An edgeR DGEList object holding data for all samples. Expression
#' matrix of dimension g,s.
#' @slot lengths  Summary of length, in bp, of transctipts for each gene. Length g.
#' @slot individual  Factors indicating which individual each sample is from. Length s.
#' @slot treatment  Factors indicating which treatment applies to each sample. Length s.
#' @slot id  Factors indicating the unique sequencing run id of each sample, eg HWI-ST625-75-D0PBDACXX-6-ATCACG. Length s.
#' @slot library_id  Factors indicating the unique library id of each sample, eg FEG365. Length s.
#' @slot sample_prep  Sample prep strategy. Length s.
#' @slot genome_type  Character indicating genome type, eg nuclear. See agalma documentation. Length g.
#' @slot molecule_type  Character indicating encoded molecule type, eg protein. See agalma documentation. Length g.
#' @slot blast_hit  Blast hit. Length g.
#' @slot rRNA  Fraction of reads in sample that are rRNA. Length s.
#' @slot protein  Fraction of reads in sample that are protein coding. Length s.
#' @slot x The counts matrix
#' @importClassesFrom DESeq2 DESeqDataSet
setClass(
	Class = "Expression",
	representation = representation(
		species = "character",
		edgeR = "DGEList",
		lengths = "vector",
		individual = "vector",
		treatment = "vector",
		id = "vector",
		library_id = "vector",
		sample_prep = "vector",
		genome_type = "vector",
		molecule_type = "vector",
		blast_hit =  "vector",
		rRNA = "vector",
		protein = "vector",
		x = "matrix",
		fpkm= "matrix",
		tpm="matrix"
	)
)


#' Construct an Expression object from list of experiment data and metadata
#' provided by Agalma.
#'
#' @param data_list A list containing the expression data
#' @return An Expression object
#' @export
Expression = function( data_list ) {
	object = methods::new( "Expression" )

	object@species = data_list$species_name



	# Parse column annotations
	object@individual = as.factor( data_list$individual )
	object@treatment  = as.factor( data_list$treatment )
	object@id  = as.factor( data_list$id )
	object@library_id = as.factor( data_list$library_id )
	object@sample_prep = data_list$sample_prep

	# Simplify sample prep names
	object@sample_prep[grep( "Illumina TruSeq", object@sample_prep )] = "Illumina TruSeq"
	object@sample_prep[grep( "Illumina mRNA-Seq", object@sample_prep )] = "Illumina mRNA-Seq"
	object@sample_prep[grep( "NEBNext", object@sample_prep )] = "NEBNext"

	# Parse row annotations
	object@genome_type = as.factor( data_list$genome_type )
	object@molecule_type = as.factor( data_list$molecule_type )

	if ( exists( 'blast_hit', where=data_list ) ){
		# old name, retained for compatability
		object@blast_hit = data_list$blast_hit
	} else {
		object@blast_hit = data_list$blast_title
	}

	# Parse fpkm values
	object@fpkm =data_list$fpkm
	rownames( object@fpkm ) = data_list$sequence_id
	colnames( object@fpkm ) = object@library_id

	# Parse tpm values
	object@tpm =data_list$tpm
	rownames( object@tpm ) = data_list$sequence_id
	colnames( object@tpm ) = object@library_id

	# Parse counts matrix
	object@x = data_list$count
	if ( is.null( dim( object@x ) ) ){
		# if there is a single sample, then x is a vector rather than an array,
		# which messes things up for row sampling later
		dim( object@x ) = c( length( object@x ), 1 )
	}
	rownames( object@x ) = data_list$sequence_id
	colnames( object@x ) = object@library_id

	# Parse the lengths, if present
	if ( exists( 'gene_length', where=data_list ) ){
		object@lengths = data_list$gene_length
	} else{
		empty_lengths = rep( NA, nrow( object@x ) )
		object@lengths = empty_lengths
	}


	# Calculate total counts, including rRNA
	totals = colSums( object@x )

	# Quantify rRNA and identify non ribosomal RNA
	rrna = object@molecule_type %in% c( 'L','S' )

	if ( sum( rrna ) == 0 ){
		object@rRNA = rep( 0, ncol( object@x ) )
	} else if ( sum( rrna ) == 1 ){
		object@rRNA = object@x[rrna,]/totals
	} else {
		object@rRNA = colSums( object@x[rrna,] ) / totals
	}

	# Identify protein coding genes
	protein_coding = ( object@molecule_type == 'P' )
	object@protein = colSums( object@x[protein_coding,] ) / totals

	# Exclude plastid genomes
	genome_keep = ( object@genome_type != 'P' ) & ( object@genome_type != 'M' )

	# Identify rows that have at least 2 libraries with count greater than 0
	passes_sampling_criterion = rowSums( object@x > 0 ) > 2

	# Create a vector of rows to keep
	keep = protein_coding & genome_keep & passes_sampling_criterion

	# Subsample matrix and row annotations
	object@x = object@x[ keep, ]
	object@lengths = object@lengths[ keep ]
	object@genome_type = object@genome_type[ keep ]
	object@molecule_type = object@molecule_type[ keep ]
	object@blast_hit = object@blast_hit[ keep ]
	object@fpkm = object@fpkm[ keep, ]
	object@tpm = object@tpm[ keep, ]

	# Prepare EdgeR DGE object
	object@edgeR = edgeR::DGEList( counts=object@x, group=object@treatment )
	object@edgeR = edgeR::calcNormFactors( object@edgeR )

	g = nrow( object@x )
	s = ncol( object@x )

	stopifnot(

		length( object@individual )    == s,
		length( object@treatment )     == s,
		# length( id ) ==  s, # not present in all jsons
		length( object@library_id )    == s,
		length( object@sample_prep )   == s,
		length( object@genome_type )   == g,
		length( object@molecule_type ) == g,
		length( object@blast_hit )     == g,
		length( object@rRNA )          == s,
		length( object@protein )       == s

	)

	object
}


#' Summarize expression libraries
#'
#' @param object An Expression object
#' @return A data frame that summarizes each library with expression data in the
#' Expression object
#' @export
setMethod( "summarize_libraries", signature( object = "Expression" ),
	function( object ) {

		# Some json files do not have an id field, need some logic to accommodate this when constructing hte table
		if ( length( object@id ) < 1 ){
			library_summary = data.frame(
				Species = rep( object@species, length( object@library_id ) ),
				Individual = object@individual,
				Treatment = object@treatment,
				Library = object@library_id,
				Preparation = object@sample_prep,
				rRNA = object@rRNA,
				Protein = object@protein,
				Reads = colSums( object@edgeR$counts )
			)
		}
		else{
			library_summary = data.frame(
				Species = rep( object@species, length( object@library_id ) ),
				Individual = object@individual,
				Treatment = object@treatment,
				Library = object@library_id,
				Preparation = object@sample_prep,
				rRNA = object@rRNA,
				Protein = object@protein,
				Reads = colSums( object@edgeR$counts ),
				Run = as.factor( sapply( object@id, get_run ) ),
				Lane = as.factor( sapply( object@id, get_lane ) )
			)
		}

		# Convert all factor columns to strings
		library_summary %<>%
			dplyr::mutate_if( is.factor, as.character )

		return( library_summary )
	}
)


#' Summarize reference sequences
#'
#' @param object An Expression object
#' @return A data frame that summarizes the reference sequences in the
#' Expression object
#' @export
setMethod( "summarize_reference", signature( object = "Expression" ),
	function( object ) {

		reference_summary = data.frame(
			Species = object@species,
			Genes = nrow( object@x )

		)

		# Convert all factor columns to strings
		reference_summary %<>%
			dplyr::mutate_if( is.factor, as.character )

		return( reference_summary )
	}
)

#' Get the species
#'
#' @param object An Expression object
#' @return Species name
#' @export
setMethod( "species", signature( object = "Expression" ),
	function( object ) {

		return( object@species )
	}
)


#' Creates a DESeqDataSet object
#'
#' @param object An Expression object
#' @param design a formula that explains the project design
#' @return A DESeqDataSet object
#' @export
setMethod( "create_DESeq2", signature( object = "Expression", design = "formula" ),
	function( object, design ) {

		colData = data.frame(
			treatment = object@treatment,
			individual = object@individual,
			row.names = object@library_id
		)

		dds = DESeq2::DESeqDataSetFromMatrix(
				countData = floor( object@x ),
				colData = colData,
				design = design
	        )

		return( dds[ rowSums( DESeq2::counts( dds ) ) > 1, ] )
	}
)



################################################################################
# Functions

#' Get run from Illumina fastq sequence header
#'
#' @param header character strings
#' @return run character strings
#' @export
get_run = function( header ) {
	header = as.character( header )
	run = ''
	if ( substring( header, 1, 1 ) == 'H' ) {
		fields = strsplit( header, "-" )
		fields = unlist( fields )
		run = paste( fields[1], fields[2], fields[3], sep = "-" )
	} else {
	  	fields = strsplit( header, "-" )
		fields = unlist( fields )
		run = paste( fields[1], fields[2], sep = "-" )
    }

	return( run )
}



#' Get lane from Illumina fastq sequence header
#'
#' @param header character strings
#' @return lane character strings
#' @export
get_lane = function( header ) {
	header = as.character( header )

	fields = strsplit( header, "-" )
	fields = unlist( fields )
	lane = paste( fields[ length( fields ) - 1 ] )

	return( lane )
}


#' Converts data frame of node annotations to NHX node labels
#'
#' @param df Data frame
#' return character A vector of node label names
dataframe_to_node_labels = function( df ){
	df$node = NULL # The nodes column is added by treeio, go ahead and remove it
	fields = names( df )
	fields = paste( ":", fields, sep="" )
	fields = paste( fields, "=", sep="" )

	labels = apply(
		df,
		1,
		function( x ) paste( c( "[&&NHX", paste( fields, x, sep="" ), "]" ), collapse='' )
	)

	return( labels )
}


#' Converts a NHX node label string to a list of named values
#'
#' @param label Character string
#' @return list A named list of values
#' @examples
#' nhx_label_to_list( "[&&NHX:Ev=S:S=58:ND=0]" )
#' @export
nhx_label_to_list = function( label ){
	label = sub( "\\[\\&\\&NHX:", "", label )
	label = sub( "\\]", "", label )
	labels = unlist( strsplit( label, ":" ) )
	fields = strsplit( labels, "=" )
	names = unlist( lapply( fields, function( x ) x[ 1 ] ) )
	values = unlist( lapply( fields, function( x ) x[ 2 ] ) )

	named_list = as.list( values )
	names( named_list ) = names
	return( named_list )
}


#' Parses text to a gene tree
#'
#' @param tree_text Text representation of a tree in nhx format with notung or phyldog fields.
#' @return phy The tree, as an ape phylo object
#' @export
parse_gene_tree = function( tree_text ){

	tree_tc = textConnection( tree_text )
	tryCatch(
		tree <- treeio::read.nhx( tree_tc ),
		error = function( c ){
			c$message = paste0( c$message, " (could not parse tree ", tree_text, ")" )
			stop( c )
		}
	)
	close( tree_tc )

	if ( "D" %in% names( tree@data ) ){
		# Notung style annotations, convert speciation annotation to that of phyldog
		colnames( tree@data )[which( names( tree@data ) == "D" )] = "Ev"
		tree@data$Ev[ tree@data$Ev == "Y" ] = "D"
		tree@data$Ev[ tree@data$Ev == "N" ] = "S"
	}
	else if ( "Ev" %in% names( tree@data ) ) {
		# Phyldog style annotations
		# Retype numeric fields
		tree@data$S =  as.numeric( tree@data$S )
		tree@data$ND = as.numeric( tree@data$ND )
	}

	# Parse some NHX fields into tree labels
	Annotations = tree@data

	# Make sure they are ordered by node number
	Annotations =
		Annotations[ order( as.numeric( Annotations$node ), na.last=FALSE ), ]

	# Retain only the annotations for internal nodes
	Annotations =
		Annotations[ ( length( tree@phylo$tip.label ) + 1 ):nrow( Annotations ), ]

	labels = dataframe_to_node_labels( Annotations )

	# Names don't necessarilly reflect node numbers, get rid of them
	# to avoid later confusion
	names( labels ) = NULL

	tree@phylo$node.label = labels


	# Parse node labels for tips and internal nodes from phylo
	phy_node_names = rep( NA, nrow( tree@data ) )
	if ( "node.label" %in% names( tree@phylo ) ){
		phy_node_names = c( tree@phylo$tip.label, tree@phylo$node.label )
	}

	# Parse sequence id, the integer after @, from the tip names
	sequence_ids = phy_node_names
	sequence_ids[ !grepl( '@', sequence_ids ) ] = NA
	sequence_ids = as.numeric( sub( '^.+@', '', sequence_ids, perl=TRUE ) )

	# Parse species, the character string before @, from the tip names
	species_names = phy_node_names
	species_names[ !grepl( '@', species_names ) ] = NA
	species_names = sub( '@.+$', '', species_names, perl=TRUE )
	species_names = sub( '_', ' ', species_names, perl=TRUE )

	# Add node depth
	node_depth = ape::node.depth( tree@phylo )

	tree@data = as.tbl(cbind(
		tree@data,
		phy_node_names=phy_node_names,
		species=species_names,
		sequence_ids=sequence_ids,
		node_depth=node_depth,
		stringsAsFactors=FALSE
		))

	return( tree )
}


#' Parse support from node name
#'
#' @param node_name A character string of format '83:N', where '83' is node support
#' and N indicates whether the node is a duplication or not
#' @return Numeric indicating node support
#' @export
node_support = function( node_name ) {

	ns = NA
	fields = strsplit( node_name, ":" )[[ 1 ]]
	if ( length( fields ) == 2 ){
		ns = as.numeric( fields[ 1 ] )
	}

	return( ns )
}


#' Parse species and sequence id from a phy
#'
#' @param phy The tree, as an ape phylo object. Tip names must have `species@@id` format.
#' @return Dataframe with one row per tip, a species column, and an id column
#' @export
get_tip_info = function( phy ) {

	tips = phy$tip.label

	tip_info = as.data.frame(
		matrix( unlist( strsplit( tips, split='@' ) ), nrow=length( tips ), byrow=T )
	)

	names( tip_info ) = c( "species", "id" )
	tip_info$id = as.numeric( as.character( tip_info$id ) ) # Change from factors to numeric

	# species names have spaces in expression data, remove underscores to make them consistent
	tip_info[,1] = sub( '_', ' ', tip_info[,1] )

	return( tip_info )

}


#' Returns true if all species are in the tree tips
#'
#' @param phy The tree, as an ape phylo object
#' @param species Vector of species names
#' @return logical
#' @export
has_species = function( phy, species ) {

	psp = get_tip_info( phy )[ , 1 ]

	present = species %in% psp

	return( all( present ) )

}


#' Takes a DGEList object, and returns a matrix of normalized counts
#'
#' @param dge edgeR DGEList object, to which normalizations have been applied
#' @return matrix of normalized counts
#' @export
apply_normalizations = function( dge ){

	# Calculate the normalization multipliers, scaled to reads per million
	norm = 1e6 / ( dge$samples[ , "lib.size" ] * dge$samples[ , "norm.factors" ] )

	# Apply them to the counts
	norm_counts = t( t( dge$counts )*norm )

	return( norm_counts )
}


#' Plots a matrix
#'
#' @param m The matrix
#' @param ... additional arguments for image
#' @export
plot_matrix = function( m, ... ) {

	nr = nrow( m )
	nc = ncol( m )
	graphics::image( 1:nc, 1:nr, t( m[ nr:1, ] ), axes=F,xlab="", ylab="", ... )
}


#' Decomposes a gene tree into a list of subtrees that have no duplication events. Assumes notung
#' style node annotations.
#'
#' @param nhx The tree, as a treeio nhx object
#' @return The subtrees as a list of ape::phylo object
#' @export
decompose_orthologs = function( nhx ){

	# Identify duplicated nodes
	duplications = which( nhx@data$Ev == "D" )

	phy = nhx@phylo

	# Get their immediate descendants, which define the clades we want to excise
	to_prune = phy$edge[ , 2 ][ phy$edge[ , 1 ] %in% duplications ]

	subtrees = hutan::decompose_tree( phy, to_prune )

	return( subtrees )

}


#' Create a data frame with summary statistics for expression
#' libraries
#'
#' @param e A list of Expression objects
#' @return A data frame of summary statistics
#' @export
summary_libraries = function( e ){

	library_summary = dplyr::bind_rows( lapply( e, summarize_libraries ) )

	library_summary =
		library_summary[ with( library_summary, order( Species, Individual, Treatment ) ), ]

	library_summary$Species =
		sub( "^(\\w)\\w+", "\\1.", library_summary$Species, perl=TRUE ) # Shorten species names

	return( library_summary )
}

#' Create a data frame with summary statistics for reference sequences
#'
#' @param e A list of Expression objects
#' @return A data frame of summary statistics
#' @export
summary_references = function( e ){
	reference_summary = dplyr::bind_rows( lapply( e, summarize_reference ) )
	reference_summary$Species = sub( "^(\\w)\\w+", "\\1.", reference_summary$Species, perl=TRUE ) # Shorten species names

	return( reference_summary )


}


#' Create a data frame with summary statistics for edges in a phyldog NHX tree
#'
#' @param nhx A treeio nhx object
#' @param default_length_val The default length of branches
#' @return A data frame of edge summary statistics
#' @export
summarize_edges = function ( nhx, default_length_val=NA ) {

	# Create a data frame of internal node annotations
	tags = nhx@data
	tags$node = as.numeric( tags$node )
	tags$S = as.numeric( tags$S )
	tags$ND = as.numeric( tags$ND )
	tags = tags[order( tags$node ),]

	parents = nhx@phylo$edge[,1]
	children = nhx@phylo$edge[,2]

	terminal = rep( FALSE, nrow( nhx@phylo$edge ) )
	terminal[ children <= length( nhx@phylo$tip.label ) ] = TRUE

	df = data.frame(
		gene_tree = rep( digest::digest( nhx ), nrow( nhx@phylo$edge ) ),
		length = nhx@phylo$edge.length,
		Ev_parent = tags$Ev[parents],
		S_parent  = as.numeric( tags$S[parents] ),
		ND_parent = as.numeric( tags$ND[parents] ),
		Ev_child = tags$Ev[children],
		S_child  = as.numeric( tags$S[children] ),
		ND_child = as.numeric( tags$ND[children] ),
		terminal = terminal,
		default_length = FALSE,
		stringsAsFactors = FALSE
	)

	if ( ! is.na( default_length_val ) ){
		df$default_length = dplyr::near( nhx@phylo$edge.length, default_length_val )
	}

	return( df )
}

#' Create a data frame with summary statistics for nodes in a phyldog NHX tree
#'
#' @param nhx A treeio nhx object
#' @param default_length_val The default length of branches, used to determine
#' if a node has descendant branches with default length
#' @return A data frame of node summary statistics
#' @export
summarize_nodes = function ( nhx, default_length_val=NA ) {

	# Create a data frame of internal node annotations
	tags = cbind(
		gene_tree= digest::digest( nhx ),
		nhx@data
	)

	# Add a boolean column that indicates if nodes are parents to
	# edges with default length
	tags$default_length = FALSE
	if ( ! is.na( default_length_val ) ){
		default_edges = dplyr::near( nhx@phylo$edge.length, default_length_val )
		parent_nodes = nhx@phylo$edge[ , 1 ]
		default_nodes = parent_nodes[ default_edges ]
		tags$default_length[ default_nodes ] = TRUE
	}

	tags %<>%
		dplyr::mutate_if( is.factor, as.character )

	return( tags )
}
