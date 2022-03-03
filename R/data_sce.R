#' Load the SCE for the batch-corrected merged data set
#'
#' @description Load the SingleCellExperiment object that holds the batch-corrected
#' reduced dimensionality results (MUT/WT = batch).
#'
#' @usage sce.filt <- load_0606T1filt()
#'
#' @seealso \code{\link{load_0606T1shared}}, \code{\link{load_0606T1merged}}
#' @return an SCE object that needs to be assigned to an object in the environment
#'
#' @export
#'
load_0606T1filt <- function(){
    out <- load_sce(which_assays = "all", sample = "0606T1")

}

#' Load the SCE for the batch-corrected merged data set
#'
#' @description Load the SingleCellExperiment object that holds the batch-corrected
#' reduced dimensionality results (MUT/WT = batch).
#'
#' @usage sce.merged <- load_0606T1merged()
#' @return an SCE object that needs to be assigned to an object in the environment
#'
#' @export
#'
load_0606T1merged <- function(){
    out <- load_sce(which_assays = "reconstructed", sample = "0606T1Merged")
    return(out)
}

#' Load the SCE for the batch-corrected merged data set
#'
#' @description Load the SingleCellExperiment object that holds the SCE
#' representing cells with clonotypes that are present in both conditions.
#' The UMAP coordinates were re-calculated on the reduced subset after removal
#' of suspected doublets (see the vignette about the filtering).
#'
#' @return an SCE object that needs to be assigned to an object in the environment
#'
#' @usage sce.shared <- load_0606T1shared()
#'
#' @export
load_0606T1shared <- function(){
    out <- load_sce(which_assays = "logcounts", sample = "0606T1Shared")
    return(out)
}


#' Filtered cells ' ' @format Named list of cell numbers following the filtering steps
#described in this vignette.
"cell_filt"

#' Filtered genes ' ' @format Named list of gene numbers following the filtering steps
#described in this vignette.
"gene_filt"


#' Load the filtered and processed SingleCellExperiment object
#'
#' @description Use this function to load the processed and filtered gene expression
#' data plus the clonotype information stored within one SingleCellExperiment
#' object.
#'
#'
#' @param which_assays can be "all" (default) or individual assays, e.g. c("logcounts",
#' "counts") etc. If space and memory are problematic, definitely limit the selection here!
#' assays that are available are: "counts", "logcounts"
#' @param ... Additional parameters passed on to \code{load_RDSdata_from_Box}, e.g.
#' \code{check_for_updates = TRUE}
#'
#' @details
#' For the entire code of the filtering and processing, see the vignette
#' \code{01_processing.Rmd}.
#'
#' The resulting SCE object contains the usual content: colData with information
#' about individidual cells, rowData with info about individual genes, reducedDims,
#' etc.
#'
#' The \code{colData} includes:
#'
#' \describe{
#' \item{Barcode:}{Each cell's barcode used for keeping track of its identity during sequencing.}
#' \item{Sample:}{'WT' or 'MUT'}
#' \item{raw_clonotype_id:}{e.g. 'clonotype94'}
#' \item{cdr3s_aa:}{The amino acid sequence of the CDR3 portion, e.g. "TRA:CIARGGGGADGLTF;TRA:CGADRNGNEKLTF;TRB:CASSLTTDREPYEQYF"}
#' \item{multiTRA:}{TRUE/FALSE entries based on whether \code{cdr3s_aa} contained more than one entry for TRA}
#' \item{multiTRB:}{TRUE/FALSE entries based on whether \code{cdr3s_aa} contained more than one entry for TRB}
#' \item{numTRA:}{Number of TRA sequences within \code{cdr3s_aa}}
#' \item{numTRB:}{Number of TRB sequences within \code{cdr3s_aa}}
#' \item{cluster:}{clustering results of all cells}
#' }
#'
#' The object also containes the coordinates from dimensionality reductions (see
#' examples for more details).
#'
#'
#' @usage sce.filt <- load_sce(which_assays = "logcounts", sample = "21LT2")
#'
#' @return A SingleCellExperiment object with cells from 'mutant' samples
#'  (stimulation with tumor antigen) and from the 'wt' sample (stimulation
#'  with an irrelevant antigen)).
#'
#'
#' @examples \dontrun{
#'
#' > library(SingleCellExperiment)
#' > sce.21LT2 <- load_sce(which_assays = "all", sample = "21LT2")
#'
#' > reducedDimNames(sce.21LT2)
#' "corrected" "TSNE"      "UMAP"
#'
#' > assayNames(sce.21LT2)
#' [1] "counts"                "logcounts"
#' }
#'
#' @return SCE object
#'
#'
#' @export
#'
load_sce <- function(which_assays = "all", sample = "Sample", ...){

    ## the Box links are noted in the text file
    fl <- system.file("extdata", paste0("sce_storage_", sample, ".txt"),
        package = "Klebanoff0606T1")
    if(fl == ""){stop(paste("sce_storage_", sample, ".txt does not exist in package 'Klebanoff0606T1'."))}

    inf <- read.table(fl,stringsAsFactors = FALSE)

    if(unique(inf$V3) != sample){stop("The sce_storage.txt file must contain a third column holding the sample name. Which should be the same as the one specified via sample = .")}

    ## DOWNLOAD AND CACHE THE FILES FROM THE BOX ==============================
    ## note: using the default cache of BioC here, we may want to change that
    ## to something more specific via the `cache_path` option of `load_RDSdata_fromBox()`

    ## load colData
    cold <- load_RDSdata_from_Box(
        shared_link = inf[inf$V1 == "colData",]$V2, data_name = paste0("KlebColData",sample), ...)

    ## load rowData
    rowd <- load_RDSdata_from_Box(
        shared_link = inf[inf$V1 == "rowData",]$V2, data_name = paste0("KlebRowData", sample) , ...)

    ## get reducedDims
    rdms <- load_RDSdata_from_Box(
        shared_link = inf[inf$V1 == "reducedDims",]$V2, data_name = paste0("KlebRedDims", sample), ... )

    ## metadata
    metd <- load_RDSdata_from_Box(
        shared_link = inf[inf$V1 == "metadata",]$V2, data_name = paste0("KlebMetadata", sample), ... )

    ## get assayData
    if(which_assays == "all"){
        ## extract corresponding assay entry from the text file
        asss <- grep("^assay:", unique(inf$V1), value = TRUE)
    }else{
        asss <- unlist(lapply(which_assays, function(x) grep(paste0(":",x,"$"), unique(inf$V1), ignore.case = TRUE, value=TRUE)))
        if(length(which_assays) == 0){
            warning("None of the assays you specified are part of the file stored in inst/extdata, i.e. we can't find the links.")
        }
    }

    assl <- list()
    for(i in asss){
        j <- gsub("^assay:","", i)
        assl[[j]] <-  load_RDSdata_from_Box(
            shared_link = inf[inf$V1 == i,]$V2, data_name = paste0("Kleb",j, sample), ...)
    }

    ## construct the SCE object =============================================
    return(SingleCellExperiment::SingleCellExperiment(assays = assl,
        colData = cold, rowData = rowd,
        metadata = metd, reducedDims = rdms))
}

