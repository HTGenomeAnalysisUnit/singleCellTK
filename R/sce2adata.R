utils_big_as.matrix <- function(
  sparseMat,
  n_slices_init=1,
  verbose=T
  ) {

  n_slices <- n_slices_init-1
  while (TRUE) {
    list_densemat = list()
    n_slices = n_slices+1
    if (verbose) message(paste0("n_slices=",n_slices))
    idx_to = 0
    for (slice in 1:n_slices) {
      if (verbose) message(paste0("converting slice ",slice,"/",n_slices))
      idx_from <- idx_to+1
      idx_to <- if (slice<n_slices) as.integer(ncol(sparseMat)*slice/n_slices) else ncol(sparseMat)
      if (verbose) message(paste0("columns ", idx_from,":", idx_to))
      densemat_sub = try(
        expr = {
          as.matrix(sparseMat[,idx_from:idx_to])
        }, silent = if (verbose) FALSE else TRUE)
      if ("try-error" %in% class(densemat_sub)) {
        break # exit to while loop
      } else {
        list_densemat[[slice]] = densemat_sub
      }
    }
    if (length(list_densemat)==n_slices) break # exit while loop
  }
  if (verbose) message("cbind dense submatrices")
  densemat <- Reduce(f=cbind, x=list_densemat)
  return(densemat)
}

#' Coverts SingleCellExperiment object from R to anndata.AnnData object in
#' Python
#'
#' The AnnData object here can be saved to .h5ad file and read into Python
#' interactive console. Mostly used senario is when you want to apply
#' reticulated Python function, which only works with an anndata.AnnData object.
#' @param SCE A SingleCellExperiment object.
#' @param useAssay Character, default `"counts"`. The name of assay of
#' interests that will be set as the primary matrix of the output AnnData.
#' Available options can be listed by `assayNames(SCE)`. Thee primary matrix
#' will be saved in `adata$X`, Other assays will be stored in `adata$obsm`
#' together with the low-dimension representations (for now).
#' @return A Python anndata.AnnData object
.sce2adata <- function(SCE, useAssay = 'counts') {
    # Transfer SCE object back to AnnData
    # Argument check first
    stopifnot(inherits(SCE, "SingleCellExperiment"))

    # Extract information that correspond to AnnData structure
    #X <- utils_big_as.matrix(t(SummarizedExperiment::assay(SCE, useAssay)))
    X <- t(SummarizedExperiment::assay(SCE, useAssay))
    message("Initialize AnnData using ", useAssay)
    AnnData <- sc$AnnData(X = X, dtype = 'float32')

    message("Copy obs data")
    obs <- as.data.frame(SummarizedExperiment::colData(SCE))
    if(length(obs) > 0){
        AnnData$obs = obs
    } else {
        AnnData$obs_names <- colnames(SCE)
    }

    message("Copy var data")
    var <- as.data.frame(SummarizedExperiment::rowData(SCE))
    if(length(var) > 0){
        AnnData$var = var
    } else {
        AnnData$var_names <- rownames(SCE)
    }
    # uns  <- S4Vectors::metadata(SCE)
    # if(length(uns) > 0){
    #     AnnData$uns <- uns
    # }

    message("Copy obsm data")   
    obsmNames <- SingleCellExperiment::reducedDimNames(SCE)
    if(length(obsmNames) > 0){
        for (i in seq_along(obsmNames)) {
            AnnData$obsm$'__setitem__'(obsmNames[i],
                                       SingleCellExperiment::reducedDim(SCE, obsmNames[i]))
        }
    }

    # Furthermore, the other assays will for now also be saved to .layers
    message("Copy additional layers data")
    allAssayNames <- SummarizedExperiment::assayNames(SCE)
    for (i in seq_along(allAssayNames)) {
        oneName <- allAssayNames[i]
        if (!oneName == useAssay) {
          message("Adding layer ", oneName)
            AnnData$layers$'__setitem__'(oneName,
                                             t(SummarizedExperiment::assay(
                                                 SCE, oneName)))
        }
    }
    return(AnnData)
}
