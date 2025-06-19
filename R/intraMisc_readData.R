# quiets concerns of R CMD check
if (getRversion() >= "2.15.1") utils::globalVariables(c("."))
#' Load Host Single-Cell or Spatial Transcriptomic Data.
#'
#' @param named.outs.dirs Default: String or string vector, representing the paths to the outs directories of the output results from Cell Ranger or Spacer Ranger.
#' @param data.type  Default: \code{"sc"}. Selected from \code{"sc"} and \code{"spt"}, sc represents single cell transcriptome data;
#' spt represents spatial transcriptome sequencing data.
#' @return seurat object
#' @export
#' @examples
#' named.outs.dirs <- c(
#'   system.file("extdata",
#'     "sc/SAMN20336351/outs/",
#'     package = "intraMisc"
#'   ),
#'   system.file("extdata",
#'     "sc/SAMN20336352/outs/",
#'     package = "intraMisc"
#'   )
#' )
#' names(named.outs.dirs) <- c("SAMN20336351", "SAMN20336352")
#' seurat.object <- read_Host(named.outs.dirs)
read_Host <- function(named.outs.dirs, data.type = "sc") {
  if (data.type == "sc") {
    seurat.object.list <- list()
    cat("Reading single cell RNA-seq data, total", length(named.outs.dirs), "samples...\n")
    for (i in seq_along(named.outs.dirs)) {
      cat("---Prosessing the", i, "/", length(named.outs.dirs), "sample:", names(named.outs.dirs)[i], "\n")
      if (file.exists(paste0(
        (named.outs.dirs)[i],
        "/filtered_feature_bc_matrix.h5"
      ))) {
        seurat.object.list[[i]] <- Seurat::Read10X_h5(paste0((named.outs.dirs)[i], "/filtered_feature_bc_matrix.h5")) %>%
          Seurat::CreateSeuratObject(project = names(named.outs.dirs)[i]) %>%
          Seurat::RenameCells(add.cell.id = names(named.outs.dirs)[i])
      } else if (file.exists(paste0(
        (named.outs.dirs)[i],
        "/filtered_feature_bc_matrix/matrix.mtx.gz"
      ))
      ) {
        seurat.object.list[[i]] <- Seurat::Read10X(paste0((named.outs.dirs)[i], "/filtered_feature_bc_matrix/")) %>%
          Seurat::CreateSeuratObject(project = names(named.outs.dirs)[i]) %>%
          Seurat::RenameCells(add.cell.id = names(named.outs.dirs)[i])
      }
    }
    cat("Merging", length(named.outs.dirs), "host seurat objects.\n")
    seurat.object <- merge(seurat.object.list[[1]], seurat.object.list[-1])
  } else if (data.type == "spt") {
    seurat.object <- Seurat::Load10X_Spatial(
      data.dir = named.outs.dirs,
      filename = "filtered_feature_bc_matrix.h5"
    ) %>%
      Seurat::RenameCells(add.cell.id = names(named.outs.dirs))
  }
  return(seurat.object)
}


#' Load Intracellular Microbiome Data.
#'
#' @param named.outs.dirs String or string vector, representing the paths to the outs directories of the output results from Cell Ranger or Spacer Ranger.
#' @param align.method Default: \code{"kraken2"}. intraAlignment alignment methods, selected from ["kraken2", "blast"]
#'
#' @return microbe object
#' @export
#'
#' @importFrom utils read.delim2
#' @importFrom stats na.omit
#' @examples
#' named.outs.dirs <- c(
#'   system.file("extdata",
#'     "sc/SAMN20336351/outs/",
#'     package = "intraMisc"
#'   ),
#'   system.file("extdata",
#'     "sc/SAMN20336352/outs/",
#'     package = "intraMisc"
#'   )
#' )
#' names(named.outs.dirs) <- c("SAMN20336351", "SAMN20336352")
#' microbe.phyloseq.object <- read_Microbe(named.outs.dirs)
read_Microbe <- function(named.outs.dirs,
                         align.method = "kraken2") {
  # named.outs.dirs <- samples.dirs
  microbe.seurat.object.list <- list()
  cat("Reading single cell Microbiome data, total", length(named.outs.dirs), "samples...\n")
  for (i in seq_along(named.outs.dirs)) {
    cat("---Prosessing the", i, "/", length(named.outs.dirs), "sample:", names(named.outs.dirs)[i], "\n")
    microbe.seurat.object.list[[i]] <- read.delim2(
      gzfile(paste0(
        (named.outs.dirs)[i], "/intraAlignment/",
        align.method,
        "/combinedMpa.",
        align.method,
        ".txt.gz"
      )),
      check.names = FALSE, row.names = 1
    ) %>%
      dplyr::filter(rownames(.) %in% grep("Bacteria.*\\|s__",
        rownames(.),
        value = TRUE
      )) %>%
      as.matrix() %>%
      Matrix::Matrix(sparse = TRUE) %>%
      {
        rownames(.) <- rownames(.) %>%
          stringr::str_replace_all("\\|", "#") %>%
          stringr::str_replace_all("_", "-") %>%
          stringr::str_replace_all(" ", "-")
        colnames(.) <- paste0(names(named.outs.dirs)[i], "_", colnames(.))
        .
      } %>%
      Seurat::CreateSeuratObject(min.cells = 0, min.features = 0, assay = "Microbe")
  }

  cat("Merging", length(named.outs.dirs), "microbe seurat objects.\n")
  microbe.seurat.object <- merge(
    microbe.seurat.object.list[[1]],
    microbe.seurat.object.list[-1]
  )
  # The original plan was saved as phyloseq object,
  # while sparse matrix storage method is not used in the otu_table() slot,
  # this will result in phyloseq objects taking up a big amount of
  # memory space when analyzing large queue data.
  # So we temporarily use Seurat object to storage the microbe data.
  return(microbe.seurat.object)
}


#' Load and Integrate Host Transcriptome and Intracellular Microbiome
#'
#' @param named.outs.dirs String or string vector, representing the paths to the outs directories of the output results from Cell Ranger or Spacer Ranger.
#' @param data.type  Default: \code{"sc"}. Selected from \code{"sc"} and \code{"spt"}, sc represents single cell transcriptome data;
#' spt represents spatial transcriptome sequencing data.
#' @param align.method Default: \code{"kraken2"}. intraAlignment alignment methods, selected from ["kraken2", "blast"]
#' @param max.dim Default: \code{30000}.
#' Return phyloseq object if cell numbers <= \strong{\code{max.dim}}, else return seurat object.
#'
#' @return seurat object with microbe data. "seurat.object@misc$microbe"
#' @export
#'
#' @examples
#' named.outs.dirs <- c(
#'   system.file("extdata",
#'     "sc/SAMN20336351/outs/",
#'     package = "intraMisc"
#'   ),
#'   system.file("extdata",
#'     "sc/SAMN20336352/outs/",
#'     package = "intraMisc"
#'   )
#' )
#' names(named.outs.dirs) <- c("SAMN20336351", "SAMN20336352")
#' seurat.object <- read_intraMisc(named.outs.dirs)
#'
read_intraMisc <- function(named.outs.dirs,
                           data.type = "sc",
                           # read_Microbe
                           align.method = "kraken2",
                           max.dim = 30000) {
  seurat.object <- read_Host(named.outs.dirs, data.type = data.type)
  seurat.object@misc$microbe$seurat.object <- read_Microbe(
    named.outs.dirs,
    align.method
  ) %>%
    subset(cells = Seurat::Cells(seurat.object))

  if (dim(seurat.object@misc$microbe$seurat.object)[2] <= max.dim) {
    # Convert microbe.seurat.object to microbe.phyloseq.object
    seurat.object@misc$microbe$phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
      seurat2phyloseq()
  }
  # add nCount_Microbe and nFeature_Microbe to host metadata
  seurat.object@meta.data$nCount_Microbe <- seurat.object@misc$microbe$seurat.object$nCount_Microbe
  seurat.object@meta.data$nFeature_Microbe <- seurat.object@misc$microbe$seurat.object$nFeature_Microbe
  return(seurat.object)
}
