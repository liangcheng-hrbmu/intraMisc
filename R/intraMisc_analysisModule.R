# quiets concerns of R CMD check
if (getRversion() >= "2.15.1") {
  utils::globalVariables(c(
    "nFeature_RNA", "nCount_RNA", "percent.mt",
    "num_cells_expressed", "mean_expression", "dispersion_empirical",
    "dispersion_fit"
  ))
}

###########################################################
#' Basic Processing Workflow for Seurat Objects.
#'
#' Include: filtering cells by the percentage of "mitochondrial genes" (Default: no filter),
#' \code{"NormalizeData > FindVariableFeatures > ScaleData > RunPCA > IntegrateLayers(harmony)"}.
#' Users can also analyze the data themselves according to the standard workflow in \link{Seurat}.
#'
#' @param seurat.object seurat object
#'
#' @param filter.cells Default: \code{FALSE}. Logical.
#' Weather to filtering cells by the percentage of "mitochondrial genes".
#'
#' @param normalization.method Default: \code{"LogNormalize"}. Details see \code{\link[Seurat]{NormalizeData}}
#'
#' @param scale.factor Default: \code{10000}. Details see \code{\link[Seurat]{NormalizeData}}
#'
#' @param selection.method Default: \code{"vst"}. Details see \code{\link[Seurat]{FindVariableFeatures}}
#'
#' @param nfeatures Default: \code{2000}. Details see \code{\link[Seurat]{FindVariableFeatures}}
#'
#' @param integrate.method Default: \code{Seurat::CCAIntegration}.
#' Details see \code{\link[Seurat]{IntegrateLayers}}
#'
#' @param integrate.new.reduction Default: \code{"cca"}. Details see \code{\link[Seurat]{IntegrateLayers}}
#'
#' @param dim.run.umap Default: \code{1:20}. Details see \code{\link[Seurat]{RunUMAP}}
#'
#' @param dim.run.tsne Default: \code{1:20}. Details see \code{\link[Seurat]{RunTSNE}}
#'
#' @param dim.find.neighbors Default: \code{1:20}. Details see \code{\link[Seurat]{FindNeighbors}}
#'
#' @param resolution.find.clusters Default: \code{0.5}. Details see \code{\link[Seurat]{FindClusters}}
#'
#' @param verbose Default: \code{FALSE}. Logical. Controls verbosity
#'
#' @return seurat object
#'
#' @export
#'
#' @examples
#' data(seurat.object)
#' seurat.object <- runSeurat(seurat.object)
#'
runSeurat <- function(seurat.object,
                      # filter cells
                      filter.cells = TRUE,
                      # NormalizeData
                      normalization.method = "LogNormalize",
                      scale.factor = 10000,
                      # FindVariableFeatures
                      selection.method = "vst",
                      nfeatures = 2000,
                      integrate.method = Seurat::CCAIntegration,
                      integrate.new.reduction = "cca",
                      dim.run.umap = 1:20,
                      dim.run.tsne = 1:20,
                      dim.find.neighbors = 1:20,
                      resolution.find.clusters = 0.5,
                      verbose = FALSE) {
  # filtering cells
  if ("RNA" %in% names(seurat.object@assays)) {
    if (filter.cells) {
      # calculated the percent of "mitochondrial genes" in each cells.
      # mitochondrial genes
      seurat.object[["percent.mt"]] <- Seurat::PercentageFeatureSet(seurat.object, pattern = "^MT-")
      seurat.object <- subset(seurat.object,
        subset = nFeature_RNA >= 100 &
          nCount_RNA >= 200 &
          percent.mt <= 10
      )
    }
    if ("data" %in% SeuratObject::DefaultLayer(seurat.object@assays$RNA)) {
      SeuratObject::DefaultLayer(seurat.object@assays$RNA) <- "counts"
    }

    for (layer_name in c("data", "scale.data")) {
      if (layer_name %in% SeuratObject::Layers(seurat.object, search = NA)) {
        SeuratObject::LayerData(seurat.object, layer = layer_name) <- NULL
      }
    }
    if ("counts" %in% SeuratObject::Layers(seurat.object, search = NA)) {
      seurat.object[["RNA"]] <- split(seurat.object[["RNA"]],
        f = seurat.object$orig.ident,
        layers = c("counts")
      )
    }
    seurat.object <- seurat.object %>%
      Seurat::NormalizeData(
        normalization.method = normalization.method,
        scale.factor = scale.factor,
        verbose = verbose
      )

    seurat.object <- seurat.object %>%
      Seurat::FindVariableFeatures(
        selection.method = selection.method,
        nfeatures = nfeatures,
        verbose = verbose
      )

    seurat.object <- seurat.object %>%
      Seurat::ScaleData(verbose = verbose) %>%
      Seurat::RunPCA(verbose = verbose)

    seurat.object <- seurat.object %>%
      Seurat::IntegrateLayers(
        method = integrate.method,
        orig.reduction = "pca",
        k.weight = 50,
        new.reduction = integrate.new.reduction,
        verbose = verbose
      )

    seurat.object <- seurat.object %>%
      Seurat::RunUMAP(
        reduction = integrate.new.reduction,
        umap.method = "uwot",
        metric = "cosine",
        dims = dim.run.umap,
        verbose = verbose
      )

    seurat.object <- seurat.object %>%
      Seurat::RunTSNE(
        reduction = integrate.new.reduction,
        dims = dim.run.tsne,
        verbose = verbose
      )

    seurat.object <- seurat.object %>%
      Seurat::FindNeighbors(
        reduction = integrate.new.reduction,
        dims = dim.find.neighbors,
        verbose = verbose
      )

    seurat.object <- seurat.object %>%
      Seurat::FindClusters(
        resolution = resolution.find.clusters,
        verbose = verbose
      )
    seurat.object <- seurat.object %>% SeuratObject::JoinLayers()
    seurat.object@misc$host$markers <- Seurat::FindAllMarkers(
      seurat.object,
      only.pos = TRUE,
      min.pct = 0.25,
      logfc.threshold = 1,
      verbose = verbose
    )
  } else if ("Spatial" %in% names(seurat.object@assays)) {
    seurat.object <- seurat.object %>%
      Seurat::SCTransform(assay = "Spatial", return.only.var.genes = FALSE, verbose = verbose)

    seurat.object <- seurat.object %>%
      Seurat::RunPCA(assay = "SCT", verbose = verbose)

    seurat.object <- seurat.object %>%
      Seurat::FindNeighbors(reduction = "pca", dims = dim.find.neighbors)

    seurat.object <- seurat.object %>%
      Seurat::FindClusters(resolution = resolution.find.clusters)

    seurat.object <- seurat.object %>%
      Seurat::RunUMAP(reduction = "pca", dims = dim.run.umap)

    seurat.object <- seurat.object %>%
      Seurat::RunTSNE(reduction = "pca", dims = dim.run.tsne)
  }


  return(seurat.object)
}


#' Single-Cell Type Annotation via annoCell (SingleR-based)
#'
#' @param seurat.object seurat object.
#' @param method Default: \code{"SingleR"}. Cell type annotation method.
#' @param clustered Default: \code{FALSE}.
#' User could set it with \code{TRUE} if \code{\link[Seurat]{FindClusters}}
#' function has been done on the seurat object.
#' @param host_organism Default: \code{"human"}.
#' Organism of host single cell RNA-seq data,
#' to determine which reference to use.
#' Selected from \code{["huamn", "mouse"]}.
#' User could also use other reference by using param \code{"custom_ref"}.
#' \describe{
#' \item{\code{"huamn"}}{Using \code{"HumanPrimaryCellAtlasData()"} as reference.
#' Details see \code{\link[celldex]{HumanPrimaryCellAtlasData}}}
#' \item{\code{"mouse"}}{Using \code{"MouseRNAseqData()"} as reference.
#' Details see \code{\link[celldex]{MouseRNAseqData}}}
#' }
#' @param custom_ref Default: \code{NULL}. \cr
#' \strong{- \code{\link{SingleR}} bulid-in references, have been migrated to the \code{celldex} package.}
#' \describe{
#' \item{\code{"HumanPrimaryCellAtlasData()"}}{\code{\link[celldex]{HumanPrimaryCellAtlasData}}.
#' Normalized expression values for 713 microarray samples from the Human Primary Cell Atlas}
#' \item{"MouseRNAseqData()"}{\code{\link[celldex]{MouseRNAseqData}}.
#' Normalized expression values of 358 bulk RNA-seq samples of sorted cell populations that can be found at GEO.}
#' \item{\code{"ImmGenData()"}}{\code{\link[celldex]{ImmGenData}}.
#' Normalized expression values of 830 microarray samples of pure mouse immune cells, generated by the Immunologic Genome Project (ImmGen).}
#' \item{...}{...}
#' }
#' \strong{- Additional references. Details see chapter 3 in \code{\link{SingleR}} package document.}
#'
#' @param ... additional params from function \strong{SingleR}.
#' Details params in \code{\link[SingleR]{SingleR}}.
#'
#' @return seurat object
#' @export
#'
#' @examples
#' # data("seurat.object")
#' # seurat.object <- runSeurat(seurat.object)
#' # seurat.object <- annoCell(seurat.object)
#'
annoCell <- function(seurat.object,
                     method = "SingleR",
                     clustered = FALSE,
                     host_organism = "human",
                     custom_ref = NULL,
                     ...) {
  cat("\033[1;31mManual cell type annotation is recommended.\033[0m\n")
  seurat.object <- seurat.object %>% SeuratObject::JoinLayers()
  if (!"data" %in% names(seurat.object@assays$RNA@layers)) {
    stop('Runing "NormalizeData" function from "Seurat" package or
         "runSeurat" function from this package to get the "data" layer first!')
  }
  if (is.null(custom_ref)) {
    if (host_organism == "human") {
      ref_use <- celldex::HumanPrimaryCellAtlasData()
    } else if (host_organism == "mouse") {

    }
  } else {
    ref_use <- custom_ref
  }
  if (clustered && (!is.null(seurat.object$seurat_clusters))) {
    pred <- SingleR::SingleR(
      test = as.matrix(seurat.object@assays$RNA$data),
      ref = ref_use,
      labels = ref_use$label.main,
      clusters = seurat.object$seurat_clusters
    )
    for (clust in unique(seurat.object$seurat_clusters)) {
      seurat.object@meta.data$singleR_cluster[seurat.object$seurat_clusters == clust] <- pred[clust, "labels"]
    }

    celltype <- data.frame(
      seurat_clusters = rownames(pred),
      singleR_cluster = pred$labels,
      cellNum = table(seurat.object$seurat_clusters)[rownames(pred)] %>% as.vector(),
      stringsAsFactors = FALSE
    )
    print(celltype)
  } else {
    pred <- SingleR::SingleR(
      test = as.matrix(seurat.object@assays$RNA$data),
      ref = ref_use,
      labels = ref_use$label.main,
      ...
    )
    seurat.object@meta.data$singleR <- pred$labels
    print(table(pred$labels) %>%
      as.data.frame() %>%
      {
        colnames(.) <- c("singleR", "cellNum")
        .
      })
  }
  return(seurat.object)
}

#' Pseudotime Analysis of Single-Cell Data via runTrajectory(Monocle3-Based).
#'
#' @param seurat.object seurat object.
#' @param method Default: \code{"monocle3"}.
#'
#' @return seurat object, with cell trajectory result in misc$host$trajectory slot.
#' @export
#' @import DDRTree
#' @examples
#' data("seurat.object")
#' # seurat.object <- runSeurat(seurat.object)
#' # seurat.object <- runTrajectory(seurat.object)
#'
runTrajectory <- function(seurat.object,
                          method = "monocle3") {
  check_package("monocle3")
  check_package("VGAM")
  # get count matrix from seurat object.
  seurat.object <- seurat.object %>% SeuratObject::JoinLayers()
  expression_data <- methods::as(
    as.matrix(SeuratObject::GetAssayData(seurat.object, layer = "counts")), "sparseMatrix"
  )

  # construct gene feature information.
  gene_metadata <- data.frame(
    gene_id = rownames(expression_data), gene_short_name = rownames(expression_data)
  )
  rownames(gene_metadata) <- rownames(expression_data)

  # construct sample(cell) metadata.
  cell_metadata <- seurat.object@meta.data

  cds <- monocle3::new_cell_data_set(
    expression_data = expression_data,
    cell_metadata = cell_metadata,
    gene_metadata = gene_metadata
  )

  ## Step 1: Normalize and pre-process the data
  cds <- monocle3::preprocess_cds(cds, num_dim = 100)

  ## Step 2: Remove batch effects with cell alignment
  cds <- monocle3::align_cds(cds, alignment_group = "orig.ident")

  ## Step 3: Reduce the dimensions using UMAP
  cds <- monocle3::reduce_dimension(cds, reduction_method = "UMAP")

  ## Step 4: Cluster the cells
  cds <- monocle3::cluster_cells(cds)

  ## Step 5: Learn a graph
  cds <- monocle3::learn_graph(cds)

  ## Step 6: Order cells
  cds <- monocle3::order_cells(cds, reduction_method = "UMAP")

  seurat.object@misc$host$trajectory$monocle3 <- cds

  return(seurat.object)
}

###########################################################
# calculate Alpha diversity

#' Alpha Diversity Analysis for Intracellular Microbiome Profiles
#'
#' @param seurat.object seurat object
#' @param measures Default :\code{c("Simpson", "Shannon")}.
#' Meaning that all available alpha-diversity measures will be included.
#' Could selected form \code{c("Observed", "Chao1", "ACE", "Shannon", "Simpson", "InvSimpson", "Fisher")}.
#' Details see \code{\link[phyloseq]{estimate_richness}} from package \code{\link{phyloseq}}.
#' @param tax_level Default: \code{"species"} taxonomy level.
#'
#' @return seurat object
#' @export
#' @examples
#' data(seurat.object)
#' seurat.object <- calcAlpha(seurat.object, tax_level = "genus")
#' head(seurat.object@meta.data)
#'
calcAlpha <- function(seurat.object,
                      measures = c("Simpson", "Shannon"),
                      tax_level = "species") {
  if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
    cat("Converting microbe seurat object to phyloseq object...\n")
    microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
      seurat2phyloseq()
    seurat.object@misc$microbe$phyloseq.object <- microbe.phyloseq.object
  } else {
    microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
  }
  cat("Calculating microbe alpha indexes...\n")
  alpha_indexes <- microbe.phyloseq.object %>%
    phyloseq_select_level(tax_level) %>%
    phyloseq::estimate_richness(measures = measures)

  colnames(alpha_indexes) <- paste0(colnames(alpha_indexes), "_", tax_level)
  seurat.object <- Seurat::AddMetaData(seurat.object, alpha_indexes)

  return(seurat.object)
}


#' Identification of Differential Microbial Features Using LEfSe.
#'
#' @param seurat.object seurat.object
#' @param group.by Default: \code{"seurat_clusters"}.
#' @param compairs Default: \code{c("0", "3")}
#'
#' @return seurat object with LEfSe analysis result in \code{seurat.object@misc$microbe$lefse} solt.
#' @export
#'
#' @examples
#' NULL
runLEfSe <- function(seurat.object,
                     group.by = "seurat_clusters",
                     compairs = c("0", "3")) {
  if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
    cat("Converting microbe seurat object to phyloseq object...\n")
    microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
      seurat2phyloseq()
    seurat.object@misc$microbe$phyloseq.object <- microbe.phyloseq.object
  } else {
    microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
  }

  otu_tab.norm <- microbe.phyloseq.object@otu_table@.Data %>%
    t() %>%
    stats::aggregate(by = list(paste0(seurat.object@meta.data$orig.ident, "-", seurat.object@meta.data[, group.by])), sum) %>%
    {
      rownames(.) <- .[, "Group.1"]
      .[, 2:ncol(.)]
    } %>%
    t() %>%
    as.data.frame() %>%
    # scale(., center = F, scale = colSums(.)) %>%
    as.data.frame() %>%
    {
      .[, !is.na(colSums(.))]
    }

  metaData <- data.frame(
    sample = colnames(otu_tab.norm),
    group.by = stringr::str_split(colnames(otu_tab.norm), "-", n = 2, simplify = TRUE)[, 2],
    row.names = colnames(otu_tab.norm)
  )
  colnames(metaData)[2] <- group.by

  metaData$Group <- metaData[, group.by]
  if (compairs[2] == "others" && (!"others" %in% metaData$Group)) {
    metaData$Group[metaData$Group != compairs[1]] <- "others"
  } else {
    metaData <- metaData[metaData$Group %in% compairs, ]
    otu_tab.norm <- otu_tab.norm[, rownames(metaData)]
  }

  taxoTab <- microbe.phyloseq.object@tax_table %>% as.data.frame()
  taxoTab$kingdom <- paste0("k__", taxoTab$kingdom)
  taxoTab$phylum <- paste0("p__", taxoTab$phylum)
  taxoTab$class <- paste0("c__", taxoTab$class)
  taxoTab$order <- paste0("o__", taxoTab$order)
  taxoTab$family <- paste0("f__", taxoTab$family)
  taxoTab$genus <- paste0("g__", taxoTab$genus)
  taxoTab$species <- paste0("s__", taxoTab$species)
  check_package("microeco")
  dataset <- microeco::microtable$new(
    sample_table = metaData,
    otu_table = otu_tab.norm,
    tax_table = taxoTab
  )

  lefse <- microeco::trans_diff$new(
    dataset = dataset,
    method = "lefse",
    group = "Group",
    alpha = 0.01,
    lefse_subgroup = NULL,
    p_adjust_method = "none",
    taxa_level = "species"
  )
  if ("lefse" %in% ls()) {
    lefse$parameters$group.by <- group.by
    lefse$parameters$compairs <- compairs
    seurat.object@misc$microbe$lefse <- lefse
  }
  return(seurat.object)
}


#' Correlation Analysis Between Microbes and Host Genes
#'
#' @param seurat.object seurat.object
#' @param host.feature Default: \code{"marker_gene"}
#' \describe{
#' \item{\code{"marker_gene"}}{Using maker genes of seurat_clusters.}
#' }
#' @param microbe.feature Default: \code{"all"}
#' \describe{
#' \item{\code{"all"}}{Using all species.}
#' }
#' @param group.by Default: \code{"seurat_clusters"}
#'
#' @return seurat object with correlation analysis result in \code{seurat.object@misc$corr$corr_result} solt.
#' @export
#'
#' @examples
#' NULL
calcCorr <- function(seurat.object,
                     host.feature = "marker",
                     microbe.feature = "all",
                     group.by = "seurat_clusters") {
  if (is.null(seurat.object@misc$microbe$phyloseq.object)) {
    cat("Converting microbe seurat object to phyloseq object...\n")
    microbe.phyloseq.object <- seurat.object@misc$microbe$seurat.object %>%
      seurat2phyloseq()
    seurat.object@misc$microbe$phyloseq.object <- microbe.phyloseq.object
  } else {
    microbe.phyloseq.object <- seurat.object@misc$microbe$phyloseq.object
  }
  if (host.feature == "marker") {
    host_matrix <- seurat.object@assays$RNA$counts[seurat.object@misc$host$markers$gene, ] %>%
      as.data.frame() %>%
      t() %>%
      stats::aggregate(by = list(paste0(seurat.object@meta.data$orig.ident, "-", seurat.object@meta.data[, group.by])), sum) %>%
      {
        rownames(.) <- .[, "Group.1"]
        .[, 2:ncol(.)]
      } %>%
      t() %>%
      as.data.frame() %>%
      scale(., center = FALSE, scale = colSums(.)) %>%
      as.data.frame() %>%
      {
        .[, !is.na(colSums(.))]
      }
  }
  if (microbe.feature == "all") {
    microbe_matrix <- microbe.phyloseq.object@otu_table@.Data %>%
      as.data.frame() %>%
      t() %>%
      stats::aggregate(by = list(paste0(seurat.object@meta.data$orig.ident, "-", seurat.object@meta.data[, group.by])), sum) %>%
      {
        rownames(.) <- .[, "Group.1"]
        .[, 2:ncol(.)]
      } %>%
      t() %>%
      as.data.frame() %>%
      # scale(., center = F, scale = colSums(.)) %>%
      as.data.frame() %>%
      {
        .[, !is.na(colSums(.))]
      }
  }
  host_matrix <- host_matrix[, intersect(colnames(host_matrix), colnames(microbe_matrix))]
  microbe_matrix <- microbe_matrix[, intersect(colnames(host_matrix), colnames(microbe_matrix))]

  A <- t(host_matrix)
  A <- A[, colSums(A) != 0]
  B <- t(microbe_matrix)

  correlation_matrix <- as.data.frame(matrix(NA, ncol = ncol(B), nrow = ncol(A)))
  p_values_matrix <- as.data.frame(matrix(NA, ncol = ncol(B), nrow = ncol(A)))

  for (i in seq_len(ncol(A))) {
    for (j in seq_len(ncol(B))) {
      correlation_test <- stats::cor.test(c(A[, i]), c(B[, j]))
      correlation_matrix[i, j] <- correlation_test$estimate
      p_values_matrix[i, j] <- correlation_test$p.value
    }
  }
  rownames(correlation_matrix) <- colnames(A)
  colnames(correlation_matrix) <- colnames(B)
  rownames(p_values_matrix) <- colnames(A)
  colnames(p_values_matrix) <- colnames(B)

  corr_result <- methods::new("corrResult",
    correlation_matrix = correlation_matrix,
    p_values_matrix = p_values_matrix
  )

  corr_result@parameters$group.by <- group.by
  corr_result@parameters$host.feature <- host.feature
  corr_result@parameters$microbe.feature <- microbe.feature
  seurat.object@misc$corr$corr_result <- corr_result
  print(corr_result)
  return(seurat.object)
}
